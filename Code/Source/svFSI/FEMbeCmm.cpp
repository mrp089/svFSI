#include "stdafx.h"
//#include "FEMbeCmm.h"
//#include "FECore/FEAnalysis.h"					// to get end time
//#include "FECore/FEModel.h"						// to get current time
#include <iostream>								// to use cin.get()
#include <iomanip>
#include <limits>
#define _USE_MATH_DEFINES						// to introduce pi constant (1/2)
#include <math.h>								// to introduce pi constant (2/2)
//#include "FECore/log.h"							// to print to log file and/or screen
#include <mat3d.h>
#include <vec2d.h>
#include <vec3d.h>
#include <tens3d.h>
#include <tens4d.h>
#include <stdafx.h>

//-----------------------------------------------------------------------------
// This function needs to return the spatial (i.e. Cauchy stress) at the material point
// which is passed as a parameter. The FEMaterialPoint class contains all the state 
// variables of the material (and its base classes).
extern"C"
{
void stress_tangent_(const double* Fe, const double* fl, const double* xgp, const double* time, double* S_out, double* CC_out)
{
	// convert deformation gradient to FEBio format
	mat3d F(Fe[0], Fe[1], Fe[2], Fe[3], Fe[4], Fe[5], Fe[6], Fe[7], Fe[8]);

	// fiber directions (currently unused)
	vec3d f1(fl[0], fl[1], fl[2]);
	vec3d f2(fl[3], fl[4], fl[5]);
	vec3d f3 = f1 ^ f2;

	// determinant of the deformation gradient
	const double J = F.det();

	// retrieve material position
	vec3d  X = {xgp[0], xgp[1], xgp[2]};
//	vec3d  X = {1.0, 1.0, 0.0};

	// get current and end times
//	double t = GetFEModel()->GetTime().currentTime;
//	double endtime = GetFEModel()->GetCurrentStep()->m_tend;
	const double t = *time;
	const double endtime = 11.0;							// 11.0 | 31.0-32.0 (TEVG)

//	std::cout<<"time "<<*time<<std::endl;
//	std::cout<<"xgp "<<xgp[0]<<" "<<xgp[1]<<" "<<xgp[2]<<std::endl;

	double eps = std::numeric_limits<double>::epsilon();
	double partialtime = endtime;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)

	double imper = 0.00;					// imper > 0 for TORTUOSITY (see Matlab script <NodesElementsAsy.m>) | 0.00 | 20.0
	double rIo = 0.6468;					// 0.6468 | 0.5678
	double hwaves = 2.0;
	double lo = 30.0;
	vec3d  Xcl = {0.0, imper/100.0*rIo*sin(hwaves*M_PI*X.z/lo), X.z};		// center line

	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};								// radial vector

	double ro = sqrt(NX*NX);

	NX /= ro;

	// retrieve local element basis directions
	vec3d N[3];

	// pointwise, consistent with mesh generated with Matlab script <NodesElementsAsy.m>
	N[2] = {0.0, imper/100.0*rIo*hwaves*M_PI/lo*cos(hwaves*M_PI*X.z/lo), 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);		// axial = d(Xcl)/d(z)
	N[1] = {-NX.y, NX.x, NX.z};																					// circumferential
	N[0] = N[2]^N[1];

	// elementwise, from input file
	// N[2] = pt.m_Q.col(0); N[1] = pt.m_Q.col(1); N[0] = pt.m_Q.col(2);							// axial, circumferential, radial

	double phieo = 0.34;								// 0.34 (CMAME | KNOCKOUTS) | 1.00 (TEVG) | 1.0/3.0 (TEVG)
	double phimo = 0.5*(1.0-phieo);
	double phico = 0.5*(1.0-phieo);

	double eta = 1.0;									// 1.0 | 1.0/3.0 (for uniform cases) | 0.714

	double mu = 89.71;
	double Get = 1.90;
	double Gez = 1.62;

	double alpha = 0.522;								// original orientation of diagonal collagen | 0.522 (CMAME | KNOCKOUTS) | 0.8713 (TEVG)

	// original homeostatic parameters (adaptive)

	// passive
	double cm = 261.4;									// 261.4 (CMAME | KNOCKOUTS) | 46.61 (TEVG)
	double dm = 0.24;
	double Gm = 1.20;
	double cc = 234.9;									// 234.9 (CMAME | KNOCKOUTS) | 328.475 (TEVG)
	double dc = 4.08;
	double Gc = 1.25;

	// orientation fractions for collagen
	double betat = 0.056;
	double betaz = 0.067;
	double betad = 0.5*(1.0 - betat - betaz);

	vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);		// original diagonal fiber direction
	vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);		// idem for symmetric

	// active
	double Tmax = 250.0 * 0.0;							// 250.0 | 50.0 | 150.0 (for uniform cases, except for contractility -> 250)
	double lamM = 1.1;
	double lam0 = 0.4;
	double CB = sqrt(log(2.0));							// such that (1-exp(-CB^2)) = 0.5
	double CS = 0.5*CB * 1.0;							// such that (1-exp( -C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0

	double KsKi = 0.35;
	double EPS  = 1.0+(1.0-1.0)*(sgr-1.0)/(endtime-1.0);

	double KfKi   = 1.0;
	double inflam = 0.0*(sgr-1.0)/(endtime-1.0);

	double aexp = 1.0;									// 1.0 (KNOCKOUTS | TEVG) | 0.0 (CMAME | TORTUOSITY)

	double delta = 0.0;

	// compute U from polar decomposition of deformation gradient tensor
	mat3ds U; mat3d R; F.right_polar(R,U);

	// right Cauchy-Green tensor and its inverse
//	mat3ds C = pt.RightCauchyGreen();
	mat3ds C = (F.transpose() * F).sym();
	mat3ds Ci = C.inverse();

	// Ge from spectral decomposition
	mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);

	// define identity tensor and some useful dyadic products of the identity tensor
	mat3dd  I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);

	// stress for elastin
	mat3ds Se = (phieo*mu*Ge*Ge).sym();						// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge

	// spatial moduli for elastin
	tens4ds ce(0.0);								// phieo/J*(FcF:GecGe:Cehat:GecGe:FTcFT) = phieo/J*(FcF:GecGe:0:GecGe:FTcFT)

	double    Jo;
	double   svo;
	mat3ds   smo;
	mat3ds   sco;
	mat3d    Fio;
	double  phic;

	// computation of the second Piola-Kirchhoff stress
	mat3ds S;
	// computation of spatial moduli
	tens4dss css;
	mat3ds sfpro;
	double eigenval[3]; vec3d eigenvec[3];
	if (t <= 1.0 + eps) {

		double Jdep = 0.9999;
		double lm = 1.0e3*mu;

		double lt = (F*N[1]).norm();
		double lz = (F*N[2]).norm();
		double lp = (F*Np).norm();
		double ln = (F*Nn).norm();

		double lmt2 = (Gm*lt)*(Gm*lt);
		double lct2 = (Gc*lt)*(Gc*lt);
		double lcz2 = (Gc*lz)*(Gc*lz);
		double lcp2 = (Gc*lp)*(Gc*lp);
		double lcn2 = (Gc*ln)*(Gc*ln);

		mat3ds tent = dyad(F*N[1]);
		mat3ds tenz = dyad(F*N[2]);
		mat3ds tenp = dyad(F*Np);
		mat3ds tenn = dyad(F*Nn);

		// passive
		mat3ds Sm = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		mat3ds Sc =	(cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );
		tens4ds cf = phimo*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
				phico*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
						2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
						2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
						2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);

		// active
		mat3ds Sa = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lt*lt)*dyad(N[1]);
		mat3ds Sx = Se + phimo * Sm + phico * Sc + phimo * Sa;
		tens4ds ca = phimo*2.0*Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*dyad1s(tent);

		S = Sx + Ci*lm*log(Jdep*J);

		mat3d u(U);

		Jo    = J;
		svo   = 1.0/3.0/J*S.dotdot(C);
		smo   = 1.0/J*(u*(Sm*u)).sym();
		sco   = 1.0/J*(u*(Sc*u)).sym();
		Fio   = F.inverse();
//		Jh    = Jo;
//		Fih   = Fio;
		phic  = phico;

		cf /= J; ca /= J;

		tens4ds c = ce + cf + ca;

		c += lm/J*(IxI-2.0*log(Jdep*J)*IoI);

		css = tens4dss(c);		// c in tens4dss form
	}

	// mrp089: how do we store step I?
	else if (t <= partialtime + eps) {

		phic = phico;																// initial guess
		double dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));		// initial tangent d(R)/d(phic)
		double Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;			// initial residue
		do {																		// local iterations to obtain phic
			phic = phic-Rphi/dRdc;													// phic
			dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));			// tangent
			Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;				// update residue
		} while (abs(Rphi) > sqrt(eps));											// && abs(Rphi/Rphi0) > sqrt(eps) && j<10
		phic = phic-Rphi/dRdc;														// converge phase -> phic (updated in material point memory)

		double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);	// phim from <J*phim/phimo=(J*phic/phico)^eta>

		// recompute remodeled original stresses for smc and collagen (from remodeled natural configurations)

		double lto = (Fio.inverse()*N[1]).norm();
		double lzo = (Fio.inverse()*N[2]).norm();
		double lpo = (Fio.inverse()*Np).norm();					// original referential stretch for deposition stretch calculation
		double lno = (Fio.inverse()*Nn).norm();					// idem for symmetric

		double lmt2 = (Gm*lto)*(Gm*lto);
		double lct2 = (Gc*lto)*(Gc*lto);
		double lcz2 = (Gc*lzo)*(Gc*lzo);
		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		double lr = (F*(Fio*N[0])).norm();						// lr -> 1 for F -> Fo
		double lt = (F*(Fio*N[1])).norm();						// lt -> 1 for F -> Fo
		double lz = (F*(Fio*N[2])).norm();						// lz -> 1 for F -> Fo

		alpha = atan(tan(alpha)*pow(lt/lz,aexp));				// update alpha
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);					// update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);					// idem for symmetric

		// passive
		mat3ds Smo = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		mat3ds Sco = (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );

		// active
		mat3ds Sao = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lto*lto)*dyad(N[1]);

		mat3ds Uo; mat3d Ro; (Fio.inverse()).right_polar(Ro,Uo);	// Uo from polar decomposition
		mat3d  uo(Uo);

		smo = 1.0/Jo*(uo*(Smo*uo)).sym();
		sco = 1.0/Jo*(uo*(Sco*uo)).sym();

		mat3ds sao = 1.0/Jo*(uo*(Sao*uo)).sym();

		// compute current stresses

		double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;				// rIrIo -> rIorIo = 1 for F -> Fo

		mat3ds sNm = phim*smo;									// phim*smhato = phim*smo
		mat3ds sNc = phic*sco;									// phic*schato = phic*sco

		mat3ds sNf = sNm + sNc;

		double Cratio = CB-CS*(EPS*pow(rIrIo,-3)-1.0);
		mat3ds sNa; sNa.zero();
		if (Cratio>0) sNa = phim*(1.0-exp(-Cratio*Cratio))/(1.0-exp(-CB*CB))*sao;

		mat3ds Ui = U.inverse();            					// inverse of U
		mat3d  ui(Ui);

		mat3ds Sf = J*(ui*sNf*ui).sym();						// J*Ui*sNf*Ui
		mat3ds Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

		mat3ds Sx = Se + Sf + Sa;

		double p = 1.0/3.0/J*Sx.dotdot(C) - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p

		S = Sx - J*p*Ci;

//		pt.m_Jh  = J;
//		pt.m_Fih = F.inverse();

		// 2nd P-K stresses
		mat3ds Se = (phieo*mu*Ge*Ge).sym();								// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge
		mat3ds Sm = J*(ui*sNm*ui).sym();						// J*Ui*sNm*Ui
		mat3ds Sc = J*(ui*sNc*ui).sym();						// J*Ui*sNc*Ui
//		mat3ds Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

//		mat3ds Sf = phim*Sm+phic*Sc;							// J*Ui*sNf*Ui

		// mrp089: is this consistent with the above definition from 2PK?
//		mat3ds Sx = Se+Sf+phim*Sa;

		// associated Cauchy stresses
		mat3ds sm = 1.0/J*(F*(Sm*F.transpose())).sym();
		mat3ds sc = 1.0/J*(F*(Sc*F.transpose())).sym();
		mat3ds sa = 1.0/J*(F*(Sa*F.transpose())).sym();
		mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();

		// compute tangent

		tens4dss IxIss = tens4dss(IxI);							// IxI in tens4dss form
		tens4dss IoIss = tens4dss(IoI);							// IoI in tens4dss form

		tens4dss Ixsx = dyad1ss(I,sx);
		tens4dss smxI = dyad1ss(sm,I);
		tens4dss scxI = dyad1ss(sc,I);
		tens4dss saxI = dyad1ss(sa,I);

		mat3ds tenr = dyad(F*(Fio*N[0]));						// Fio needed for consistency (from computation of lr)
		mat3ds tent = dyad(F*(Fio*N[1]));
		mat3ds tenz = dyad(F*(Fio*N[2]));

		tens4dss Ixnrr = dyad1ss(I,tenr);
		tens4dss Ixntt = dyad1ss(I,tent);

		// contribution due to constant Cauchy stresses at constituent level
		tens4dss cfss(0.0);

		sfpro.zero();
		sfpro(0,0) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[0]);
		sfpro(1,1) = eigenvec[1]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[1]);
		sfpro(2,2) = eigenvec[2]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);
		sfpro(0,1) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[1]);
		sfpro(1,2) = eigenvec[1]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);
		sfpro(0,2) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);

		vec3d Fxeigenvec[3];

		Fxeigenvec[0] = F*eigenvec[0];
		Fxeigenvec[1] = F*eigenvec[1];
		Fxeigenvec[2] = F*eigenvec[2];

		for (int i=0; i<3; i++) {

			mat3ds ten1 = dyad(Fxeigenvec[i]);

			for (int j=0; j<3; j++) {

				double component = sfpro(i,j) / pow(eigenval[i],3) / eigenval[j];

				mat3ds ten2 = dyads(Fxeigenvec[i],Fxeigenvec[j]);

				cfss -= component*dyad1ss(ten2,ten1);

				for (int k=0; k<3; k++) {

					if (k == i) continue;

					mat3ds ten3 = dyads(Fxeigenvec[j],Fxeigenvec[k]);
					mat3ds ten4 = dyads(Fxeigenvec[k],Fxeigenvec[i]);

					component = sfpro(i,j) / eigenval[i] / eigenval[j] / eigenval[k] / (eigenval[i] + eigenval[k]);

					cfss -= component*dyad1ss(ten3,ten4);
				}
			}
		}

		double dphiRm = phimo*eta*pow(J/Jo*phic/phico,eta-1.0)/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);
		double dphiRc = phico/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);

		cfss += dphiRm*(smxI+saxI) + dphiRc*scxI;

		// contribution due to the ratio of vasocontrictors to vasodilators in the active stress
		tens4dss saoxnrr = dyad1ss((R*sao*R.transpose()).sym(),tenr);
		tens4dss saoxntt = dyad1ss((R*sao*R.transpose()).sym(),tent);

		// 1/J * FoF : [ J * phim * 1/(1.0-exp(-CB*CB)) * (Ui*sao*Ui) x d(1-exp(-Cratio^2))/d(C/2) ] : (Ft)o(Ft)
		tens4dss cass = phim * 6.0*Cratio*CS*EPS*pow(rIrIo,-4)*exp(-Cratio*Cratio)/(1.0-exp(-CB*CB)) * (ro/rIo/lt*saoxntt-(ro-rIo)/rIo/lr*saoxnrr);

		// contribution due to change in Cauchy stresses at constituent level (orientation only, for now)
		tens4dss cpnss(0.0);

//		double lpo = (Fio.inverse()*Np).norm();					// original referential stretch for deposition stretch calculation
//		double lno = (Fio.inverse()*Nn).norm();					// idem for symmetric
//
//		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
//		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric
//
//		double lz = (F*(Fio*N[2])).norm();						// lz -> 1 for F -> Fo

		alpha = atan(tan(alpha)*pow(lt/lz,aexp));				// update alpha
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);					// update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);					// idem for symmetric

		double scphato = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);	// constant stress magnitude at constituent level
		double scnhato = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);

		vec3d dNpdta = (N[1]-N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		vec3d dNndta = (N[1]+N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);

		mat3ds ten1 = 1.0/Jo*dyads(R*(Uo*dNpdta),R*(Uo*Np));					// FoF : (Ui)o(Ui) : d(NpxNp)/d(tan(alpha)), with Jo and Uo needed for consistency (from computation of sco)
		mat3ds ten2 = 1.0/Jo*dyads(R*(Uo*dNndta),R*(Uo*Nn));

		mat3ds ten3 = aexp*tan(alpha)*(1.0/(lt*lt)*tent-1.0/(lz*lz)*tenz);		// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)

		cpnss += (phic*betad) * scphato * dyad1ss(ten1,ten3);					// 1/J * FoF : [ J * phicp * scphato * (Ui)o(Ui) : 2*d(NpxNp)/d(C) ] : (Ft)o(Ft)
		cpnss += (phic*betad) * scnhato * dyad1ss(ten2,ten3);

		tens4dss cess = tens4dss(ce);							// ce in tens4dss form

		css = cess + cfss + cass + cpnss;

		css += 1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddotss(IxIss,css))
									 + svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam)*(IxIss-2.0*IoIss)
									 - 3.0*svo/(1.0-delta)*KsKi*EPS*pow(rIrIo,-4)*(ro/rIo/lt*Ixntt-(ro-rIo)/rIo/lr*Ixnrr);

	}

	// convert to vector for FORTRAN
	typedef double (*ten2)[3];
	typedef double (*ten4)[3][3][3];

	ten2 S2 = (ten2) S_out;
	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			S2[j][i] = S(i, j);

	ten4 C4 = (ten4) CC_out;
	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			for (int k=0; k < 3; k++)
				for (int l=0; l < 3; l++)
					C4[l][k][j][i] = css(i, j, k, l);

// check for nans
	for (int i=0; i<9; i++)
	{
		if (std::isnan(S_out[i]))
		{
			std::cout<<scientific<<std::setprecision(2);

			std::cout<<"t "<<*time<<std::endl;

			std::cout<<"gpx";
			for (int i=0; i<3; i++)
				std::cout<<" "<<xgp[i];
			std::cout<<std::endl;

			std::cout<<"J "<<J<<std::endl;

			std::cout<<"Fe";
			for (int i=0; i<9; i++)
			{
				if (i%3 == 0)
					std::cout<<std::endl;
				std::cout<<"\t"<<Fe[i];
			}
			std::cout<<std::endl;

			std::cout<<"S";
			for (int i=0; i<9; i++)
			{
				if (i%3 == 0)
					std::cout<<std::endl;
				std::cout<<"\t"<<S_out[i];
			}
			std::cout<<std::endl;

			std::cout<<"C";
			for (int i=0; i < 3; i++)
				for (int j=0; j < 3; j++)
					for (int k=0; k < 3; k++)
						for (int l=0; l < 3; l++)
							std::cout<<" "<<css(i, j, k, l);
			std::cout<<std::endl;

			std::terminate();
		}
	}

}
}

//void print_mat(const std::string name, const double* mat) {
//	std::cout<<name;
//	for (int i=0; i<9; i++)
//	{
//		if (i%3 == 0)
//			std::cout<<std::endl;
//		std::cout<<"\t"<<std::showpos<<std::to_string(mat[i]);
//	}
//	std::cout<<std::endl;
//
//}
