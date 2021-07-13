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
	mat3d F(Fe[0], Fe[3], Fe[6], Fe[1], Fe[4], Fe[7], Fe[2], Fe[5], Fe[8]);

	// radial (from file)
	vec3d f_rad(fl[0], fl[1], fl[2]);

	// axial (constant)
	vec3d f_axi(0.0, 0.0, 1.0);

	// circumferential (from cross product)
	vec3d f_cir = f_rad ^ f_axi;

	// determinant of the deformation gradient
	const double J = F.det();

	// get current and end times
	const double t = *time;
	const double endtime = 11.0;							// 11.0 | 31.0-32.0 (TEVG)

	const double eps = std::numeric_limits<double>::epsilon();
	const double partialtime = endtime;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	const double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)

	const double imper = 0.00;					// imper > 0 for TORTUOSITY (see Matlab script <NodesElementsAsy.m>) | 0.00 | 20.0
	const double rIo = 0.6468;					// 0.6468 | 0.5678
	const double hwaves = 2.0;
	const double lo = 30.0;
	const double ro = sqrt(f_rad*f_rad);

	// retrieve local element basis directions
	vec3d N[3];

	// pointwise, consistent with mesh generated with Matlab script <NodesElementsAsy.m>
	N[0] = f_rad;
	N[1] = f_cir;
	N[2] = f_axi;

	const double phieo = 0.34;								// 0.34 (CMAME | KNOCKOUTS) | 1.00 (TEVG) | 1.0/3.0 (TEVG)
	const double phimo = 0.5*(1.0-phieo);
	const double phico = 0.5*(1.0-phieo);

	const double eta = 1.0;									// 1.0 | 1.0/3.0 (for uniform cases) | 0.714

	const double mu = 89.71;
	const double Get = 1.90;
	const double Gez = 1.62;

	double alpha = 0.522;								// original orientation of diagonal collagen | 0.522 (CMAME | KNOCKOUTS) | 0.8713 (TEVG)

	// original homeostatic parameters (adaptive)

	// passive
	const double cm = 261.4;									// 261.4 (CMAME | KNOCKOUTS) | 46.61 (TEVG)
	const double dm = 0.24;
	const double Gm = 1.20;
	const double cc = 234.9;									// 234.9 (CMAME | KNOCKOUTS) | 328.475 (TEVG)
	const double dc = 4.08;
	const double Gc = 1.25;

	// orientation fractions for collagen
	const double betat = 0.056;
	const double betaz = 0.067;
	const double betad = 0.5*(1.0 - betat - betaz);

	const vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);		// original diagonal fiber direction
	const vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);		// idem for symmetric

	// active
	const double Tmax = 250.0 * 0.0;							// 250.0 | 50.0 | 150.0 (for uniform cases, except for contractility -> 250)
	const double lamM = 1.1;
	const double lam0 = 0.4;
	const double CB = sqrt(log(2.0));							// such that (1-exp(-CB^2)) = 0.5
	const double CS = 0.5*CB * 1.0;							// such that (1-exp( -C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0

	const double KsKi = 0.35;
	const double EPS  = 1.0+(1.0-1.0)*(sgr-1.0)/(endtime-1.0);

	const double KfKi   = 1.0;
	const double inflam = 0.0*(sgr-1.0)/(endtime-1.0);

	const double aexp = 1.0;									// 1.0 (KNOCKOUTS | TEVG) | 0.0 (CMAME | TORTUOSITY)

	const double delta = 0.0;

	// compute U from polar decomposition of deformation gradient tensor
	mat3ds U; mat3d R; F.right_polar(R,U);

	// right Cauchy-Green tensor and its inverse
//	mat3ds C = pt.RightCauchyGreen();
	const mat3ds C = (F.transpose() * F).sym();
	const mat3ds Ci = C.inverse();

	// Ge from spectral decomposition
	const mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);

	// define identity tensor and some useful dyadic products of the identity tensor
	const mat3dd  I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);

	// stress for elastin
	const mat3ds Se = (phieo*mu*Ge*Ge).sym();						// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge

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

		const double Jdep = 0.9999;
		const double lm = 1.0e3*mu;

		const double lt = (F*N[1]).norm();
		const double lz = (F*N[2]).norm();
		const double lp = (F*Np).norm();
		const double ln = (F*Nn).norm();

		const double lmt2 = (Gm*lt)*(Gm*lt);
		const double lct2 = (Gc*lt)*(Gc*lt);
		const double lcz2 = (Gc*lz)*(Gc*lz);
		const double lcp2 = (Gc*lp)*(Gc*lp);
		const double lcn2 = (Gc*ln)*(Gc*ln);

		const mat3ds tent = dyad(F*N[1]);
		const mat3ds tenz = dyad(F*N[2]);
		const mat3ds tenp = dyad(F*Np);
		const mat3ds tenn = dyad(F*Nn);

		// passive
		const mat3ds Sm = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		const mat3ds Sc =	(cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );
		const tens4ds cf = phimo*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
				phico*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
						2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
						2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
						2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);

		// active
		const mat3ds Sa = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lt*lt)*dyad(N[1]);
		const mat3ds Sx = Se + phimo * Sm + phico * Sc + phimo * Sa;
		const tens4ds ca = phimo*2.0*Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*dyad1s(tent);

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
//
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

		const double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);	// phim from <J*phim/phimo=(J*phic/phico)^eta>

		// recompute remodeled original stresses for smc and collagen (from remodeled natural configurations)

		const double lto = (Fio.inverse()*N[1]).norm();
		const double lzo = (Fio.inverse()*N[2]).norm();
		const double lpo = (Fio.inverse()*Np).norm();					// original referential stretch for deposition stretch calculation
		const double lno = (Fio.inverse()*Nn).norm();					// idem for symmetric

		const double lmt2 = (Gm*lto)*(Gm*lto);
		const double lct2 = (Gc*lto)*(Gc*lto);
		const double lcz2 = (Gc*lzo)*(Gc*lzo);
		const double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		const double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		const double lr = (F*(Fio*N[0])).norm();						// lr -> 1 for F -> Fo
		const double lt = (F*(Fio*N[1])).norm();						// lt -> 1 for F -> Fo
		const double lz = (F*(Fio*N[2])).norm();						// lz -> 1 for F -> Fo

		alpha = atan(tan(alpha)*pow(lt/lz,aexp));				// update alpha
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);					// update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);					// idem for symmetric

		// passive
		const mat3ds Smo = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		const mat3ds Sco = (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );

		// active
		const mat3ds Sao = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lto*lto)*dyad(N[1]);

		mat3ds Uo; mat3d Ro; (Fio.inverse()).right_polar(Ro,Uo);	// Uo from polar decomposition
		const mat3d  uo(Uo);

		smo = 1.0/Jo*(uo*(Smo*uo)).sym();
		sco = 1.0/Jo*(uo*(Sco*uo)).sym();

		const mat3ds sao = 1.0/Jo*(uo*(Sao*uo)).sym();

		// compute current stresses

		const double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;				// rIrIo -> rIorIo = 1 for F -> Fo

		const mat3ds sNm = phim*smo;									// phim*smhato = phim*smo
		const mat3ds sNc = phic*sco;									// phic*schato = phic*sco

		const mat3ds sNf = sNm + sNc;

		const double Cratio = CB-CS*(EPS*pow(rIrIo,-3)-1.0);
		mat3ds sNa; sNa.zero();
		if (Cratio>0) sNa = phim*(1.0-exp(-Cratio*Cratio))/(1.0-exp(-CB*CB))*sao;

		const mat3ds Ui = U.inverse();            					// inverse of U
		const mat3d  ui(Ui);

		const mat3ds Sf = J*(ui*sNf*ui).sym();						// J*Ui*sNf*Ui
		const mat3ds Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

		const mat3ds Sx = Se + Sf + Sa;

		const double p = 1.0/3.0/J*Sx.dotdot(C) - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p

		S = Sx - J*p*Ci;

//		pt.m_Jh  = J;
//		pt.m_Fih = F.inverse();

		// 2nd P-K stresses
		const mat3ds Se = (phieo*mu*Ge*Ge).sym();								// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge
		const mat3ds Sm = J*(ui*sNm*ui).sym();						// J*Ui*sNm*Ui
		const mat3ds Sc = J*(ui*sNc*ui).sym();						// J*Ui*sNc*Ui
//		mat3ds Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

//		mat3ds Sf = phim*Sm+phic*Sc;							// J*Ui*sNf*Ui

		// mrp089: is this consistent with the above definition from 2PK?
//		mat3ds Sx = Se+Sf+phim*Sa;

		// associated Cauchy stresses
		const mat3ds sm = 1.0/J*(F*(Sm*F.transpose())).sym();
		const mat3ds sc = 1.0/J*(F*(Sc*F.transpose())).sym();
		const mat3ds sa = 1.0/J*(F*(Sa*F.transpose())).sym();
		const mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();

		// compute tangent

		tens4dss IxIss = tens4dss(IxI);							// IxI in tens4dss form
		tens4dss IoIss = tens4dss(IoI);							// IoI in tens4dss form

		const tens4dss Ixsx = dyad1ss(I,sx);
		const tens4dss smxI = dyad1ss(sm,I);
		const tens4dss scxI = dyad1ss(sc,I);
		const tens4dss saxI = dyad1ss(sa,I);

		const mat3ds tenr = dyad(F*(Fio*N[0]));						// Fio needed for consistency (from computation of lr)
		const mat3ds tent = dyad(F*(Fio*N[1]));
		const mat3ds tenz = dyad(F*(Fio*N[2]));

		const tens4dss Ixnrr = dyad1ss(I,tenr);
		const tens4dss Ixntt = dyad1ss(I,tent);

		// contribution due to constant Cauchy stresses at constituent level
		const tens4dss cfss(0.0);

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

			const mat3ds ten1 = dyad(Fxeigenvec[i]);

			for (int j=0; j<3; j++) {

				double component = sfpro(i,j) / pow(eigenval[i],3) / eigenval[j];

				const mat3ds ten2 = dyads(Fxeigenvec[i],Fxeigenvec[j]);

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

		const double dphiRm = phimo*eta*pow(J/Jo*phic/phico,eta-1.0)/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);
		const double dphiRc = phico/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);

		cfss += dphiRm*(smxI+saxI) + dphiRc*scxI;

		// contribution due to the ratio of vasocontrictors to vasodilators in the active stress
		const tens4dss saoxnrr = dyad1ss((R*sao*R.transpose()).sym(),tenr);
		const tens4dss saoxntt = dyad1ss((R*sao*R.transpose()).sym(),tent);

		// 1/J * FoF : [ J * phim * 1/(1.0-exp(-CB*CB)) * (Ui*sao*Ui) x d(1-exp(-Cratio^2))/d(C/2) ] : (Ft)o(Ft)
		const tens4dss cass = phim * 6.0*Cratio*CS*EPS*pow(rIrIo,-4)*exp(-Cratio*Cratio)/(1.0-exp(-CB*CB)) * (ro/rIo/lt*saoxntt-(ro-rIo)/rIo/lr*saoxnrr);

		// contribution due to change in Cauchy stresses at constituent level (orientation only, for now)
		const tens4dss cpnss(0.0);

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

		const double scphato = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);	// constant stress magnitude at constituent level
		const double scnhato = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);

		const vec3d dNpdta = (N[1]-N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		const vec3d dNndta = (N[1]+N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);

		const mat3ds ten1 = 1.0/Jo*dyads(R*(Uo*dNpdta),R*(Uo*Np));					// FoF : (Ui)o(Ui) : d(NpxNp)/d(tan(alpha)), with Jo and Uo needed for consistency (from computation of sco)
		const mat3ds ten2 = 1.0/Jo*dyads(R*(Uo*dNndta),R*(Uo*Nn));

		const mat3ds ten3 = aexp*tan(alpha)*(1.0/(lt*lt)*tent-1.0/(lz*lz)*tenz);		// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)

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

	std::cout<<scientific<<std::setprecision(2);

//	std::cout<<"t "<<*time<<std::endl;
//
//	std::cout<<"gpx";
//	for (int i=0; i<3; i++)
//		std::cout<<" "<<xgp[i];
//	std::cout<<std::endl;

//	std::cout<<"f1 "<<f1.x<<" "<<f1.y<<" "<<f1.z<<std::endl;
//	std::cout<<"f2 "<<f2.x<<" "<<f2.y<<" "<<f2.z<<std::endl;
//	std::cout<<"f3 "<<f3.x<<" "<<f3.y<<" "<<f3.z<<std::endl;
//
//	std::cout<<"J "<<J<<std::endl;
//
//	std::cout<<"Fe";
//	for (int i=0; i<9; i++)
//	{
//		if (i%3 == 0)
//			std::cout<<std::endl;
//		std::cout<<"\t"<<Fe[i];
//	}
//	std::cout<<std::endl;

//	std::cout<<"N";
//	for (int i=0; i<3; i++)
//	{
//		std::cout<<"\t"<<N[i].x;
//		std::cout<<"\t"<<N[i].y;
//		std::cout<<"\t"<<N[i].z;
//		std::cout<<std::endl;
//	}
//	std::cout<<std::endl;


//	std::cout<<"S";
//	for (int i=0; i<9; i++)
//	{
//		if (i%3 == 0)
//			std::cout<<std::endl;
//		std::cout<<"\t"<<S_out[i];
//	}
//	std::cout<<"\n"<<std::endl;
//	//
//	std::cout<<"C";
//	for (int i=0; i < 3; i++)
//		for (int j=0; j < 3; j++)
//			for (int k=0; k < 3; k++)
//				for (int l=0; l < 3; l++)
//					std::cout<<" "<<css(i, j, k, l);
//	std::cout<<"\n"<<std::endl;

//	check for nans
	for (int i=0; i<9; i++)
		if (std::isnan(S_out[i]))
			std::terminate();
}
}


void stress_tangent_stvk_(const double* Fe, const double* fl, const double* xgp, const double* time, double* S_out, double* CC_out)
{
	// convert deformation gradient to FEBio format
	mat3d F(Fe[0], Fe[3], Fe[6], Fe[1], Fe[4], Fe[7], Fe[2], Fe[5], Fe[8]);

	// material parameters
	const double young = 1000.0;
	const double nu = 0.4;

	// lame parameters
	const double mu = young / ( 2.0*(1.0+nu) );
	const double lambda = nu*young / ((1.0+nu)*(1.0-2.0*nu));

	// define identity tensor and some useful dyadic products of the identity tensor
	const mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4dss cIxI = tens4dss(IxI);
	tens4dss cIoI = tens4dss(IoI);

	// green-lagrange strain
	mat3ds E = 0.5 * ((F.transpose() * F).sym()- I);

	// stress
	mat3ds S = lambda*E.tr()*I + 2.0 * mu * E;

	// tangent
	tens4dss css = lambda * cIxI + 2.0 * mu * cIoI;

	// convert to vector for FORTRAN
	typedef double (*ten2)[3];
	typedef double (*ten4)[3][3][3];

	ten2 S2 = (ten2) S_out;
	ten4 C4 = (ten4) CC_out;

	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			S2[j][i] = S(i, j);

	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			for (int k=0; k < 3; k++)
				for (int l=0; l < 3; l++)
					C4[l][k][j][i] = css(i, j, k, l);
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
