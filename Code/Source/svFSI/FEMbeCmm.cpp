#include "stdafx.h"
//#include "FEMbeCmm.h"
//#include "FECore/FEAnalysis.h"					// to get end time
//#include "FECore/FEModel.h"						// to get current time
#include <iostream>								// to use cin.get()
#include <iomanip>
#include <limits>
#define _USE_MATH_DEFINES						// to introduce pi constant (1/2)
#include <math.h>								// to introduce pi constant (2/2)
#include <cmath>
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

void stress_tangent_(const double* Fe, const double* fl, const double* time, double* eVWP, double* grInt, double* S_out, double* CC_out, int* gp)
{
	// convert deformation gradient to FEBio format
	mat3d F(Fe[0], Fe[3], Fe[6], Fe[1], Fe[4], Fe[7], Fe[2], Fe[5], Fe[8]);

	// right Cauchy-Green tensor and its inverse
	const mat3ds C = (F.transpose() * F).sym();
	const mat3ds Ci = C.inverse();

	// compute U from polar decomposition of deformation gradient tensor
	mat3ds U; mat3d R; F.right_polar(R,U);
	double eigenval[3]; vec3d eigenvec[3];
	U.eigen2(eigenval,eigenvec);
	mat3d u(U);

	// determinant of the deformation gradient
	const double J = F.det();
	if (J < 0.0)
		std::terminate();

	// couple wss
	const bool coup_wss = false;

	// set example
	enum Example { none, aneurysm, tortuosity, stenosis };
	Example example = aneurysm;
	const bool example_asym = true;
	
	double KsKi = 0.35;
	// double KsKi = 0.0;
	// double KsKi = 1.0;

	const double curve = 0.0;

	// length multiplier
	const double mult = 1.0;

	// get current time
	double t;
	if (coup_wss)
		// time from partitioned coupling (passed through input file)
		t = eVWP[7];
	else
		// time from monolithic coupling (svFSI provided)
		t = *time;

	// time step size
	const double dt = 1.0;

	// number of time steps between updates
	const int n_t_update = 1;

	// number of time steps for prestressing
	const int n_t_pre = 1;

	// number of time steps total
	const int n_t_end = 101;

	const double pretime = n_t_pre * dt;
	const double endtime = n_t_end * dt;							// 11.0 | 31.0-32.0 (TEVG)

	const double eps = std::numeric_limits<double>::epsilon();
	const double epst = 1.0e-12;
	const double partialtime = endtime;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	const double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)

//	const double imper = 0.00;					// imper > 0 for TORTUOSITY (see Matlab script <NodesElementsAsy.m>) | 0.00 | 20.0
	const double rIo = 0.6468;					// 0.6468 | 0.5678
	const double hwaves = 2.0;
	const double lo = 15.0 * mult;

	const vec3d  X(eVWP[3], eVWP[4], eVWP[5]);
//	const vec3d  Xcl(0.0, imper/100.0*rIo*sin(hwaves*M_PI*X.z/lo), X.z);		// center line
	const vec3d  Xcl(curve / 2.0 * (1.0 - cos(2.0 * M_PI * X.z / lo)), 0.0, X.z);		// center line
	vec3d NX(X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z);								// radial vector

	const double ro = sqrt(NX*NX);

	NX /= ro;

	// retrieve local element basis directions
	vec3d N[3];

	// WSS
	const double tau = eVWP[6];

	// WSS of previous time step
	const double tau_old = eVWP[12];

	// dWSS
	const vec3d dtau(eVWP[9],eVWP[10], eVWP[11]);
//	std::cout<<eVWP[9]<<" "<<eVWP[10]<<" "<<eVWP[11]<<std::endl;

	// pointwise, consistent with mesh generated with Matlab script <NodesElementsAsy.m>
//	const vec3d   n2(0.0, imper/100.0*rIo*hwaves*M_PI/lo*cos(hwaves*M_PI*X.z/lo), 1.0);
	const vec3d n2(curve * M_PI / lo * sin(2.0 * M_PI * X.z / lo), 0.0, 1.0);
	const vec3d n1(-NX.y, NX.x, NX.z);
	N[2] = n2;
	N[2] = N[2]/sqrt(N[2]*N[2]);		// axial = d(Xcl)/d(z)
	N[1] = n1;							// circumferential
	N[0] = N[2]^N[1];					// radial
	
	const double azimuth = acos(-NX.y);							// azimuth wrt axis -Y

	double phieo = 0.34;								// 0.34 (CMAME | KNOCKOUTS) | 1.00 (TEVG) | 1.0/3.0 (TEVG)
	double phimo = 0.5*(1.0-phieo);
	double phico = 0.5*(1.0-phieo);
	double J_star = 1.0;

	const double eta = 1.0;									// 1.0 | 1.0/3.0 (for uniform cases) | 0.714

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
	const double betat = 0.056;
	const double betaz = 0.067;
	const double betad = 0.5*(1.0 - betat - betaz);

	vec3d Np = N[1]*sin(alpha)+N[2]*cos(alpha);		// original diagonal fiber direction
	vec3d Nn = N[1]*sin(alpha)-N[2]*cos(alpha);		// idem for symmetric

	// active
	const double Tmax = 250.0 * 0.0;							// 250.0 | 50.0 | 150.0 (for uniform cases, except for contractility -> 250)
	const double lamM = 1.1;
	const double lam0 = 0.4;
	double CB = sqrt(log(2.0));							// such that (1-exp(-CB^2)) = 0.5
	double CS = 0.5*CB * 1.0;							// such that (1-exp( -C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0

	// time factor [0, 1]
	double f_time = (sgr - pretime) / (endtime - pretime);

	const double EPS  = 1.0+(1.0-1.0)*f_time;

	double KfKi   = 0.0;
	const double inflam = 0.0*f_time;

	double aexp = 0.0;									// 1.0 (KNOCKOUTS | TEVG) | 0.0 (CMAME | TORTUOSITY)

	const double delta = 0.0;

	// hyperelastic incompressibility
	const double Jdep = 0.999;
	const double lm = 1.0e3*mu;
	const double kappa = 1.0e0 * mu;

	// Lagrange multiplier
	double p;
	double p0;

	// WSS ratio
	double tau_ratio;

	// intramural stress
	double svh;

	// stimulus residual
	double ups = 0.0;

	// todo: do outside and once per time step
	// select G&R mode
	enum GR_Mode { prestress, gr, elastic };
	GR_Mode mode;

	// prestressing
	if (0.0 <= t && t <= pretime + epst)
		mode = prestress;
	// last step of prestressing: do one gr update to store history values
	else
		mode = gr;

	// examples from fig. 8, doi.org/10.1016/j.cma.2020.113156
	if (example == aneurysm and mode == gr) {
		// apply transfer function to advance time more equally
		const double t_fac = 2.0;
		f_time = tanh(t_fac * f_time) / tanh(t_fac);

		// no fiber reorientation
		aexp = 0.0;

		// location of aneurysm (= middle)
		const double z_om = lo/2.0;

		double theta_od;
		double phi_e_hm = 0.65;
		const int vza = 5;
		int vzc;
		const double z_od = lo/3.0/mult;
		double nc;
		if (example_asym)
		{
			// 8d
//			z_od = lo/3.0/mult;
//			vz = 5;
			theta_od = M_PI/3.0;

			// z_od = lo/4.0/mult;
			// theta_od = 0.9;
			// vza = 4;
			vzc = 5;
			// nc = 5.0;
		}
		else
		{
			// 8b, 8c
			// z_od = lo/4.0/mult;
			// vza = 2;
		}

		// axial factor (0, 1]
		const double f_axi = exp(-pow(abs((X.z - z_om) / z_od), vza));

		// azimuth factor (0, 1]
		double f_cir = 1.0;
		if (example_asym)
			f_cir = exp(-pow(abs((azimuth - M_PI) / (theta_od)), vzc));
			// f_cir = (nc + cos(azimuth - M_PI)) / (1.0 + nc);
			// f_cir = (nc + (1.0 + cos(pow(azimuth / M_PI - 1.0, vzc) * M_PI)) / 2.0) / (1.0 + nc);

		mu   *= 1.0 - f_time * f_axi * f_cir * phi_e_hm;
		KsKi *= 1.0 - f_time * f_axi;
	}

	if (example == tortuosity){
		// set constants to example
		phieo = 0.331;
		phimo = 0.265;
		phico = 0.404;

		mu = 88.46;									// 88.46

		alpha = 0.665;								// 0.615 (110) / 0.665 (120)

		cm = 419.1;									// 419.1
		dm = 0.105;									// 0.105
		Gm = 1.198;									// 1.198
		cc = 1278.;									// 1278.
		dc = 11.10;									// 11.10
		Gc = 1.092;									// 1.092

		KsKi = 0.0;
		KfKi = 0.0;
		aexp = 0.0;									// 1.0 (KNOCKOUTS) | 0.0 (CMAME and TORTUOSITY)

		// change parameters over time
		if (mode == gr){
			const double muuni = 0.850*mu;						// 0.825 (110) / 0.850 (120)
			const double alpout = 1.00*alpha;
			const double alpin  = 1.32*alpha;						// 1.55 (110) / 1.32 (120)

			mu = mu + (muuni - mu) * f_time;
			alpha = alpout + (alpin - alpout) * f_time;
		}
	}
	if (example == stenosis and mode == gr)
	{
		// location (= middle)
		const double z_om = lo/2.0;

		const double phi_e_hm = 2.0;
		const int vza = 2;
		const double z_od = lo/16.0/mult;

		// axial factor (0, 1]
		const double f_axi = exp(-pow(abs((X.z - z_om) / z_od), vza));

		// azimuth factor (0, 1]
		double f_cir = 1.0;
		if (example_asym)
		{
			double theta_od = 0.75;
			const int vzc = 6;
			f_cir = exp(-pow(abs((azimuth - M_PI) / (M_PI * theta_od)), vzc));
		}

		mu   *= 1.0 + f_time * f_axi * f_cir * phi_e_hm;
		KsKi *= 1.0 - f_time * f_axi * f_cir * phi_e_hm;
	}

	// Ge from spectral decomposition
	const mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);

	// stress for elastin
	const mat3ds Se = (phieo*mu*Ge*Ge).sym();						// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge

	// passive stresses
	mat3ds Sm;
	mat3ds Sc;

	// active stress
	mat3ds Sa;

	// define identity tensor and some useful dyadic products of the identity tensor
	const mat3dd  I(1.0);
	const tens4ds IxI = dyad1s(I);
	const tens4ds IoI = dyad4s(I);
	const tens4dmm IxIss = tens4dmm(IxI);							// IxI in tens4dmm form
	const tens4dmm IoIss = tens4dmm(IoI);							// IoI in tens4dmm form

	// spatial moduli for elastin
	tens4ds ce(0.0);								// phieo/J*(FcF:GecGe:Cehat:GecGe:FTcFT) = phieo/J*(FcF:GecGe:0:GecGe:FTcFT)

	// computation of the second Piola-Kirchhoff stress
	mat3ds S(0.0);

	// computation of spatial moduli
	tens4dmm css;
	mat3ds sfpro;

	// retrieve internal variables 2*(1+1+1+1+6+6+9) = 50
	double   Jo;
	double   Jh;
	double  svo;
	// double  svh;
	double phic = phico;
	double phich;
	double tauo;
	double tauh;
	double po;
	mat3ds  smo;
	mat3ds  smh;
	mat3ds  sco;
	mat3ds  sch;
	mat3d   Fio;
	mat3d   Fih;

	// retrieve prestress
	if (mode == gr or mode == elastic)
	{
		Jo   = grInt[0];
		svo  = grInt[1];
		phic = grInt[2];
		tauo = grInt[3];
		po   = grInt[4];
		int k = 5;
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				smo(i,j) = grInt[k];
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				sco(i,j) = grInt[k];
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
			{
				Fio(i,j) = grInt[k];
				k++;
			}
	}

	// retrieve frozen g&r
	if (mode == elastic)
	{
		Jh    = grInt[25];
		svh   = grInt[26];
		phich = grInt[27];
		tauh  = grInt[28];
		int k = 29;
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				smh(i,j) = grInt[k];
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				sch(i,j) = grInt[k];
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
			{
				Fih(i,j) = grInt[k];
				k++;
			}
	}

	// select material evaluation according to G&R mode
	switch(mode)
	{
	case prestress:
	{
		// compute stress
		const double lt = (F*N[1]).norm();
		const double lz = (F*N[2]).norm();
		const double lp = (F*Np).norm();
		const double ln = (F*Nn).norm();

		const double lmt2 = (Gm*lt)*(Gm*lt);
		const double lct2 = (Gc*lt)*(Gc*lt);
		const double lcz2 = (Gc*lz)*(Gc*lz);
		const double lcp2 = (Gc*lp)*(Gc*lp);
		const double lcn2 = (Gc*ln)*(Gc*ln);

		// passive
		const mat3ds Sm = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		const mat3ds Sc =	(cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );

		// active
		const mat3ds Sa = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lt*lt)*dyad(N[1]);

		const mat3ds Sx = Se + phimo * Sm + phico * Sc + phimo * Sa;

		// Lagrange multiplier during prestress
		po = -lm*log(Jdep*J);

		S = Sx + Ci*lm*log(Jdep*J);

		// compute tangent
		const mat3ds tent = dyad(F*N[1]);
		const mat3ds tenz = dyad(F*N[2]);
		const mat3ds tenp = dyad(F*Np);
		const mat3ds tenn = dyad(F*Nn);

		// passive
		tens4ds cf = phimo*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
				phico*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
						2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
						2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
						2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);

		// active
		tens4ds ca = phimo*2.0*Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*dyad1s(tent);

		cf /= J; ca /= J;

		tens4ds c = ce + cf + ca;

		c += lm/J*(IxI-2.0*log(Jdep*J)*IoI);

		css = tens4dmm(c);		// c in tens4dmm form
		break;
	}
	case gr:
	{
		// compute stress
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

		// assume eta = 1
		const double phic0 = (1.0 - Jo/J * phieo) / (1.0 + phimo / phico);
		const double phim0 = (1.0 - Jo/J * phieo) / (1.0 + phico / phimo);

		// Jacobian from mass fractions
		J_star = Jo * phieo / (1.0 - phic - phim);

		// unneccessary linearizations of J
		// const double Ja =  Jo * phieo;
		// const double Jb =  1.0 / (1.0 + phimo / phico) + 1.0 / (1.0 + phico / phimo);
		// const double J_star0 = Ja / (1.0 - (1.0 - Ja/J) * Jb);
		// const double dJ_star_dJ = pow(J_star0, 2) * Jb / pow(J, 2);

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
		mat3d  uo(Uo);

		smo = 1.0/Jo*(uo*(Smo*uo)).sym();
		sco = 1.0/Jo*(uo*(Sco*uo)).sym();

		const mat3ds sao = 1.0/Jo*(uo*(Sao*uo)).sym();

		// compute current stresses

		const double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;				// rIrIo -> rIorIo = 1 for F -> Fo

		if (coup_wss)
			tau_ratio = tau/tauo;
		else
			tau_ratio = pow(rIrIo,-3);

		mat3ds sNm = phim*smo;									// phim*smhato = phim*smo
		mat3ds sNc = phic*sco;									// phic*schato = phic*sco

		const mat3ds sNf = sNm + sNc;

//		const double Cratio = CB-CS*(EPS*pow(rIrIo,-3)-1.0);
		const double Cratio = CB-CS*(EPS*tau_ratio-1.0);
		mat3ds sNa; sNa.zero();
		if (Cratio>0) sNa = phim*(1.0-exp(-Cratio*Cratio))/(1.0-exp(-CB*CB))*sao;

		const mat3ds Ui = U.inverse();            					// inverse of U
		const mat3d  ui(Ui);

		const mat3ds Sf = J*(ui*sNf*ui).sym();						// J*Ui*sNf*Ui
		const mat3ds Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

		const mat3ds Sx = Se + Sf + Sa;

		// compute current stresses
		sNm = smo;										// phim*smhato = phim*smo
		sNc = sco;										// phic*schato = phic*sco

		// 2nd P-K stresses
		const mat3ds Sm = J*(ui*sNm*ui).sym();						// J*Ui*sNm*Ui
		const mat3ds Sc = J*(ui*sNc*ui).sym();						// J*Ui*sNc*Ui

		// associated Cauchy stresses
		const mat3ds sm = 1.0/J*(F*(Sm*F.transpose())).sym();
		const mat3ds sc = 1.0/J*(F*(Sc*F.transpose())).sym();
		const mat3ds sa = 1.0/J*(F*(Sa*F.transpose())).sym();
		const mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();

		const tens4dmm Ixsx = dyad1mm(I,sx);
		const tens4dmm smxI = dyad1mm(sm,I);
		const tens4dmm scxI = dyad1mm(sc,I);
		const tens4dmm saxI = dyad1mm(sa,I);

		const mat3ds tenr = dyad(F*(Fio*N[0]));						// Fio needed for consistency (from computation of lr)
		const mat3ds tent = dyad(F*(Fio*N[1]));
		const mat3ds tenz = dyad(F*(Fio*N[2]));

		const tens4dmm Ixnrr = dyad1mm(I,tenr);
		const tens4dmm Ixntt = dyad1mm(I,tent);

		// contribution due to constant Cauchy stresses at constituent level
		tens4dmm cfss(0.0);

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

				cfss -= component*dyad1mm(ten2,ten1);

				for (int k=0; k<3; k++) {

					if (k == i) continue;

					mat3ds ten3 = dyads(Fxeigenvec[j],Fxeigenvec[k]);
					mat3ds ten4 = dyads(Fxeigenvec[k],Fxeigenvec[i]);

					component = sfpro(i,j) / eigenval[i] / eigenval[j] / eigenval[k] / (eigenval[i] + eigenval[k]);

					cfss -= component*dyad1mm(ten3,ten4);
				}
			}
		}

		const double dphiRm = phimo*eta*pow(J/Jo*phic/phico,eta-1.0)/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);
		const double dphiRc = phico/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);

		cfss += dphiRm*(smxI+saxI) + dphiRc*scxI;

		// contribution due to the ratio of vasocontrictors to vasodilators in the active stress
		const tens4dmm saoxnrr = dyad1mm((R*sao*R.transpose()).sym(),tenr);
		const tens4dmm saoxntt = dyad1mm((R*sao*R.transpose()).sym(),tent);

		// contribution due to change in Cauchy stresses at constituent level (orientation only, for now)
		tens4dmm cpnss(0.0);

		const double scphato = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);	// constant stress magnitude at constituent level
		const double scnhato = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);

		const vec3d dNpdta = (N[1]-N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		const vec3d dNndta = (N[1]+N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);

		const mat3ds ten1 = 1.0/Jo*dyads(R*(Uo*dNpdta),R*(Uo*Np));					// FoF : (Ui)o(Ui) : d(NpxNp)/d(tan(alpha)), with Jo and Uo needed for consistency (from computation of sco)
		const mat3ds ten2 = 1.0/Jo*dyads(R*(Uo*dNndta),R*(Uo*Nn));

		const mat3ds ten3 = aexp*tan(alpha)*(1.0/(lt*lt)*tent-1.0/(lz*lz)*tenz);		// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)

		cpnss += (phic*betad) * scphato * dyad1mm(ten1,ten3);					// 1/J * FoF : [ J * phicp * scphato * (Ui)o(Ui) : 2*d(NpxNp)/d(C) ] : (Ft)o(Ft)
		cpnss += (phic*betad) * scnhato * dyad1mm(ten2,ten3);

		const tens4dmm cess = tens4dmm(ce);							// ce in tens4dmm form

//		const double p = 1.0/3.0/J*Sx.dotdot(C) - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p
		css = cess + cfss + cpnss;
		tens4dmm css_ref = J*css.pp(F.inverse());

		svh = 1.0/3.0/J*Sx.dotdot(C);

		const double delta_sig = svh / svo - 1.0;
		const double delta_tau = tau_ratio - 1.0;
		ups = delta_sig - KsKi * delta_tau;
		p0 = kappa*ups;

		const mat3ds d_svh = 1.0/3.0/J/svo * (Sx + css_ref.dot(C) - Sx.dotdot(C)/2.0 * Ci);
		// const mat3ds d_svh = (1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddot(IxIss,css))).dot(C)/svo;
		const mat3ds d_tau = -3.0/2.0*pow(rIrIo,-4) * (ro/rIo/lt*tent - (ro-rIo)/rIo/lr*tenr);
		const mat3ds Sp = 2.0*kappa*ups * (d_svh - KsKi * d_tau);

//		const double p = 1.0/3.0/J*Sx.dotdot(C) - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p
		p = svh - svo/(1.0-delta)*(1.0+KsKi*(EPS*tau_ratio-1.0)-KfKi*inflam);		// Ups = 1 -> p
		// p = po;
		// p = po - kappa * (J - J_star) * (1 - dJ_star_dJ); // second part is always 0
		// p = po + svh;

		S = Sx - J*p*Ci;
		// S = Sx + Sp;
		// S = Sx * svo/svh;
		// S = Sx - kappa*ups * J*Ci;
		// S = Sx + Ci*lm*log(Jdep*(ups + 1.0));

		// const double a = 1.0e-2;
		// const double a = 0.0;
		// S = Sx * (1+a) - (svh * (1+a) - svo * (1.0+KsKi*(tau_ratio-1.0))) * J*Ci;

		css += 1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddot(IxIss,css));
		css += svo/(1.0-delta)*(1.0+KsKi*(EPS*tau_ratio-1.0)-KfKi*inflam)*(IxIss-2.0*IoIss);
	
		// wss linearization
		if (!coup_wss)
		{
			css += 3.0 * pow(rIrIo,-4) * phim * 2.0 * Cratio * CS * EPS * exp(-Cratio * Cratio) / (1.0-exp(-CB*CB)) * (ro/rIo/lt*saoxntt-(ro-rIo)/rIo/lr*saoxnrr);
			css -= 3.0 * pow(rIrIo,-4) * svo/(1.0-delta) * KsKi * EPS                                               * (ro/rIo/lt*Ixntt  -(ro-rIo)/rIo/lr*Ixnrr);
		}

		break;
	}
	case elastic:
	{
		// compute stress
		double phimh = phimo/(Jh/Jo)*pow(Jh/Jo*phich/phico,eta);	// evolved homeostatic phimh from <Jh*phimh/phimo=(Jh*phich/phico)^eta>

		double lrh = (Fih.inverse()*(Fio*N[0])).norm();				// lrh -> 1 for Fh -> Fo
		double lth = (Fih.inverse()*(Fio*N[1])).norm();				// lth -> 1 for Fh -> Fo
		double lzh = (Fih.inverse()*(Fio*N[2])).norm();				// lzh -> 1 for Fh -> Fo

		alpha = atan(tan(alpha)*pow(lth/lzh,aexp));					// remodeled (evolved homeostatic) alpha at h (alphah)
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);						// update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);						// idem for symmetric

		mat3ds Uo; mat3d Ro; (Fio.inverse()).right_polar(Ro,Uo);	// Uo from polar decomposition
		mat3ds Uh; mat3d Rh; (Fih.inverse()).right_polar(Rh,Uh);	// Uh from polar decomposition

		double lt = (F*(Uh.inverse()*(Uo*N[1]))).norm();
		double lz = (F*(Uh.inverse()*(Uo*N[2]))).norm();
		double lp = (F*(Uh.inverse()*(Uo*Np))).norm();
		double ln = (F*(Uh.inverse()*(Uo*Nn))).norm();

		double lmt2 = (Gm*lt)*(Gm*lt);
		double lct2 = (Gc*lt)*(Gc*lt);
		double lcz2 = (Gc*lz)*(Gc*lz);
		double lcp2 = (Gc*lp)*(Gc*lp);
		double lcn2 = (Gc*ln)*(Gc*ln);

		// compute stress
		// passive | consistent with: Sth = Jh*(Ui*sNt*Ui), sNt = phih/phio*sto, sto = 1.0/Jo*(Uo*(Sto*Uo)), Sto = phio*c*(lt2-1.0)*exp(d*(lt2-1.0)*(lt2-1.0))*(Gt*Gt)*dyad(N[1])

		mat3ds Sm = (Jh/Jo*phimh)*(cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(Uh.inverse()*(Uo*N[1])));
		mat3ds Sc = (Jh/Jo*phich)*(cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo*N[1]))*betat +
								   cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo*N[2]))*betaz +
								   cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo* Np ))*betad +
								   cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo* Nn ))*betad );

		mat3ds Sx = Se + Sm + Sc;

		S = Sx + Ci*lm*log(Jdep*J/Jh);

		// compute tangent
		mat3ds tent = dyad(F*(Uh.inverse()*(Uo*N[1])));
		mat3ds tenz = dyad(F*(Uh.inverse()*(Uo*N[2])));
		mat3ds tenp = dyad(F*(Uh.inverse()*(Uo*Np)));
		mat3ds tenn = dyad(F*(Uh.inverse()*(Uo*Nn)));

		// passive
		tens4ds cf = (Jh/Jo*phimh)*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
				(Jh/Jo*phich)*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
						2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
						2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
						2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);

		cf /= J;

		tens4ds c = ce + cf;

		c += lm/J*(IxI-2.0*log(Jdep*J/Jh)*IoI);

		css = tens4dmm(c);		// c in tens4dss form

		break;
	}
	}

	// pull back to reference configuration
	tens4dmm css_ref = J*css.pp(F.inverse());

	// convert to vector for FORTRAN
	typedef double (*ten_2nd)[3];
	typedef double (*ten_4th)[3][3][3];

	ten_2nd S2 = (ten_2nd) S_out;
	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
		{
			S2[j][i] = S(i, j);
			if (std::isnan(S(i, j)))
				std::terminate();
		}

	ten_4th C4 = (ten_4th) CC_out;
	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			for (int k=0; k < 3; k++)
				for (int l=0; l < 3; l++)
				{
					C4[l][k][j][i] = css_ref(i, j, k, l);
					if (std::isnan(css_ref(i, j, k, l)))
						std::terminate();
				}

	// store prestress state
	if (mode == prestress)
	{
		Jo = J;
		svo = 1.0/3.0/J*S.dotdot(C);
		phic = phico;
		smo = 1.0/J*(u*(Sm*u)).sym();
		sco = 1.0/J*(u*(Sc*u)).sym();
		Fio = F.inverse();

		grInt[0] = Jo;
		grInt[1] = svo;
		grInt[2] = phic;
		grInt[3] = tau;
		grInt[4] = po;
		int k = 5;
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				grInt[k] = smo(i,j);
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=i; j<3; j++)
			{
				grInt[k] = sco(i,j);
				k++;
			}
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
			{
				grInt[k] = Fio(i,j);
				k++;
			}
		// grInt[28] = tau;
		
		k = 26;
		grInt[k]      = J;
		grInt[k + 1]  = 1.0/3.0/J*S.dotdot(C);
		grInt[k + 2]  = phico;
		grInt[k + 3]  = 1.0;
		grInt[k + 4]  = p;
		grInt[k + 5]  = grInt[k + 3] - 1.0; // delta tau
		grInt[k + 6]  = grInt[k + 1] / grInt[1] - 1.0; // delta sigma
		grInt[k + 7]  = KsKi; // kski = delta sigma / deltau tau
		grInt[k + 8]  = grInt[k + 6] - KsKi * grInt[k + 5]; // ups -> 0
		grInt[k + 9]  = p0;
		grInt[k + 10] = p0 / p;
		grInt[k + 11] = phic;
	}
	// store g&r state
	if (mode == gr)
	{
		int k = 26;
		grInt[k]      = J;
		grInt[k + 1]  = 1.0/3.0/J*S.dotdot(C);
		grInt[k + 2]  = phico;
		grInt[k + 3]  = tau_ratio;
		grInt[k + 4]  = p;
		grInt[k + 5]  = grInt[k + 3] - 1.0; // delta tau
		grInt[k + 6]  = grInt[k + 1] / grInt[1] - 1.0; // delta sigma
		grInt[k + 7]  = grInt[k + 6] / grInt[k + 5]; // kski = delta sigma / deltau tau
		grInt[k + 8]  = grInt[k + 6] - KsKi * grInt[k + 5]; // ups -> 0
		grInt[k + 9]  = p0;
		grInt[k + 10] = p0 / p;
		grInt[k + 11] = phic;
		// Fih = F.inverse();
		// grInt[25] = J;
		// grInt[26] = svo;
		// grInt[27] = phico;
		// int k = 29;
		// for (int i=0; i<3; i++)
		// 	for (int j=i; j<3; j++)
		// 	{
		// 		grInt[k] = smo(i,j);
		// 		k++;
		// 	}
		// for (int i=0; i<3; i++)
		// 	for (int j=i; j<3; j++)
		// 	{
		// 		grInt[k] = sco(i,j);
		// 		k++;
		// 	}
		// for (int i=0; i<3; i++)
		// 	for (int j=0; j<3; j++)
		// 	{
		// 		grInt[k] = Fih(i,j);
		// 		k++;
		// 	}
	}
}

void stress_tangent_stvk_(const double* Fe, const double* fl, const double* time, double* eVWP, double* grInt, double* S_out, double* CC_out)
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
	tens4dmm cIxI = tens4dmm(IxI);
	tens4dmm cIoI = tens4dmm(IoI);

	// green-lagrange strain
	mat3ds E = 0.5 * ((F.transpose() * F).sym()- I);

	// stress
	mat3ds S = lambda*E.tr()*I + 2.0 * mu * E;

	// tangent
	tens4dmm css = lambda * cIxI + 2.0 * mu * cIoI;

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
