//=============================================================================
// This FEBio plugin implements a mechanobiologically equilibrated G&R material model
// based on a constrained mixture approach.
//
// Ref: Latorre M, Humphrey JD (2020) Fast, rate-independent, finite element implementation of
//      a 3D constrained mixture model of soft tissue growth and remodeling. CMAME 368, 113156.
//
// Author: Marcos Latorre (marcos.latorre@yale.edu) - 2020
//=============================================================================

//#include "stdafx.h"
//#include "FEMbeCmm.h"
//#include "FECore/FEModel.h"                             // to get current time

#include <iostream>
#include <limits>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//#include <bits/stdc++.h>
#include <mat3d.h>
#include <vec2d.h>
#include <vec3d.h>
#include <stdafx.h>

extern"C"
{
  void stress_2pk_(const double* Fe, const double* fl, double* S_out, double* CC_out)
  {
	  // convert deformation gradient to FEBio format
	  mat3d F(Fe[0], Fe[1], Fe[2], Fe[3], Fe[4], Fe[5], Fe[6], Fe[7], Fe[8]);

	  // fiber directions
	  vec3d f1(fl[0], fl[1], fl[2]);
	  vec3d f2(fl[3], fl[4], fl[5]);
	  vec3d f3 = f1 ^ f2;

	  // determinant of the deformation gradient
	  const double J = F.det();

	  for (int i=0; i<6; i++)
		  std::cout<<i<<" "<<fl[i]<<std::endl;

	  double eps = std::numeric_limits<double>::epsilon();                // machine epsilon (for floating point arithmetic)

	  // get current G&R time
	  double sgr = 0.0;//GetFEModel()->GetTime().currentTime;

//	  // retrieve material position
//	  vec3d  X = 0.0;//pt.m_r0;
//	  vec2d NX = {X.x,X.y};                                               // radial vector
//
//	  double ro = sqrt(NX*NX);                                            // radius in reference configuration (for a cylinder)
	  double ro = 1.0;
//
//	  // retrieve local element basis directions from input file
	  vec3d N[3];
//	  N[2] = pt.m_Q.col(0); // axial
//	  N[1] = pt.m_Q.col(1); // circumferential
//	  N[0] = pt.m_Q.col(2); // radial
	  N[2] = f3; // axial
	  N[1] = f2; // circumferential
	  N[0] = f1; // radial

	  // original homeostatic material parameters (ToBeDone: should be read from input file...)

	  double phieo = 0.34;                                                // mass fraction elastin
	  double phimo = 0.5*(1.0-phieo);                                     // smooth muscle cells (smc)
	  double phico = 0.5*(1.0-phieo);                                     // collagen

	  double ce = 89.71;                                                  // elastin shear modulus

	  double alpha = 0.522;                                               // original orientation of diagonal collagen

	  double cm = 261.4;                                                  // c1 smc
	  double dm = 0.24;                                                   // c2 smc
	  double Gm = 1.20;                                                   // deposition stretch smc
	  double cc = 234.9;                                                  // c1 collagen
	  double dc = 4.08;                                                   // c2 collagen
	  double Gc = 1.25;                                                   // deposition stretch collagen

	  vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);                        // original diagonal fiber direction
	  vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);                        // idem for symmetric family

	  double Get = 1.90;                                                  // circumferential deposition stretch for elastin
	  double Gez = 1.62;                                                  //      axial      deposition stretch for elastin

	  double betat = 0.056;                                               // orientation fractions for collagen, circumferential
	  double betaz = 0.067;                                               // axial
	  double betad = 0.5*(1.0-betat-betaz);                               // diagonal (both)

	  mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);   // Ge from spectral decomposition

	  double KsKi = 1.0;                                                  // shear/intramural stress gain ratio
	  double EPS  = 1.0;                                                  // blood flow rate ratio

	  double eta = 1.0;                                                   // smc/collagen turnover ratio

	  // compute U from polar decomposition of deformation gradient tensor
	  mat3ds U;
	  mat3d R;
	  F.right_polar(R,U);

	  // right Cauchy-Green tensor and its inverse
	  mat3ds C = (F.transpose() * F).sym();
	  mat3ds Ci = C.inverse();

	  // stress for elastin
	  mat3ds Se;// = phieo*ce*Ge*Ge;                                         // phieo*Ge*Sehat*Ge = phieo*Ge*(ce*I)*Ge

	  // computation of the second Piola-Kirchhoff stress
	  mat3ds S;



	  // retrieve initial Jacobian, vol. stress, smc rotated stress, collagen rotated stress, inverse of Fo, collagen mass fraction

	  double    Jo;
	  double   svo;
	  mat3ds   smo;
	  mat3ds   sco;
	  mat3d    Fio;
	  double  phic;
	  // stored in material point memory for tangent computation phase

	  // STAGE I = ELASTIC PRE-LOADING = ORIGINAL HOMEOSTATIC STATE

	  if (sgr <= 1.0 + eps) {

		  double Jdep = 0.9999;                                           // "deposition volume ratio"
		  double lm = 1.0e3*ce;                                           // bulk modulus for volumetric penalty (nearly incompressibility)

		  double lt = (F*N[1]).norm();                                    // circumferential stretch (from reference configuration)
		  double lz = (F*N[2]).norm();                                    // axial
		  double lp = (F*Np).norm();                                      // diagonal 1
		  double ln = (F*Nn).norm();                                      // diagonal 2

		  double lmt2 = (Gm*lt)*(Gm*lt);                                  // smc circumferential stretch squared
		  double lct2 = (Gc*lt)*(Gc*lt);                                  // circumferential collagen stretch squared
		  double lcz2 = (Gc*lz)*(Gc*lz);                                  //      axial      collagen stretch squared
		  double lcp2 = (Gc*lp)*(Gc*lp);                                  //    diagonal 1   collagen stretch squared
		  double lcn2 = (Gc*ln)*(Gc*ln);                                  //    diagonal 2   collagen stretch squared

		  // second Piola-Kirchhoff stress

		  mat3ds Sm = phimo * (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));           // smc
		  mat3ds Sc = phico * (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				  cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				  cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				  cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );    // collagen

		  mat3ds Sx = Se + Sm + Sc;                                       // elastin + smc + collagen + ...

		  S = Sx + Ci*lm*log(Jdep*J);                                     // ... + volumetric (penalty) contribution

		  // store initial Jacobian, vol. stress, smc rotated stress, collagen rotated stress, inverse of Fo, collagen mass fraction

		  Jo   = J;
		  svo  = 1.0/3.0/J*S.dotdot(C);
		  smo  = (U*(Sm*U)).sym();
		  smo *= 1.0/J;
		  sco  = (U*(Sc*U)).sym();
		  sco *= 1.0/J;
		  Fio  = F.inverse();
		  phic = phico;
	  }

	  // STAGE II = MECHANOBIOLOGICALLY EQUILIBRATED G&R COMPUTATION / EVOLUTION

	  else {
		  // solve nonlinear residual equation resulting from <SUM(phi)=1> and <J*phim/phimo=(J*phic/phico)^eta>: Rphi = phieo + phimo*pow(J/Jo*phic/phico,eta) + J/Jo*phic - J/Jo = 0

		  phic = phico;                                                            // initial guess
		  double dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));   // initial tangent d(R)/d(phic)
		  double Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;       // initial residue
		  do {                                                                     // local iterations to obtain phic
			  phic = phic-Rphi/dRdc;                                               // phic
			  dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));      // tangent
			  Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;          // update residue
		  } while (abs(Rphi) > sqrt(eps));                                         // && abs(Rphi/Rphi0) > sqrt(eps) && j<10
		  phic = phic-Rphi/dRdc;                                                   // converge phase -> phic (stored in material point memory for tangent phase)

		  double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);                     // phim from <J*phim/phimo=(J*phic/phico)^eta>

		  double lr = (F*(Fio*N[0])).norm();                               // lr -> 1 for F -> Fo
		  double lt = (F*(Fio*N[1])).norm();                               // lt -> 1 for F -> Fo

		  // COMPUTE CAUCHY STRESS

		  double rIo = 0.6468;                                             // initial inner radius (TBD: read it from input file...)
		  double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;                        // relative change in radius: rIrIo -> rIorIo = 1 for F -> Fo

		  mat3ds sNm = phim/phimo*smo;                                     // phim*smhato = phim*(smo/phimo) = (phim/phimo)*smo
		  mat3ds sNc = phic/phico*sco;                                     // phic*schato = phic*(sco/phico) = (phic/phico)*sco

		  mat3ds sNf = sNm + sNc;                                          // smc + collagen

		  mat3ds Ui = U.inverse();                                         // inverse of U

		  mat3ds Sf = (Ui*sNf*Ui).sym();                                       // J*Ui*sNf*Ui
		  Sf *= J;

		  mat3ds Sx = Se+Sf;                                               // elastin + smc + collagen + ...

		  double p = 1.0/3.0/J*Sx.dotdot(C) - svo*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0));     // Ups = 1 -> p (evolving Lagrange multiplier)

		  S = Sx - J*p*Ci;                                                 // ... + mechanobiologically equilibrated contribution from p
	  }

	  mat3ds s = 1.0/J*((F*(S*F.transpose()))).sym();                      // push forward S

	  // the Cauchy stress is returned

	  // convert to vector for FORTRAN
	  S_out[0] = s.xx();
	  S_out[1] = s.xy();
	  S_out[2] = s.xz();
	  S_out[3] = s.xy();
	  S_out[4] = s.yy();
	  S_out[5] = s.yz();
	  S_out[6] = s.xz();
	  S_out[7] = s.yz();
	  S_out[8] = s.zz();
  }
}


////-----------------------------------------------------------------------------
//// This function needs to return the spatial (i.e. Cauchy stress) at the material point
//// which is passed as a parameter. The FEMaterialPoint class contains all the state 
//// variables of the material (and its base classes).
//mat3ds FEMbeCmm::Stress(FEMaterialPoint& mp)
//{
//    // The FEMaterialPoint classes are stored in a linked list. The specific material
//    // point data needed by this function can be accessed using the ExtractData member.
//    // In this case, we want to FEElasticMaterialPoint data since it stores the deformation
//    // information that is needed to evaluate the stress.
//    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
//
//    // We'll need the deformation gradient and its determinant in this function.
//    mat3d &F = pt.m_F;
//    double J = pt.m_J;
//
//    double eps = std::numeric_limits<double>::epsilon();                // machine epsilon (for floating point arithmetic)
//
//    // get current G&R time
//    double sgr = GetFEModel()->GetTime().currentTime;
//
//    // retrieve material position
//    vec3d  X = pt.m_r0;
//    vec2d NX = {X.x,X.y};                                               // radial vector
//
//    double ro = sqrt(NX*NX);                                            // radius in reference configuration (for a cylinder)
//
//    // retrieve local element basis directions from input file
//    vec3d N[3];
//    N[2] = pt.m_Q.col(0); N[1] = pt.m_Q.col(1); N[0] = pt.m_Q.col(2);   // axial, circumferential, radial
//
//    // original homeostatic material parameters (ToBeDone: should be read from input file...)
//
//    double phieo = 0.34;                                                // mass fraction elastin
//    double phimo = 0.5*(1.0-phieo);                                     // smooth muscle cells (smc)
//    double phico = 0.5*(1.0-phieo);                                     // collagen
//
//    double ce = 89.71;                                                  // elastin shear modulus
//
//    double alpha = 0.522;                                               // original orientation of diagonal collagen
//
//    double cm = 261.4;                                                  // c1 smc
//    double dm = 0.24;                                                   // c2 smc
//    double Gm = 1.20;                                                   // deposition stretch smc
//    double cc = 234.9;                                                  // c1 collagen
//    double dc = 4.08;                                                   // c2 collagen
//    double Gc = 1.25;                                                   // deposition stretch collagen
//
//    vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);                        // original diagonal fiber direction
//    vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);                        // idem for symmetric family
//
//    double Get = 1.90;                                                  // circumferential deposition stretch for elastin
//    double Gez = 1.62;                                                  //      axial      deposition stretch for elastin
//
//    double betat = 0.056;                                               // orientation fractions for collagen, circumferential
//    double betaz = 0.067;                                               // axial
//    double betad = 0.5*(1.0-betat-betaz);                               // diagonal (both)
//
//    mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);   // Ge from spectral decomposition
//
//    double KsKi = 1.0;                                                  // shear/intramural stress gain ratio
//    double EPS  = 1.0;                                                  // blood flow rate ratio
//
//    double eta = 1.0;                                                   // smc/collagen turnover ratio
//
//    // compute U from polar decomposition of deformation gradient tensor
//    mat3ds U; mat3d R; F.right_polar(R,U);
//
//    // right Cauchy-Green tensor and its inverse
//    mat3ds C = pt.RightCauchyGreen();
//    mat3ds Ci = C.inverse();
//    
//    // stress for elastin
//    mat3ds Se = phieo*ce*Ge*Ge;                                         // phieo*Ge*Sehat*Ge = phieo*Ge*(ce*I)*Ge
//
//    // computation of the second Piola-Kirchhoff stress
//    mat3ds S;
//
//    // STAGE I = ELASTIC PRE-LOADING = ORIGINAL HOMEOSTATIC STATE
//
//    if (sgr <= 1.0 + eps) {
//
//        double Jdep = 0.9999;                                           // "deposition volume ratio"
//        double lm = 1.0e3*ce;                                           // bulk modulus for volumetric penalty (nearly incompressibility)
//        
//        double lt = (F*N[1]).norm();                                    // circumferential stretch (from reference configuration)
//        double lz = (F*N[2]).norm();                                    // axial
//        double lp = (F*Np).norm();                                      // diagonal 1
//        double ln = (F*Nn).norm();                                      // diagonal 2
//        
//        double lmt2 = (Gm*lt)*(Gm*lt);                                  // smc circumferential stretch squared
//        double lct2 = (Gc*lt)*(Gc*lt);                                  // circumferential collagen stretch squared
//        double lcz2 = (Gc*lz)*(Gc*lz);                                  //      axial      collagen stretch squared
//        double lcp2 = (Gc*lp)*(Gc*lp);                                  //    diagonal 1   collagen stretch squared
//        double lcn2 = (Gc*ln)*(Gc*ln);                                  //    diagonal 2   collagen stretch squared
//
//        // second Piola-Kirchhoff stress
//        
//        mat3ds Sm = phimo * (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));           // smc
//        mat3ds Sc = phico * (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
//                             cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
//                             cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
//                             cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );    // collagen
//        
//        mat3ds Sx = Se + Sm + Sc;                                       // elastin + smc + collagen + ...
//        
//        S = Sx + Ci*lm*log(Jdep*J);                                     // ... + volumetric (penalty) contribution
//        
//        // store initial Jacobian, vol. stress, smc rotated stress, collagen rotated stress, inverse of Fo, collagen mass fraction
//
//        pt.m_Jo   = J;
//        pt.m_svo  = 1.0/3.0/J*S.dotdot(C);
//        pt.m_smo  = 1.0/J*(U*(Sm*U));
//        pt.m_sco  = 1.0/J*(U*(Sc*U));
//        pt.m_Fio  = F.inverse();
//        pt.m_phic = phico;
//    }
//
//    // STAGE II = MECHANOBIOLOGICALLY EQUILIBRATED G&R COMPUTATION / EVOLUTION
//
//    else {
//
//        // retrieve initial Jacobian, vol. stress, smc rotated stress, collagen rotated stress, inverse of Fo, collagen mass fraction
//
//        double    Jo = pt.m_Jo;
//        double   svo = pt.m_svo;
//        mat3ds   smo = pt.m_smo;
//        mat3ds   sco = pt.m_sco;
//        mat3d    Fio = pt.m_Fio;
//        double &phic = pt.m_phic;                                        // stored in material point memory for tangent computation phase
//
//        // solve nonlinear residual equation resulting from <SUM(phi)=1> and <J*phim/phimo=(J*phic/phico)^eta>: Rphi = phieo + phimo*pow(J/Jo*phic/phico,eta) + J/Jo*phic - J/Jo = 0
//        
//        phic = phico;                                                            // initial guess
//        double dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));   // initial tangent d(R)/d(phic)
//        double Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;       // initial residue
//        do {                                                                     // local iterations to obtain phic
//            phic = phic-Rphi/dRdc;                                               // phic
//            dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));      // tangent
//            Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;          // update residue
//        } while (abs(Rphi) > sqrt(eps));                                         // && abs(Rphi/Rphi0) > sqrt(eps) && j<10
//        phic = phic-Rphi/dRdc;                                                   // converge phase -> phic (stored in material point memory for tangent phase)
//
//        double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);                     // phim from <J*phim/phimo=(J*phic/phico)^eta>
//
//        double lr = (F*(Fio*N[0])).norm();                               // lr -> 1 for F -> Fo
//        double lt = (F*(Fio*N[1])).norm();                               // lt -> 1 for F -> Fo
//
//        // COMPUTE CAUCHY STRESS
//        
//        double rIo = 0.6468;                                             // initial inner radius (TBD: read it from input file...)
//        double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;                        // relative change in radius: rIrIo -> rIorIo = 1 for F -> Fo
//
//        mat3ds sNm = phim/phimo*smo;                                     // phim*smhato = phim*(smo/phimo) = (phim/phimo)*smo
//        mat3ds sNc = phic/phico*sco;                                     // phic*schato = phic*(sco/phico) = (phic/phico)*sco
//
//        mat3ds sNf = sNm + sNc;                                          // smc + collagen
//
//        mat3ds Ui = U.inverse();                                         // inverse of U
//
//        mat3ds Sf = J*(Ui*sNf*Ui);                                       // J*Ui*sNf*Ui
//        
//        mat3ds Sx = Se+Sf;                                               // elastin + smc + collagen + ...
//
//        double p = 1.0/3.0/J*Sx.dotdot(C) - svo*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0));     // Ups = 1 -> p (evolving Lagrange multiplier)
//        
//        S = Sx - J*p*Ci;                                                 // ... + mechanobiologically equilibrated contribution from p
//    }
//
//    mat3ds s = 1.0/J*((F*(S*F.transpose()))).sym();                      // push forward S
//
//    // the Cauchy stress is returned
//    return s;
//}
//
////-----------------------------------------------------------------------------
//// This function calculates the spatial tangent tensor. 
//// It takes one parameter, the FEMaterialPoint and returns a tens4ds object
//// which is a fourth-order tensor with minor symmetries only.
//tens4dss FEMbeCmm::Tangent(FEMaterialPoint& mp)
//{
//    // As in the Stress function, we need the data from the FEElasticMaterialPoint
//    // class to calculate the tangent.
//    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
//
//    // Get the deformation gradient and its determinant
//    mat3d &F = pt.m_F;
//    double J = pt.m_J;
//
//    double eps = std::numeric_limits<double>::epsilon();                // machine epsilon (for floating point arithmetic)
//
//    // get current G&R time
//    double sgr = GetFEModel()->GetTime().currentTime;
//
//    // retrieve material position
//    vec3d  X = pt.m_r0;
//    vec2d NX = {X.x,X.y};                                               // radial vector
//
//    double ro = sqrt(NX*NX);                                            // radius in reference configuration (for a cylinder)
//
//    // retrieve local element basis directions from input file
//    vec3d N[3];
//    N[2] = pt.m_Q.col(0); N[1] = pt.m_Q.col(1); N[0] = pt.m_Q.col(2);   // axial, circumferential, radial
//
//    // original homeostatic material parameters (TBD: should be read from input file...)
//
//    double phieo = 0.34;                                                // mass fraction elastin
//    double phimo = 0.5*(1.0-phieo);                                     // smooth muscle cells (smc)
//    double phico = 0.5*(1.0-phieo);                                     // collagen
//
//    double ce = 89.71;                                                  // elastin modulus
//
//    double alpha = 0.522;                                               // original orientation of diagonal collagen
//
//    double cm = 261.4;                                                  // c1 smc
//    double dm = 0.24;                                                   // c2 smc
//    double Gm = 1.20;                                                   // deposition stretch smc
//    double cc = 234.9;                                                  // c1 collagen
//    double dc = 4.08;                                                   // c2 collagen
//    double Gc = 1.25;                                                   // deposition stretch collagen
//
//    vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);                        // original diagonal fiber direction
//    vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);                        // idem for symmetric family
//
//    double Get = 1.90;                                                  // circumferential deposition stretch for elastin
//    double Gez = 1.62;                                                  //      axial      deposition stretch for elastin
//
//    double betat = 0.056;                                               // orientation fractions for collagen, circumferential
//    double betaz = 0.067;                                               // axial
//    double betad = 0.5*(1.0-betat-betaz);                               // diagonal (both)
//
//    mat3ds Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);    // Ge from spectral decomposition
//
//    double KsKi = 1.0;                                                  // shear/intramural stress gain ratio
//    double EPS  = 1.0;                                                  // blood flow rate ratio
//
//    double eta = 1.0;                                                   // smc/collagen turnover ratio
//
//    // define identity tensor and other dyadic products of the identity tensor
//    mat3dd  I(1.0);
//    tens4ds IxI = dyad1s(I);
//    tens4ds IoI = dyad4s(I);
//
//    // define right Cauchy-Green tensor, its inverse, and other requisite dyadic products
//    mat3ds  C  = pt.RightCauchyGreen();
//    mat3ds  Ci = C.inverse();
//    tens4ds CixCi = dyad1s(Ci);
//    tens4ds CioCi = dyad4s(Ci);
//
//    // spatial moduli for elastin (zero for neo-Hookean)
//    tens4ds ce(0.0);                                                    // phieo/J*(FcF:GecGe:Cehat:GecGe:FTcFT) = phieo/J*(FcF:GecGe:0:GecGe:FTcFT)
//
//    // computation of spatial moduli
//    mat3ds sfpro;
//    double eigenval[3]; vec3d eigenvec[3];
//
//    // STAGE I = ELASTIC PRE-LOADING = ORIGINAL HOMEOSTATIC STATE
//
//    if (sgr <= 1.0 + eps) {
//
//        double Jdep = 0.9999;                                           // "deposition volume ratio"
//        double lm = 1.0e3*ce;                                           // bulk modulus for volumetric penalty (nearly incompressibility)
//        
//        double lt = (F*N[1]).norm();                                    // circumferential stretch (from reference configuration)
//        double lz = (F*N[2]).norm();                                    // axial
//        double lp = (F*Np).norm();                                      // diagonal 1
//        double ln = (F*Nn).norm();                                      // diagonal 2
//        
//        double lmt2 = (Gm*lt)*(Gm*lt);                                  // smc circumferential stretch squared
//        double lct2 = (Gc*lt)*(Gc*lt);                                  // circumferential collagen stretch squared
//        double lcz2 = (Gc*lz)*(Gc*lz);                                  //      axial      collagen stretch squared
//        double lcp2 = (Gc*lp)*(Gc*lp);                                  //    diagonal 1   collagen stretch squared
//        double lcn2 = (Gc*ln)*(Gc*ln);                                  //    diagonal 2   collagen stretch squared
//
//
//        mat3ds tent = dyad(F*N[1]);                                     // dyadic products
//        mat3ds tenz = dyad(F*N[2]);
//        mat3ds tenp = dyad(F*Np);
//        mat3ds tenn = dyad(F*Nn);
//
//        tens4ds cf = phimo*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +    // smc
//                     phico*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
//                            2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
//                            2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
//                            2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);    // collagen
//
//        cf /= J;                                                        // divide by Jacobian
//
//        tens4ds c = ce + cf;                                            // elastin + smc + collagen + ...
//
//        // the spatial tangent tensor is returned
//        c += lm/J*(IxI-2.0*log(Jdep*J)*IoI);                            // ... + volumetric (penalty) contribution
//    }
//
//    // STAGE II = MECHANOBIOLOGICALLY EQUILIBRATED G&R COMPUTATION / EVOLUTION
//
//    else
//
//        // retrieve initial Jacobian, vol. stress, smc rotated stress, collagen rotated stress, inverse of Fo, collagen mass fraction
//
//        double   Jo = pt.m_Jo;
//        double  svo = pt.m_svo;
//        mat3ds  smo = pt.m_smo;
//        mat3ds  sco = pt.m_sco;
//        mat3d   Fio = pt.m_Fio;
//        double phic = pt.m_phic;                                        // converged phic from stress computation phase
//
//        double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);            // phim from <J*phim/phimo=(J*phic/phico)^eta>
//
//        // compute U, eigenvalues, eigenvectors, and inverse from polar decomposition of deformation gradient tensor
//        mat3ds U; mat3d R; F.right_polar(R,U);
//
//        U.eigen2(eigenval,eigenvec);
//        mat3ds Ui = U.inverse(); 
//
//        // compute current rotated stresses
//        mat3ds sNm = phim/phimo*smo;                                    // phim*smhato = phim*(smo/phimo) = (phim/phimo)*smo
//        mat3ds sNc = phic/phico*sco;                                    // phic*schato = phic*(sco/phico) = (phic/phico)*sco
//
//        double lr = (F*(Fio*N[0])).norm();                              // lr -> 1 for F -> Fo
//        double lt = (F*(Fio*N[1])).norm();                              // lt -> 1 for F -> Fo
//
//        double rIo = 0.6468;                                            // initial inner radius (TBD: read it from input file...)
//        double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;                       // relative change in radius: rIrIo -> rIorIo = 1 for F -> Fo
//
//        // 2nd P-K stresses
//        mat3ds Se = phieo*ce*Ge*Ge;                                     // phieo*Ge*Sehat*Ge = phieo*Ge*(ce*I)*Ge
//        mat3ds Sm = J*(ui*sNm*ui).sym();                                // J*Ui*sNm*Ui
//        mat3ds Sc = J*(ui*sNc*ui).sym();                                // J*Ui*sNc*Ui
//        mat3ds Sx = Se+Sm+Sc;                                           // elastin + smc + collagen
//
//        // associated Cauchy stresses
//        mat3ds sm = 1.0/J*(F*(Sm*F.transpose())).sym();
//        mat3ds sc = 1.0/J*(F*(Sc*F.transpose())).sym();
//        mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();
//
//        double p = 1.0/3.0/J*Sx.dotdot(C)-svo*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0));    // Ups = 1 -> p (evolving Lagrange multiplier)
//
//        // COMPUTE SPATIAL TANGENT
//
//        tens4ds Ixsx = dyad1s(I,sx);
//        tens4ds smxI = dyad1s(sm,I);
//        tens4ds scxI = dyad1s(sc,I);
//        
//        mat3ds tenr = dyad(F*(Fio*N[0]));                               // Fio needed for numerical consistency (from computation of lr) 
//        mat3ds tent = dyad(F*(Fio*N[1]));
//        
//        tens4ds Ixnrr = dyad1s(I,tenr);
//        tens4ds Ixntt = dyad1s(I,tent);
//
//        // contribution due to constant rotated Cauchy stresses at constituent level (see Appendix D in Latorre & Humphrey, CMAME 2020)
//
//        tens4ds cf(0.0);
//
//        // components of stress in Lagrangian strain basis
//        sfpro.zero();
//        sfpro(0,0) = eigenvec[0]*((sNm+sNc)*eigenvec[0]);
//        sfpro(1,1) = eigenvec[1]*((sNm+sNc)*eigenvec[1]);
//        sfpro(2,2) = eigenvec[2]*((sNm+sNc)*eigenvec[2]);
//        sfpro(0,1) = eigenvec[0]*((sNm+sNc)*eigenvec[1]);
//        sfpro(1,2) = eigenvec[1]*((sNm+sNc)*eigenvec[2]);
//        sfpro(0,2) = eigenvec[0]*((sNm+sNc)*eigenvec[2]);
//
//        vec3d Fxeigenvec[3];
//
//        // push forward material eigenvectors
//        Fxeigenvec[0] = F*eigenvec[0];
//        Fxeigenvec[1] = F*eigenvec[1];
//        Fxeigenvec[2] = F*eigenvec[2];
//
//        // compute spatial tangent contribution in spectral form
//        for (int i=0; i<3; i++) {
//
//            mat3ds ten1 = dyad(Fxeigenvec[i]);
//
//            for (int j=0; j<3; j++) {
//
//                double component = sfpro(i,j) / pow(eigenval[i],3) / eigenval[j];
//
//                mat3ds ten2 = dyads(Fxeigenvec[i],Fxeigenvec[j]);
//
//                cf -= component*dyad1ss(ten2,ten1);
//
//                for (int k=0; k<3; k++) {
//
//                    if (k == i) continue;
//
//                    mat3ds ten3 = dyads(Fxeigenvec[j],Fxeigenvec[k]);
//                    mat3ds ten4 = dyads(Fxeigenvec[k],Fxeigenvec[i]);
//
//                    component = sfpro(i,j) / eigenval[i] / eigenval[j] / eigenval[k] / (eigenval[i] + eigenval[k]);
//
//                    cf -= component*dyad1ss(ten3,ten4);
//                }
//            }
//        }
//
//        // additional contributions to tangent by smc and collagen
//
//        // derivatives of referential mass fractions with 
//        double dphiRm = phimo*eta*pow(J/Jo*phic/phico,eta-1.0)/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);
//        double dphiRc = phico/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico);
//
//        cf += dphiRm/phim*smxI + dphiRc/phic*scxI;
//
//        // total tangent
//
//        c = ce + cf;                                                    // elastin + smc + collagen + ...
//        
//        c += 1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddotss(IxIss,css))
//           + svo*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam)*(IxIss-2.0*IoIss)
//           - 3.0*svo*KsKi*EPS*pow(rIrIo,-4)*(ro/rIo/lt*Ixntt-(ro-rIo)/rIo/lr*Ixnrr);    // ... + mechanobiologically equilibrated contribution from p
//    }
//
//    // the spatial tangent tensor is returned
//    return c;
//}
