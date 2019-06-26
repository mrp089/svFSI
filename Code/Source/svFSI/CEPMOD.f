!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!-----------------------------------------------------------------------
!
!     This module defines data structures for cardiac electrophysiology
!     model equation. It also interfaces with individual modules for
!     the cellular activation model.
!
!-----------------------------------------------------------------------

      MODULE CEPMOD
      USE APMOD
      USE FNMOD
      USE TTPMOD
      USE BOMOD
      IMPLICIT NONE

!     Type of cardiac electrophysiology models: Aliev-Panfilov model,
!     tenTusscher-Panfilov model, Bueno-Orovio model
      INTEGER, PARAMETER :: cepModel_NA = 100, cepModel_AP = 101,
     2   cepModel_FN = 102, cepModel_TTP = 103, cepModel_BO = 104

!     Time integration scheme: Forward-Euler, Runge-Kutta 4th order,
!     Crank-Nicholson
      INTEGER, PARAMETER :: tIntType_NA  = 200, tIntType_FE = 201,
     2   tIntType_RK4 = 202, tIntType_CN2 = 203

!     Time integration scheme and related parameters
      TYPE odeType
!        Time integration method type
         INTEGER :: tIntType = tIntType_NA
!        Max. iterations for Newton-Raphson method
         INTEGER :: maxItr
!        Absolute tolerance
         REAL(KIND=8) :: absTol
!        Relative tolerance
         REAL(KIND=8) :: relTol
      END TYPE odeType

!     External stimulus type
      TYPE stimType
!        start time
         REAL(KIND=8) :: Ts
!        duration of stimulus
         REAL(KIND=8) :: Td
!        cycle length
         REAL(KIND=8) :: CL
!        stimulus amplitude
         REAL(KIND=8) :: A
      END TYPE stimType

!     Cardiac electrophysiology model type
      TYPE cepModelType
!        Type of cardiac electrophysiology model
         INTEGER :: cepType = cepModel_NA
!        Number of unknowns
         INTEGER :: nX
!        Number of fiber directions
         INTEGER :: nFn
!        Myocardium zone id
         INTEGER :: imyo
!        Time step for integration
         REAL(KIND=8) :: dt
!        Isotropic conductivity
         REAL(KIND=8) :: Diso = 0D0
!        Constant for stretch-activated-currents
         REAL(KIND=8) :: Ksac
!        Anisotropic conductivity
         REAL(KIND=8), ALLOCATABLE :: Dani(:)
!        External stimulus
         TYPE(stimType) :: Istim
!        Time integration options
         TYPE(odeType) :: odeS
      END TYPE cepModelType

!     Cardiac electromechanics model type
      TYPE cemModelType
!        Whether electrophysiology and mechanics are coupled
         LOGICAL :: cpld = .FALSE.
!        Whether active stress formulation is employed
         LOGICAL :: aStress = .FALSE.
!        Whether active strain formulation is employed
         LOGICAL :: aStrain = .FALSE.
!        Local variable integrated in time
!          := activation force for active stress model
!          := fiber stretch for active strain model
         REAL(KIND=8), ALLOCATABLE :: Ya(:)
      END TYPE cemModelType

!     Cardiac electromechanics type
      TYPE(cemModelType) :: cem

      END MODULE CEPMOD
!#######################################################################
