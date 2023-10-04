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
!--------------------------------------------------------------------
!
!     This routines is for solving nonlinear structural mechanics
!     problem (pure displacement-based formulation).
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_dSOLID(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys, iFn, nFn
      REAL(KIND=RKIND) w, Jac, grInt(nGrInt), ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:), N(:),
     3   Nx(:,:), lR(:,:), lK(:,:,:), lVWP(:,:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     STRUCT: dof = nsd
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN), N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN),
     4   lK(dof*dof,eNoN,eNoN), lVWP(nvwp,eNoN))


!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_struct .AND.
     2       cPhys .NE. phys_gr) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

      !   CALL EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptr, lK, lR)
         CALL EVAL_dSOLID_FD(e, lM, Ag, Yg, Dg, ptr, lK, lR)

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      ! WRITE(*,*) lK(1,:,:)
      ! WRITE(*,*) lR
      ! CALL EXIT(0)

      DEALLOCATE(ptr, xl, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, N, Nx,
     2   lR, lK, lVWP)

      RETURN
      END SUBROUTINE CONSTRUCT_dSOLID

      SUBROUTINE EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptr, lK, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys, iFn, nFn
      REAL(KIND=RKIND) w, Jac, grInt(nGrInt), ksix(nsd,nsd)

      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,lM%eNoN),
     2   lK(dof*dof,lM%eNoN,lM%eNoN)

      INTEGER(KIND=IKIND), INTENT(INOUT) :: ptr(lM%eNoN)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:), N(:),
     3   Nx(:,:), lVWP(:,:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     STRUCT: dof = nsd
      ALLOCATE(xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN), N(eNoN), Nx(nsd,eNoN),
     4   lVWP(nvwp,eNoN))

!        Create local copies
         fN   = 0._RKIND
         pS0l = 0._RKIND
         ya_l = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)

!           Variable wall - SCHWARZ July 2021---------------------------
!           Calculate local wall property
            IF (useVarWall) lVWP(:,a) = vWP0(:,Ac)
!           ------------------------------------------------------------
            IF (ALLOCATED(lM%fN)) THEN
               DO iFn=1, nFn
                  fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
               END DO
            END IF
            IF (ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
            IF (cem%cpld) ya_l(a) = cem%Ya(Ac)
         END DO

!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g) * Jac
            N = lM%N(:,g)

!           retrieve g&r internal variables
            grInt(:) = 0._RKIND
            IF (ALLOCATED(lM%grVn)) grInt(1:nGrInt) = lM%grVo(:,g,e)

            pSl = 0._RKIND
            IF (nsd .EQ. 3) THEN
               CALL STRUCT3D(eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2            pS0l, pSl, ya_l, lR, lK, grInt, lVWP)

            ELSE IF (nsd .EQ. 2) THEN
               CALL STRUCT2D(eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2            pS0l, pSl, ya_l, lR, lK, grInt)

            END IF

!           Prestress
            IF (pstEq) THEN
               DO a=1, eNoN
                  Ac = ptr(a)
                  pSn(:,Ac) = pSn(:,Ac) + w*N(a)*pSl(:)
                  pSa(Ac)   = pSa(Ac)   + w*N(a)
               END DO
            END IF

!           Update g&r variables
            IF (ALLOCATED(lM%grVn)) lM%grVo(:,g,e) = grInt(1:nGrInt)
            IF (ALLOCATED(lM%grVn)) lM%grVn(:,g,e) = grInt(1:nGrInt)

         END DO ! g: loop

!     ptr, lR, lK,
      DEALLOCATE(xl, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, N, Nx, lVWP)

      RETURN
      END SUBROUTINE EVAL_dSOLID


!####################################################################


      SUBROUTINE EVAL_dSOLID_FD(e, lM, Ag, Yg, Dg, ptr, lK, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(INOUT) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo), lR(dof,lM%eNoN), lK(dof*dof,lM%eNoN,lM%eNoN)

      INTEGER(KIND=IKIND) e, Ac, eNoN
      REAL(KIND=RKIND) w, Jac, grInt(nGrInt), ksix(nsd,nsd)

      INTEGER(KIND=IKIND), INTENT(INOUT) :: ptr(lM%eNoN)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptrtmp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lKtmp(:,:,:), lRp(:,:), dlR(:,:),
     2   du
      INTEGER(KIND=IKIND) ii, jj, kk, ll, dd
      REAL(KIND=RKIND) :: fa, fy, fd

      eNoN = lM%eNoN

!     STRUCT: dof = nsd
      ALLOCATE(ptrtmp(lM%eNoN), lRp(dof,eNoN), lKtmp(dof*dof,eNoN,eNoN))

!     time integration factors
      fd = eq(cEq)%af*eq(cEq)%beta*dt*dt
      fy = eq(cEq)%af*eq(cEq)%gam*dt
      fa = eq(cEq)%am

!     initialize
      lK = 0._RKIND

!     central evaluation
      CALL EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptr, lK, lR)

!     loop nodes
      DO ii=1,dof
        DO jj=1,eNoN
          du = 1.E-8_RKIND

!         global node id
          Ac = lM%IEN(jj,e)

!         perturb displacement vector
          Dg(ii,Ac) = Dg(ii,Ac) + du

          CALL EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore displacement vector
          Dg(ii,Ac) = Dg(ii,Ac) - du
!
!         calculate finite difference
          dlR = fd * (lRp - lR) / du

!         perturb velocity vector
          Yg(ii,Ac) = Yg(ii,Ac) + du

          CALL EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore velocity vector
          Yg(ii,Ac) = Yg(ii,Ac) - du
!
!         calculate finite difference
          dlR = dlR + fy * (lRp - lR) / du

!         perturb acceleration vector
          Ag(ii,Ac) = Ag(ii,Ac) + du

          CALL EVAL_dSOLID(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore acceleration vector
          Ag(ii,Ac) = Ag(ii,Ac) - du
!
!         calculate finite difference
          dlR = dlR + fa * (lRp - lR) / du

!         assign to tangent matrix
          DO kk=1,dof
            dd = (kk - 1) * dof + ii
            DO ll=1,lM%eNoN
              lK(dd,ll,jj) = dlR(kk,ll)
            END DO
          END DO

        END DO
      END DO
      RETURN
      END SUBROUTINE EVAL_dSOLID_FD


!####################################################################
      SUBROUTINE STRUCT3D(eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2   pS0l, pSl, ya_l, lR, lK, grInt, lVWP)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(3,eNoN),
     3   fN(3,nFn), pS0l(6,eNoN), ya_l(eNoN), lVWP(nvwp,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: pSl(6)
      REAL(KIND=RKIND), INTENT(INOUT) :: grInt(nGrInt), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j, k, ii, jj, dd
      REAL(KIND=RKIND) :: rho, dmp, T1, amd, afl, ya_g, fb(3), ud(3),
     2   NxSNx, BmDBm, F(3,3), S(3,3), P(3,3), Dm(6,6), DBm(6,3),
     3   Bm(6,3,eNoN), S0(3,3), eVWP(nvwp), p_equi, stim

!     Define parameters
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp     = eq(cEq)%dmn(cDmn)%prop(damping)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      amd     = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     Inertia, body force and deformation tensor (F)
      ud     = -rho*fb
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      S0     = 0._RKIND
      ya_g   = 0._RKIND
      eVWP   = 0._RKIND
      p_equi = 0._RKIND

      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))
         ud(3) = ud(3) + N(a)*(rho*(al(k,a)-bfl(3,a)) + dmp*yl(k,a))

         F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
         F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
         F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
         F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
         F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
         F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
         F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)

!     Variable wall - SCHWARZ July 2021---------------------------------
!     Calculate local wall property
         IF (useVarWall) eVWP(:) = eVWP(:) + N(a)*lVWP(:,a)
!     ------------------------------------------------------------------

         S0(1,1) = S0(1,1) + N(a)*pS0l(1,a)
         S0(2,2) = S0(2,2) + N(a)*pS0l(2,a)
         S0(3,3) = S0(3,3) + N(a)*pS0l(3,a)
         S0(1,2) = S0(1,2) + N(a)*pS0l(4,a)
         S0(2,3) = S0(2,3) + N(a)*pS0l(5,a)
         S0(3,1) = S0(3,1) + N(a)*pS0l(6,a)

         ya_g    = ya_g + N(a)*ya_l(a)

!        interpolate lagrange multiplier
         p_equi = p_equi + N(a) * dl(4,a)
      END DO
      S0(2,1) = S0(1,2)
      S0(3,2) = S0(2,3)
      S0(1,3) = S0(3,1)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
!     Voigt notationa (Dm)
      CALL GETPK2CC(eq(cEq)%dmn(cDmn), F, nFn, fN, ya_g, grInt, S, Dm,
     2              eVWP, p_equi, stim)

!     Prestress
      pSl(1) = S(1,1)
      pSl(2) = S(2,2)
      pSl(3) = S(3,3)
      pSl(4) = S(1,2)
      pSl(5) = S(2,3)
      pSl(6) = S(3,1)
      S = S + S0

!     1st Piola-Kirchhoff tensor (P)
      P    = MATMUL(F, S)

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)
         Bm(1,3,a) = Nx(1,a)*F(3,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)
         Bm(2,3,a) = Nx(2,a)*F(3,2)

         Bm(3,1,a) = Nx(3,a)*F(1,3)
         Bm(3,2,a) = Nx(3,a)*F(2,3)
         Bm(3,3,a) = Nx(3,a)*F(3,3)

         Bm(4,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(4,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
         Bm(4,3,a) = (Nx(1,a)*F(3,2) + F(3,1)*Nx(2,a))

         Bm(5,1,a) = (Nx(2,a)*F(1,3) + F(1,2)*Nx(3,a))
         Bm(5,2,a) = (Nx(2,a)*F(2,3) + F(2,2)*Nx(3,a))
         Bm(5,3,a) = (Nx(2,a)*F(3,3) + F(3,2)*Nx(3,a))

         Bm(6,1,a) = (Nx(3,a)*F(1,1) + F(1,3)*Nx(1,a))
         Bm(6,2,a) = (Nx(3,a)*F(2,1) + F(2,3)*Nx(1,a))
         Bm(6,3,a) = (Nx(3,a)*F(3,1) + F(3,3)*Nx(1,a))
      END DO

!     Local residue and tangent matrices
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(1,a)*P(1,1) +
     2      Nx(2,a)*P(1,2) + Nx(3,a)*P(1,3))
         lR(2,a) = lR(2,a) + w*(N(a)*ud(2) + Nx(1,a)*P(2,1) +
     2      Nx(2,a)*P(2,2) + Nx(3,a)*P(2,3))
         lR(3,a) = lR(3,a) + w*(N(a)*ud(3) + Nx(1,a)*P(3,1) +
     2      Nx(2,a)*P(3,2) + Nx(3,a)*P(3,3))

         DO b=1, eNoN
!           Geometric stiffness
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(2,a)*S(2,1)*Nx(1,b) +
     2              Nx(3,a)*S(3,1)*Nx(1,b) + Nx(1,a)*S(1,2)*Nx(2,b) +
     3              Nx(2,a)*S(2,2)*Nx(2,b) + Nx(3,a)*S(3,2)*Nx(2,b) +
     4              Nx(1,a)*S(1,3)*Nx(3,b) + Nx(2,a)*S(2,3)*Nx(3,b) +
     5              Nx(3,a)*S(3,3)*Nx(3,b)
            T1 = amd*N(a)*N(b) + afl*NxSNx

!           Material Stiffness (Bt*D*B)
            DBm = MATMUL(Dm, Bm(:,:,b))

            BmDBm = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     2              Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            lK(1,a,b) = lK(1,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2              Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     2              Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            lK(2,a,b) = lK(2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2              Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     2              Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            lK(3,a,b) = lK(3,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2              Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     2              Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            lK(4,a,b) = lK(dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2              Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     2              Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            lK(5,a,b) = lK(dof+2,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2              Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     2              Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            lK(6,a,b) = lK(dof+3,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2              Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     2              Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            lK(7,a,b) = lK(2*dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2              Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     2              Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            lK(8,a,b) = lK(2*dof+2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2              Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     2              Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            lK(9,a,b) = lK(2*dof+3,a,b) + w*(T1 + afl*BmDBm)
         END DO
      END DO

      IF (eq(cEq)%dmn(cDmn)%phys .EQ. phys_gr) THEN
!        residual
         DO a=1, eNoN
            lR(4,a) = lR(4,a) + N(a) * (p_equi - stim)
         END DO
      END IF

      RETURN
      END SUBROUTINE STRUCT3D
!####################################################################
      SUBROUTINE STRUCT2D(eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2   pS0l, pSl, ya_l, lR, lK, grInt)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(2,eNoN),
     3   fN(2,nFn), pS0l(3,eNoN), ya_l(eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: pSl(3)
      REAL(KIND=RKIND), INTENT(INOUT) :: grInt(nGrInt), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j
      REAL(KIND=RKIND) :: rho, dmp, T1, amd, afl, ya_g, fb(2), ud(2),
     2   NxSNx, BmDBm, F(2,2), S(2,2), P(2,2), Dm(3,3), DBm(3,2),
     3   Bm(3,2,eNoN), S0(2,2)

!     Define parameters
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp     = eq(cEq)%dmn(cDmn)%prop(damping)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      amd     = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i       = eq(cEq)%s
      j       = i + 1

!     Inertia, body force and deformation tensor (F)
      ud     = -rho*fb
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      S0     = 0._RKIND
      ya_g   = 0._RKIND
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))

         F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
         F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)

         S0(1,1) = S0(1,1) + N(a)*pS0l(1,a)
         S0(2,2) = S0(2,2) + N(a)*pS0l(2,a)
         S0(1,2) = S0(1,2) + N(a)*pS0l(3,a)

         ya_g    = ya_g + N(a)*ya_l(a)
      END DO
      S0(2,1) = S0(1,2)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
!     Voigt notation
      CALL GETPK2CC(eq(cEq)%dmn(cDmn), F, nFn, fN, ya_g, grInt, S, Dm)

!     Prestress
      pSl(1) = S(1,1)
      pSl(2) = S(2,2)
      pSl(3) = S(1,2)
      S = S + S0

!     1st Piola-Kirchhoff tensor (P)
      P = MATMUL(F, S)

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)

         Bm(3,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(3,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
      END DO

!     Local residue and tangent matrices
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(1,a)*P(1,1) +
     2      Nx(2,a)*P(1,2))
         lR(2,a) = lR(2,a) + w*(N(a)*ud(2) + Nx(1,a)*P(2,1) +
     2      Nx(2,a)*P(2,2))

         DO b=1, eNoN
!           Geometric stiffness
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(2,a)*S(2,1)*Nx(1,b) +
     2              Nx(1,a)*S(1,2)*Nx(2,b) + Nx(2,a)*S(2,2)*Nx(2,b)
            T1 = amd*N(a)*N(b) + afl*NxSNx

!           Material stiffness (Bt*D*B)
            DBm(1,1) = Dm(1,1)*Bm(1,1,b) + Dm(1,2)*Bm(2,1,b) +
     2         Dm(1,3)*Bm(3,1,b)
            DBm(1,2) = Dm(1,1)*Bm(1,2,b) + Dm(1,2)*Bm(2,2,b) +
     2         Dm(1,3)*Bm(3,2,b)

            DBm(2,1) = Dm(2,1)*Bm(1,1,b) + Dm(2,2)*Bm(2,1,b) +
     2         Dm(2,3)*Bm(3,1,b)
            DBm(2,2) = Dm(2,1)*Bm(1,2,b) + Dm(2,2)*Bm(2,2,b) +
     2         Dm(2,3)*Bm(3,2,b)

            DBm(3,1) = Dm(3,1)*Bm(1,1,b) + Dm(3,2)*Bm(2,1,b) +
     2         Dm(3,3)*Bm(3,1,b)
            DBm(3,2) = Dm(3,1)*Bm(1,2,b) + Dm(3,2)*Bm(2,2,b) +
     2         Dm(3,3)*Bm(3,2,b)

            BmDBm = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2         Bm(3,1,a)*DBm(3,1)
            lK(1,a,b) = lK(1,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2         Bm(3,1,a)*DBm(3,2)
            lK(2,a,b) = lK(2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2         Bm(3,2,a)*DBm(3,1)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2         Bm(3,2,a)*DBm(3,2)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*(T1 + afl*BmDBm)
         END DO
      END DO

      RETURN
      END SUBROUTINE STRUCT2D
!####################################################################
      SUBROUTINE BSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(3)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: af, h, Jac, wl, Ku, F(3,3), Fi(3,3), nFi(3),
     2   NxFi(3,eNoN)

      af     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i      = eq(cEq)%s
      j      = i + 1
      k      = j + 1

      h      = 0._RKIND
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      DO a=1, eNoN
         h       = h      + N(a)*hl(a)
         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(1,3)  = F(1,3) + Nx(3,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)
         F(2,3)  = F(2,3) + Nx(3,a)*dl(j,a)
         F(3,1)  = F(3,1) + Nx(1,a)*dl(k,a)
         F(3,2)  = F(3,2) + Nx(2,a)*dl(k,a)
         F(3,3)  = F(3,3) + Nx(3,a)*dl(k,a)
      END DO
      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

      nFi(1) = nV(1)*Fi(1,1) + nV(2)*Fi(2,1) + nV(3)*Fi(3,1)
      nFi(2) = nV(1)*Fi(1,2) + nV(2)*Fi(2,2) + nV(3)*Fi(3,2)
      nFi(3) = nV(1)*Fi(1,3) + nV(2)*Fi(2,3) + nV(3)*Fi(3,3)

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1) + Nx(3,a)*Fi(3,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2) + Nx(3,a)*Fi(3,2)
         NxFi(3,a) = Nx(1,a)*Fi(1,3) + Nx(2,a)*Fi(2,3) + Nx(3,a)*Fi(3,3)
      END DO

      wl = w*Jac*h
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)
         lR(3,a) = lR(3,a) - wl*N(a)*nFi(3)

         DO b=1, eNoN
            Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lK(2,a,b)     = lK(2,a,b)     + Ku
            lK(dof+1,a,b) = lK(dof+1,a,b) - Ku

            Ku = wl*af*N(a)*(nFi(3)*NxFi(1,b) - nFi(1)*NxFi(3,b))
            lK(3,a,b)       = lK(3,a,b)       + Ku
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) - Ku

            Ku = wl*af*N(a)*(nFi(3)*NxFi(2,b) - nFi(2)*NxFi(3,b))
            lK(dof+3,a,b)   = lK(dof+3,a,b)   + Ku
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) - Ku
         END DO
      END DO

      RETURN
      END SUBROUTINE BSTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE BSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, a, b
      REAL(KIND=RKIND) :: af, h, Jac, wl, Ku, F(2,2), Fi(2,2), nFi(2),
     2   NxFi(2,eNoN)

      af     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i      = eq(cEq)%s
      j      = i + 1

      h      = 0._RKIND
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      DO a=1, eNoN
         h       = h      + N(a)*hl(a)
         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)
      END DO
      Jac = MAT_DET(F, 2)
      Fi  = MAT_INV(F, 2)

      nFi(1) = nV(1)*Fi(1,1) + nV(2)*Fi(2,1)
      nFi(2) = nV(1)*Fi(1,2) + nV(2)*Fi(2,2)

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2)
      END DO

      wl = w*Jac*h
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)

         DO b=1, eNoN
            Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lK(2,a,b)     = lK(2,a,b)     + Ku
            lK(dof+1,a,b) = lK(dof+1,a,b) - Ku
         END DO
      END DO

      RETURN
      END SUBROUTINE BSTRUCT2D
!####################################################################
