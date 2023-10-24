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

      ALLOCATE(ptrtmp(lM%eNoN), lRp(dof,eNoN), lKtmp(dof*dof,eNoN,eNoN))

!     time integration factors
      fd = eq(cEq)%af*eq(cEq)%beta*dt*dt
      fy = eq(cEq)%af*eq(cEq)%gam*dt
      fa = eq(cEq)%am

!     initialize
      lK = 0._RKIND

!     central evaluation
      CALL EVAL_dSOLID_GR(e, lM, Ag, Yg, Dg, ptr, lK, lR)

!     loop nodes
      DO ii=1,dof
            DO jj=1,eNoN
            du = 1.E-8_RKIND

!         global node id
            Ac = lM%IEN(jj,e)

!         perturb displacement vector
            Dg(ii,Ac) = Dg(ii,Ac) + du

            CALL EVAL_dSOLID_GR(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore displacement vector
            Dg(ii,Ac) = Dg(ii,Ac) - du

!         calculate finite difference
            dlR = fd * (lRp - lR) / du

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


      SUBROUTINE EVAL_dSOLID_GR(e, lM, Ag, Yg, Dg, ptr, lK, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys, iFn, nFn
      REAL(KIND=RKIND) wd, wp, Jac, grInt(nGrInt), ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,lM%eNoN),
     2   lK(dof*dof,lM%eNoN,lM%eNoN)

      INTEGER(KIND=IKIND), INTENT(INOUT) :: ptr(lM%eNoN)
      INTEGER(KIND=IKIND) ifs

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:), 
     3   lVWP(:,:), Nd(:), Nxd(:,:), Np(:), Nxp(:,:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     STRUCT: dof = nsd
      ALLOCATE(xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN), lVWP(nvwp,eNoN),
     4   Nd(eNoN), Nxd(nsd,eNoN), Np(eNoN), Nxp(nsd,eNoN))

      CALL GETTHOODFS(fs, lM, 0, 3)
      
!     Create local copies
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

         IF (useVarWall) lVWP(:,a) = vWP0(:,Ac)
         IF (ALLOCATED(lM%fN)) THEN
            DO iFn=1, nFn
               fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
            END DO
         END IF
         IF (ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
         IF (cem%cpld) ya_l(a) = cem%Ya(Ac)
      END DO

!     Gauss integration
      lR = 0._RKIND
      lK = 0._RKIND

!     Loop shape functions
      DO g=1, fs(1)%nG
!        displacement and pressure shape functions
         CALL GNN(eNoN, nsd, fs(1)%Nx(:,:,g), xl, Nxd, Jac, ksix)
         IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
         Nxp(:,:) = fs(2)%Nx(:,:,1)
         
         wd = fs(1)%w(g) * Jac
         wp = fs(2)%w(1) * Jac
         Nd = fs(1)%N(:,g)
         Np = fs(2)%N(:,1)

!        retrieve g&r internal variables
         grInt(:) = 0._RKIND
         IF (ALLOCATED(lM%grVn)) grInt(1:nGrInt) = lM%grVo(:,g,e)

         pSl = 0._RKIND
         IF (nsd .EQ. 3) THEN
            CALL STRUCT3D_GR(eNoN, nFn, wd, wp, Nd, Np, Nxd, Nxp,
     2                       al, yl, dl, bfl, fN, pS0l, pSl, 
     3                       ya_l, lR, lK, grInt, lVWP, ifs)
         END IF

!        Update g&r variables
         IF (ALLOCATED(lM%grVn)) lM%grVo(:,g,e) = grInt(1:nGrInt)
         IF (ALLOCATED(lM%grVn)) lM%grVn(:,g,e) = grInt(1:nGrInt)
      END DO ! g: loop

!     ptr, lR, lK,
      DEALLOCATE(xl, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, 
     2            Nd, Nxd, Np, Nxp, lVWP)

      RETURN
      END SUBROUTINE EVAL_dSOLID_GR


      SUBROUTINE STRUCT3D_GR(eNoN, nFn, wd, wp, Nd, Np, Nxd, Nxp,
     2 al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK, grInt, lVWP, ifs)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn, ifs
      REAL(KIND=RKIND), INTENT(IN) :: wd, wp, Nd(eNoN), Np(eNoN),
     3   Nxd(3,eNoN), Nxp(3,eNoN),
     4   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(3,eNoN),
     5   fN(3,nFn), pS0l(6,eNoN), ya_l(eNoN), lVWP(nvwp,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: pSl(6)
      REAL(KIND=RKIND), INTENT(INOUT) :: grInt(nGrInt), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j, k, ii, jj, dd
      REAL(KIND=RKIND) :: rho, dmp, T1, amd, afl, ya_g, fb(3), ud(3),
     2   NxSNx, BmDBm, F(3,3), S(3,3), P(3,3), Dm(6,6), DBm(6,3),
     3   Bm(6,3,eNoN), eVWP(nvwp), p_equi, stim

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
      ud     = 0._RKIND
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      ya_g   = 0._RKIND
      eVWP   = 0._RKIND
      p_equi = 0._RKIND

      DO a=1, eNoN
         F(1,1) = F(1,1) + Nxd(1,a)*dl(i,a)
         F(1,2) = F(1,2) + Nxd(2,a)*dl(i,a)
         F(1,3) = F(1,3) + Nxd(3,a)*dl(i,a)
         F(2,1) = F(2,1) + Nxd(1,a)*dl(j,a)
         F(2,2) = F(2,2) + Nxd(2,a)*dl(j,a)
         F(2,3) = F(2,3) + Nxd(3,a)*dl(j,a)
         F(3,1) = F(3,1) + Nxd(1,a)*dl(k,a)
         F(3,2) = F(3,2) + Nxd(2,a)*dl(k,a)
         F(3,3) = F(3,3) + Nxd(3,a)*dl(k,a)

!        Calculate local wall property
         IF (useVarWall) eVWP(:) = eVWP(:) + Nd(a)*lVWP(:,a)

!        interpolate lagrange multiplier
         p_equi = p_equi + Np(a) * dl(4,a)
      END DO

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
!     Voigt notationa (Dm)
      CALL GETPK2CC(eq(cEq)%dmn(cDmn), F, nFn, fN, ya_g, grInt, S, Dm,
     2              eVWP, ifs, p_equi, stim)

!     1st Piola-Kirchhoff tensor (P)
      P    = MATMUL(F, S)

      DO a=1, eNoN
         Bm(1,1,a) = Nxd(1,a)*F(1,1)
         Bm(1,2,a) = Nxd(1,a)*F(2,1)
         Bm(1,3,a) = Nxd(1,a)*F(3,1)

         Bm(2,1,a) = Nxd(2,a)*F(1,2)
         Bm(2,2,a) = Nxd(2,a)*F(2,2)
         Bm(2,3,a) = Nxd(2,a)*F(3,2)

         Bm(3,1,a) = Nxd(3,a)*F(1,3)
         Bm(3,2,a) = Nxd(3,a)*F(2,3)
         Bm(3,3,a) = Nxd(3,a)*F(3,3)

         Bm(4,1,a) = (Nxd(1,a)*F(1,2) + F(1,1)*Nxd(2,a))
         Bm(4,2,a) = (Nxd(1,a)*F(2,2) + F(2,1)*Nxd(2,a))
         Bm(4,3,a) = (Nxd(1,a)*F(3,2) + F(3,1)*Nxd(2,a))

         Bm(5,1,a) = (Nxd(2,a)*F(1,3) + F(1,2)*Nxd(3,a))
         Bm(5,2,a) = (Nxd(2,a)*F(2,3) + F(2,2)*Nxd(3,a))
         Bm(5,3,a) = (Nxd(2,a)*F(3,3) + F(3,2)*Nxd(3,a))

         Bm(6,1,a) = (Nxd(3,a)*F(1,1) + F(1,3)*Nxd(1,a))
         Bm(6,2,a) = (Nxd(3,a)*F(2,1) + F(2,3)*Nxd(1,a))
         Bm(6,3,a) = (Nxd(3,a)*F(3,1) + F(3,3)*Nxd(1,a))
      END DO

!     Local residue and tangent matrices
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wd*(Nd(a)*ud(1) + Nxd(1,a)*P(1,1) +
     2      Nxd(2,a)*P(1,2) + Nxd(3,a)*P(1,3))
         lR(2,a) = lR(2,a) + wd*(Nd(a)*ud(2) + Nxd(1,a)*P(2,1) +
     2      Nxd(2,a)*P(2,2) + Nxd(3,a)*P(2,3))
         lR(3,a) = lR(3,a) + wd*(Nd(a)*ud(3) + Nxd(1,a)*P(3,1) +
     2      Nxd(2,a)*P(3,2) + Nxd(3,a)*P(3,3))

         DO b=1, eNoN
!           Geometric stiffness
            NxSNx = Nxd(1,a)*S(1,1)*Nxd(1,b) + 
     2              Nxd(2,a)*S(2,1)*Nxd(1,b) +
     3              Nxd(3,a)*S(3,1)*Nxd(1,b) + 
     4              Nxd(1,a)*S(1,2)*Nxd(2,b) +
     5              Nxd(2,a)*S(2,2)*Nxd(2,b) + 
     6              Nxd(3,a)*S(3,2)*Nxd(2,b) +
     7              Nxd(1,a)*S(1,3)*Nxd(3,b) + 
     8              Nxd(2,a)*S(2,3)*Nxd(3,b) +
     9              Nxd(3,a)*S(3,3)*Nxd(3,b)
            T1 = amd*Nd(a)*Nd(b) + afl*NxSNx

!           Material Stiffness (Bt*D*B)
            DBm = MATMUL(Dm, Bm(:,:,b))

            BmDBm = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     2              Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            lK(1,a,b) = lK(1,a,b) + wd*(T1 + afl*BmDBm)

            BmDBm = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2              Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     2              Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            lK(2,a,b) = lK(2,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2              Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     2              Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            lK(3,a,b) = lK(3,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2              Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     2              Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            lK(4,a,b) = lK(dof+1,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2              Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     2              Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            lK(5,a,b) = lK(dof+2,a,b) + wd*(T1 + afl*BmDBm)

            BmDBm = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2              Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     2              Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            lK(6,a,b) = lK(dof+3,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2              Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     2              Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            lK(7,a,b) = lK(2*dof+1,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2              Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     2              Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            lK(8,a,b) = lK(2*dof+2,a,b) + wd*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2              Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     2              Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            lK(9,a,b) = lK(2*dof+3,a,b) + wd*(T1 + afl*BmDBm)
         END DO
      END DO

      IF (eq(cEq)%dmn(cDmn)%phys .EQ. phys_gr) THEN
!        residual
         DO a=1, eNoN
            lR(4,a) = lR(4,a) + wp * Np(a) * (p_equi - stim)
         END DO
      END IF

      RETURN
      END SUBROUTINE STRUCT3D_GR