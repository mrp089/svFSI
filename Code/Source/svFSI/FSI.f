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
!     This is for constructing FSI equations on fluid and solid
!     domains.
!
!--------------------------------------------------------------------

!     evaluates fsi in an element
      SUBROUTINE EVAL_FSI(e, lM, Ag, Yg, Dg, ptr, lK, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: ptr(lM%eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,lM%eNoN),
     2                                   lK(dof*dof,lM%eNoN,lM%eNoN)

      LOGICAL :: vmsStab
      INTEGER(KIND=IKIND) a, g, l, Ac, eNoN, cPhys, iFn, nFn,
     2                    iFluid, iSolid1, iSolid2, iInt
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd), grInt(nGrInt), diff, tol,
     2                 du, ii, jj, swss
      TYPE(fsType) :: fs(2)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:),
     3   lKd(:,:,:), lVWP(:,:), lRtau(:,:),
     4   lRtmp(:,:), lRp(:,:), lRtau_fd(:,:), lR0(:,:), lK0(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:), wss(:,:), wsse(:), dwss(:,:), sdwss(:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
      END IF

!     l = 3, if nsd==2 ; else 6;
      l = nsymd

      ALLOCATE(xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN),
     4   lKd(dof*nsd,eNoN,eNoN), lVWP(nvwp,eNoN), lRtau(dof,eNoN),
     5   lRtmp(dof,eNoN), lRtau_fd(dof,eNoN), lR0(dof,eNoN),
     6   lRp(dof,eNoN), lK0(dof*dof,eNoN,eNoN), wss(maxnsd,msh(1)%nNo),
     7   wsse(msh(1)%nEl), dwss(maxnsd,msh(1)%nNo), sdwss(maxnsd))

!     calculate wss for g&r
!     fixme: select fluid/solid mesh automatically
      IF (ALLOCATED(lM%grVn)) THEN
!        post-process wss
         wss = 0._RKIND
         wsse = 0._RKIND
         dwss = 0._RKIND
         CALL BPOST(msh(1), wss, wsse, dwss, Yg, Dg, outGrp_WSS)

!       is wss dependency correct?

!        fixme: run in parallel (evaluate for this element only)
!        todo: work with non-uniform meshes
!        map from fluid to solid interface
         DO iFluid=1, msh(1)%gnNo
            DO iSolid1=1, msh(2)%gnNo
!              check where fluid and solid nodes intersect
               IF (msh(1)%gN(iFluid) .EQ. msh(2)%gN(iSolid1)) THEN
!                 get wss norm on fluid interface
                  swss = SQRT(NORM(wss(:,iFluid)))
                  sdwss = dwss(:,iFluid)

!                 assign wss to solid interface
                  vWP0(7,iSolid1) = swss
                  vWP0(10:12,iSolid1) = sdwss(:)

!                 get solid interface id
                  iInt = vWP0(9,msh(2)%gN(iSolid1))

!                 assign wss to all points with same interface id
                  DO iSolid2=1, msh(2)%gnNo
                     IF (vWP0(9,msh(2)%gN(iSolid2)) .EQ. iInt) THEN
                        vWP0(7,msh(2)%gN(iSolid2)) = swss
                        vWP0(10:12,msh(2)%gN(iSolid2)) = sdwss
!                        WRITE(*,*) swss

!                       store wss of previously converged time step
                        IF (eq(cEq)%itr .EQ. 1) THEN
                          vWP0(13,msh(2)%gN(iSolid2)) = swss
!                          WRITE(*,*) vWP0(:,a)
                       END IF
                     END IF
                  END DO
               END IF
            END DO
         END DO
      END IF

         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

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

!        For FSI, fluid domain should be in the current configuration
         IF (cPhys .EQ. phys_fluid) THEN
            xl(:,:) = xl(:,:) + dl(nsd+2:2*nsd+1,:)
         END IF

!        Initialize residue and tangents
         lR  = 0._RKIND
         lK  = 0._RKIND
         lKd = 0._RKIND
         lK0 = 0._RKIND
!          lRp = 0._RKIND
!          lRtmp = 0._RKIND
!          lRtau_fd = 0._RKIND
!          lRtau = 0._RKIND

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN),
     2      Nwxx(l,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND
         Nqx      = 0._RKIND
         Nwxx     = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g),
     2            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
            END IF
            w = fs(1)%w(g) * Jac

!           retrieve g&r internal variables
            grInt(:) = 0._RKIND
            IF (ALLOCATED(lM%grVn)) grInt(1:nGrInt) = lM%grVo(:,g,e)

            IF (nsd .EQ. 3) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_lElas)
                  CALL LELAS3D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
     2               bfl, pS0l, pSl, lR, lK, lVWP)

               CASE (phys_struct)
                  CALL STRUCT3D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
     2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK,
     3               grInt, lVWP, lRtau)
!
!!                 check Jacobian with finite differences
!                  tol = 1.E-6_RKIND
!
!!                 perturb wss
!                  lVWP(7,:) = lVWP(7,:) + tol
!
!                  CALL STRUCT3D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
!     2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lRp, lK0,
!     3               grInt, lVWP, lRtmp)
!
!!                 restore wss
!                  lVWP(7,:) = lVWP(7,:) - tol

!               Update g&r variables
                IF (ALLOCATED(lM%grVn)) lM%grVo(:,g,e) = grInt(1:nGrInt)
                IF (ALLOCATED(lM%grVn)) lM%grVn(:,g,e) = grInt(1:nGrInt)

               CASE (phys_ustruct)
                  CALL USTRUCT3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
     2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
     3               dl, bfl, fN, ya_l, lR, lK, lKd)

               END SELECT

            ELSE IF (nsd .EQ. 2) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_lElas)
                  CALL LELAS2D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
     2               bfl, pS0l, pSl, lR, lK)

               CASE (phys_struct)
                  CALL STRUCT2D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
     2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
     2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
     3               dl, bfl, fN, ya_l, lR, lK, lKd)

               END SELECT
            END IF
         END DO ! g: loop
!
!!         calculate finite difference
!          lRtau_fd = (lRp - lR) / tol
!
!          WRITE(*,*) ""
!          WRITE(*,*) "analytical"
!          WRITE(*,*) lRtau(1:3,:)
!          WRITE(*,*) ""
!          WRITE(*,*) "FD"
!          WRITE(*,*) lRtau_fd(1:3,:)
!          WRITE(*,*) ""

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 2)

!        Gauss integration 2
         DO g=1, fs(2)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(2)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
     3               dl, bfl, lR, lK, lKd)
               END SELECT

            ELSE IF (nsd .EQ. 2) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
     3               dl, bfl, lR, lK, lKd)
               END SELECT
            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nwxx, Nqx)

      DEALLOCATE(xl, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lKd)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

      RETURN
      END SUBROUTINE EVAL_FSI
!####################################################################


!     evaluates fsi in an element with finite difference tangent matrix
      SUBROUTINE EVAL_FSI_FD(e, lM, Ag, Yg, Dg, ptr, lK, lR, lKfd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e
      REAL(KIND=RKIND), INTENT(INOUT) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: ptr(lM%eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,lM%eNoN),
     2  lK(dof*dof,lM%eNoN,lM%eNoN), lKfd(dof*dof,lM%eNoN,lM%eNoN)

      REAL(KIND=RKIND), ALLOCATABLE :: lRp(:,:), lRm(:,:),
     2  lKtmp(:,:,:), du, dlR(:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptrtmp(:)
      INTEGER(KIND=IKIND) ii, jj, kk, ll, Ac, dd
      REAL(KIND=RKIND) :: fa, fy, fd

      ALLOCATE(lRp(dof,lM%eNoN), lRm(dof,lM%eNoN),
     3         lKtmp(dof*dof,lM%eNoN,lM%eNoN),
     2         dlR(dof,lM%eNoN), ptrtmp(lM%eNoN))

!     time integration factors
      fd = eq(cEq)%af*eq(cEq)%beta*dt*dt
      fy = eq(cEq)%af*eq(cEq)%gam*dt
      fa = eq(cEq)%am

!     initialize
      lKfd = 0._RKIND

!     central evaluation
      CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptr, lK, lR)

!     loop nodes
      DO ii=1,dof
        DO jj=1,lM%eNoN
          du = 1.E-8_RKIND

!         global node id
          Ac = lM%IEN(jj,e)

!         perturb displacement vector
          Dg(ii,Ac) = Dg(ii,Ac) + du

          CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore displacement vector
          Dg(ii,Ac) = Dg(ii,Ac) - du

!         calculate finite difference
          dlR = fd * (lRp - lR) / du

!         perturb velocity vector
          Yg(ii,Ac) = Yg(ii,Ac) + du

          CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore velocity vector
          Yg(ii,Ac) = Yg(ii,Ac) - du

!         calculate finite difference
          dlR = dlR + fy * (lRp - lR) / du

!         perturb acceleration vector
          Ag(ii,Ac) = Ag(ii,Ac) + du

          CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptrtmp, lKtmp, lRp)

!         restore acceleration vector
          Ag(ii,Ac) = Ag(ii,Ac) - du

!         calculate finite difference
          dlR = dlR + fa * (lRp - lR) / du

!         assign to tangent matrix
          DO kk=1,dof
            dd = (kk - 1) * dof + ii
            DO ll=1,lM%eNoN
              lKfd(dd,ll,jj) = dlR(kk,ll)
            END DO
          END DO

        END DO
      END DO

      END SUBROUTINE EVAL_FSI_FD


!     evaluates elements and performs assembly
      SUBROUTINE CONSTRUCT_FSI(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      REAL(KIND=RKIND), ALLOCATABLE :: ptr(:), lR(:,:), lK(:,:,:),
     2                                 lKfd(:,:,:), diff(:,:,:)
      INTEGER(KIND=IKIND) e, eNoN, cPhys, ii, jj, kk, ll, dd

      eNoN = lM%eNoN

      ALLOCATE(ptr(eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     2         lKfd(dof*dof,eNoN,eNoN), diff(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
        cDmn  = DOMAIN(lM, cEq, e)
        cPhys = eq(cEq)%dmn(cDmn)%phys
        IF ((cPhys .NE. phys_fluid)  .AND.
     2      (cPhys .NE. phys_lElas)  .AND.
     3      (cPhys .NE. phys_struct) .AND.
     4      (cPhys .NE. phys_ustruct)) CYCLE

        IF (cPhys .EQ. phys_struct) THEN
!          CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptr, lK, lR)
          CALL EVAL_FSI_FD(e, lM, Ag, Yg, Dg, ptr, lK, lR, lKfd)
          lK = lKfd

!          DO ii=1,dof
!            DO jj=1,dof
!              DO kk=1,eNoN
!                DO ll=1,eNoN
!                  dd = (jj - 1) * dof + ii
!                  IF (.NOT.ISZERO(lK(dd,kk,ll))) THEN
!                    diff(dd,kk,ll) = ABS((lK(dd,kk,ll) - lKfd(dd,kk,ll))
!     2                      / lK(dd,kk,ll))
!                  ELSE
!                    diff(dd,kk,ll) = 0._RKIND
!                  END IF
!                END DO
!              END DO
!            END DO
!          END DO

!          WRITE(*,*) e
!          WRITE(*,*) MAXVAL(diff)

!          IF (e.EQ.1) THEN
!            WRITE(*,*) lKfd
!          END IF

!            IF (MAXVAL(diff) .GT. 1.E-3_RKIND) THEN
!              WRITE(*,*) e
!              DO kk=1,eNoN
!                DO ll=1,eNoN
!                  WRITE(*,*) kk,ll
!                    DO ii=1,3
!                      jj = (ii - 1) * 4
!!                      WRITE(*,*) lK(jj+1:jj+3,kk,ll)
!!                      WRITE(*,*) lKfd(jj+1:jj+3,kk,ll)
!                      WRITE(*,*) diff(jj+1:jj+3,kk,ll)
!                    END DO
!                    WRITE(*,*) ""
!                  END DO
!              END DO
!            END IF
        ELSE
          CALL EVAL_FSI(e, lM, Ag, Yg, Dg, ptr, lK, lR)
        END IF

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            IF (cPhys .EQ. phys_ustruct) err = "Cannot assemble "//
     2         "USTRUCT using Trilinos"
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
!            IF (cPhys .EQ. phys_ustruct) THEN
!               CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lK, lR)
!            ELSE
               CALL DOASSEM(eNoN, ptr, lK, lR)
!            END IF
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, lR, lK)

      END SUBROUTINE CONSTRUCT_FSI
