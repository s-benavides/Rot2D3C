!=================================================================
! PSEUDOSPECTRAL subroutines
!=================================================================

!*****************************************************************
      SUBROUTINE derivk2(a,b,dir)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =0 just copies a to b
!
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: dir
      INTEGER :: i,j

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax)) THEN
                  b(j,i) = im*ka(i)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE IF (dir.eq.2) THEN
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax)) THEN
                  b(j,i) = im*ka(j)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
      ELSE  ! copy
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax)) THEN
                  b(j,i) = a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO 
      ENDIF

      RETURN
      END SUBROUTINE derivk2

!*****************************************************************
      SUBROUTINE laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,n
            b(j,i) = -ka2(j,i)*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak2

!*****************************************************************
      SUBROUTINE inv_laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,n
            IF (ka2(j,i).gt.tiny) THEN
            b(j,i) = -a(j,i)/ka2(j,i)
            ELSE
            b(j,i) = 0.0d0
            ENDIF
         END DO
      END DO

      RETURN
      END SUBROUTINE inv_laplak2


!*****************************************************************
      SUBROUTINE curl_2D_z(a,b,c)
!-----------------------------------------------------------------
!
! Parameters
!     :  input matrix
!     b:  input matrix
!     c: output matrix
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE var
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      INTEGER :: i,j
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax)) THEN
               c(j,i) = im*( ka(i)*b(j,i) - ka(j)*a(j,i) )
            ELSE
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE curl_2D_z



!#################################################################
!#####################    NONLINEARITIES   #######################
!#################################################################

!*****************************************************************
      SUBROUTINE poisson(a,b,c)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: Poisson bracket {A,B} [output]
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO

!
! Computes dA/dy.dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/dble(n)**4
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).ge.kmax).or.(ka2(j,i).le.tiny)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE pmult(a,b,c)
!-----------------------------------------------------------------
!
! Pointwise multiplication of the scalar fields A and B
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: pointwise multiplication of a*b [output]
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

      CALL derivk2(a,c1,0)
      CALL derivk2(b,c2,0)

      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r1(i,j)*r2(i,j))/dble(n)**4
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).ge.kmax).or.(ka2(j,i).le.tiny)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO


      END SUBROUTINE pmult

!#################################################################
!########################   ANALYSIS   ###########################
!#################################################################

!*****************************************************************
      SUBROUTINE mom_calc(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the <|a|^kin >
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the moment
!     kin: =1 computes the first moment <|a|>
!          =2 computes the second moment <|a|^2>
!          =3 computes the third moment <|a|^3>,
!          etc.
!

      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION  :: b
      DOUBLE PRECISION  :: bloc
      DOUBLE PRECISION  :: tmp
      INTEGER :: kin,two
      INTEGER :: i,j

      bloc = 0.0d0
      tmp = 1.0d0/dble(n)**4

!
! Computes the kin'th moment of the scalar field
!
        DO i = ista,iend
           two = 2
           if (i.eq.1) two = 1
               DO j = 1,n
                  if (ka2(j,i).ge.tiny) then
                  bloc = bloc+two*abs(a(j,i))**kin*tmp
                  endif
               END DO
        END DO
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mom_calc

!###############################


!*****************************************************************
      SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2D,
! and the mean square current density or vorticity. 
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the energy
!     kin: =0 computes the square of the scalar field
!          =1 computes the energy
!          =2 computes the current or vorticity
!
      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
      DOUBLE PRECISION  :: b
      DOUBLE PRECISION  :: bloc
      DOUBLE PRECISION  :: tmp
      INTEGER :: kin,two
      INTEGER :: i,j

      bloc = 0.0d0
      tmp = 1.0d0/dble(n)**4

!
! Computes the square of the scalar field * k^(2kin)
!
        IF (kin.ge.0) THEN
            DO i = ista,iend
               two = 2
               if (i.eq.1) two = 1  
               DO j = 1,n
                  bloc = bloc+two*abs(a(j,i))**2*ka2(j,i)**kin*tmp
               END DO
            END DO
        ELSE
        DO i = ista,iend
               two = 2
               if (i.eq.1) two = 1
               DO j = 1,n
                  if (ka2(j,i).ge.tiny) then
                  bloc = bloc+two*abs(a(j,i))**2*ka2(j,i)**kin*tmp
                  endif
               END DO
            END DO
        ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!###############################
!*****************************************************************
      SUBROUTINE inerprod(a,b,kin,rslt)
!-----------------------------------------------------------------
! Parameters
!     a  : first  input matrix
!     b  : second input matrix
!     kin: = multiplies by the laplacian to this power

      USE kes
      USE grid
      USE mpivars
      USE ali

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION  :: tmp,tmq
      DOUBLE PRECISION  :: rslt
      INTEGER :: kin,two
      INTEGER :: i,j

      tmp = 0.0d0
      tmq = 1./dble(n)**4
      IF (kin.ge.0) THEN
         DO i = ista,iend
            two = 2
            if (i.eq.1) two = 1 
            DO j = 1,n
!            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
            IF ((ka2(j,i).le.kmax)) THEN
            tmp = tmp+two*(ka2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
       ELSE
         DO i = ista,iend
            two = 2
            if (i.eq.1) two = 1
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
!            IF ((ka2(j,i).le.kmax)) THEN
            tmp = tmp+two*(ka2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
       ENDIF
      CALL MPI_REDUCE(tmp,rslt,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rslt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  

      
      RETURN
      END SUBROUTINE inerprod
!###############################

!*****************************************************************
      SUBROUTINE restress_out(ps,vz,fp,fz,ps_out,vz_out,kin)
!-----------------------------------------------------------------
!
! Calculates Reynolds averages and stresses, depending on kin:
!      kin ==0 : computes <ps> and <vz> (average = y-average in space)    
!      kin ==1 : computes <fp> and <fz> (average = y-average in space)    
!      kin ==2 : computes <[psi',omega']> and <[psi',vz']> (average = y-average in space)    
!
      USE kes
      USE grid
      USE mpivars
      USE fft
      USE ali

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::ps,vz,fp,fz
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::ps_f,vz_f,fp_f,fz_f
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::C1,C2,C3
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)   ::ps_m,vz_m,fp_m,fz_m
      DOUBLE COMPLEX, DIMENSION(n,ista:iend),intent(out) ::ps_out,vz_out
      INTEGER :: i,j,kin

!!!!!!!!!!!!!!!!  MEAN - FLUCTUATIONS DECOMP !!!!!!!!!!!!!!!!
      ! First separate z-mean and fluctuations
      
      ! Calculate means
      do i = ista,iend
        do j = 1,n
          if (j.eq.1) then
            ! mean
            ps_m(j,i) = ps(j,i)
            vz_m(j,i) = vz(j,i)
            fp_m(j,i) = fp(j,i)
            fz_m(j,i) = fz(j,i)
            ! fluctuations
            ps_f(j,i) = 0.0d0
            vz_f(j,i) = 0.0d0
            fp_f(j,i) = 0.0d0
            fz_f(j,i) = 0.0d0
          else
            ! mean
            ps_m(j,i) = 0.0d0
            vz_m(j,i) = 0.0d0
            fp_m(j,i) = 0.0d0
            fz_m(j,i) = 0.0d0
            ! fluctuations
            ps_f(j,i) = ps(j,i)
            vz_f(j,i) = vz(j,i)
            fp_f(j,i) = fp(j,i)
            fz_f(j,i) = fz(j,i)
          endif
        end do
      end do

      if (kin.eq.0) then
          ps_out = ps_m
          vz_out = vz_m
      else if (kin.eq.1) then
          ps_out = fp_m
          vz_out = fz_m
      else if (kin.eq.2) then
          ! Calculate [psi',vz'] and other nonlinear terms
          CALL laplak2(ps_f,C1)               ! make - W_2D
          CALL poisson(ps_f,C1,C2)            ! -curl(u_2D x w_2D) (correct sign for lhs)
          CALL poisson(ps_f,vz_f,C3)          ! curl(u_2D x vz) (correct sign for rhs)

        do i = ista,iend
         do j = 1,n
           if (j.eq.1) then
           ! mean
            ps_out(j,i) = -C2(j,i)
            vz_out(j,i) = C3(j,i)
          else
            ps_out(j,i) = 0.0d0
            vz_out(j,i) = 0.0d0
          endif
         end do
        end do
      endif 

      RETURN
      END SUBROUTINE restress_out

!*****************************************************************
      SUBROUTINE cond_check(ps,vz,fp,fz,time,nn,nu,mm,hnu,nnv,nuv,mmv,hnuv,kup,omega,restress_calc)
!-----------------------------------------------------------------
!
! Condition check for the conservation of energy etc in HD 2D
!
      USE kes
      USE grid
      USE mpivars
      USE fft
      USE ali
      USE var

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::ps,vz,fp,fz
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::C1,C2,C3,C4
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: R1
      DOUBLE PRECISION    :: Ep1,Ep2,Epv,Eph
      DOUBLE PRECISION    :: Ek1,Ek2
      DOUBLE PRECISION    :: Dsp,Hds,Dspw,Hdsw
      DOUBLE PRECISION    :: Dspv,Hdsv,Dsph1,Hdsh1,Dsph2,Hdsh2
      DOUBLE PRECISION    :: injp1,injp2,injz,injh1,injh2
      DOUBLE PRECISION    :: omega,coup,maxv,minv,maxw,minw
      DOUBLE PRECISION    :: nu,hnu,nuv,hnuv
      DOUBLE PRECISION    :: b1,p1,b2,p2,pk,z1,zp1,z2,zp2
      DOUBLE PRECISION    :: Efk, Efp, kup, kmn, Efpz
      DOUBLE PRECISION    :: KE_ps_mean,KE_vz_mean,Inj_ps_mean,Inj_vz_mean 
      DOUBLE PRECISION    :: Diss_ps_mean, Diss_vz_mean, HDiss_ps_mean,HDiss_vz_mean
      DOUBLE PRECISION    :: NL_ps_mean,NL_vz_mean 
      DOUBLE PRECISION    :: tmq,tmp,tmp0,tmp1,tmp2,tmp3,tmp4,two,time
      DOUBLE COMPLEX      :: NL,z1_c,z2_c, z1_c_sum, z2_c_sum
      INTEGER :: i,j,nn,mm,nnv,mmv,cnt,cnt_tmp
      LOGICAL, INTENT(in) :: restress_calc

!! ENERGY
!       Horizontal 
      CALL energy(ps,Ep1,1)        ! |u|^2
      CALL inerprod(ps,fp,1,injp1) ! energy injection
      CALL energy(ps,Dsp,nn+1)     ! Dissipation 
      CALL energy(ps,Hds,1-mm)
!       vz
      CALL energy(vz,Epv,0)  ! |vz|^2
      CALL inerprod(vz,fz,0,injz) ! energy injection
      CALL energy(vz,Dspv,nnv)     ! Dissipation
      CALL energy(vz,Hdsv,-mmv)

!       extra: rotation coupling
      CALL derivk2(ps,C1,2)   ! vx
      CALL inerprod(vz,C1,0,coup)   ! coupling term

!! ENSTROPHY
      CALL energy(ps,Ep2,2)  ! |w|^2
      CALL inerprod(ps,fp,2,injp2) ! enstrophy injection
      CALL energy(ps,Dspw,nn+2)
      CALL energy(ps,Hdsw,2-mm)

!! HELICITY
      CALL inerprod(vz,ps,1,Eph) ! -vz*W_z
      CALL inerprod(vz,fp,1,injh1) ! helicity injection 1
      CALL inerprod(ps,fz,1,injh2) ! helicity injection 2
      CALL inerprod(vz,ps,nn+1,Dsph1) ! helicity dissipation 1: vz *nabla^2*nn (nabla^2 psi)
      CALL inerprod(vz,ps,1-mm,Hdsh1) ! helicity hypodissipation 1: vz 
      CALL inerprod(ps,vz,nnv+1,Dsph2) ! helicity dissipation 2
      CALL inerprod(ps,vz,1-mmv,Hdsh2) ! helicity hypodissipation 2

!! MEAN ENERGY BALANCE
      IF (restress_calc) THEN
          CALL restress_out(ps,vz,fp,fz,C1,C2,0) !ps_mean, vz_mean
          CALL energy(C1,KE_ps_mean,1)
          CALL energy(C2,KE_vz_mean,0)
          CALL restress_out(ps,vz,fp,fz,C3,C4,1) !fp_mean, fz_mean
          CALL inerprod(C1,C3,1,Inj_ps_mean) ! forcing injection
          CALL inerprod(C2,C4,0,inj_vz_mean) ! forcing injection
          CALL energy(C1,Diss_ps_mean,nn+1)     ! Dissipation 
          CALL energy(C2,Diss_vz_mean,nnv)     ! Dissipation 
          CALL energy(C1,HDiss_ps_mean,1-mm)
          CALL energy(C2,HDiss_vz_mean,-mmv)
          CALL restress_out(ps,vz,fp,fz,C3,C4,2) !<[psi',omega']>,<[psi',vz']>
          CALL inerprod(C1,C3,0,NL_ps_mean) ! NL injection (kin = 0 because [psi',omega'] is for omega, so for energy we just need psi*[psi',omega'])
          CALL inerprod(C2,C4,0,NL_vz_mean) ! NL injection
      ENDIF

!!!!!!!!!!! Computes the energy at largest scale!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(n)**4
      tmp1=0.0d0
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
!            IF ((ka2(j,i).le.(2.01)).and.(ka2(j,i).ge.0)) THEN
            IF ((kmn.gt.0).and.(kmn.le.1)) THEN
               tmp1 = tmp1+two*ka2(j,i)*abs(ps(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(tmp1,Ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!!!!
!!!!!!!!!!! Computes the energy at second-largest scale!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(n)**4
      tmp2=0.0d0
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
!            IF ((ka2(j,i).le.(4.01)).and.(ka2(j,i).ge.(2.01))) THEN
            IF ((kmn.gt.1).and.(kmn.le.2)) THEN
               tmp2 = tmp2+two*ka2(j,i)*abs(ps(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(tmp2,Ek2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!!!!


!!!!!!!!!!! Computes the energy at KUP scale!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(n)**4
      tmp1=0.0d0
      tmp2=0.0d0
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            IF ((ka2(j,i).le.(2.01)*kup*kup).and.(ka2(j,i).ge.kup*kup )) THEN
               tmp1 = tmp1+two*ka2(j,i)*abs(ps(j,i))**2*tmp
               tmp2 = tmp2+two*abs(vz(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(tmp1,Efp,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp2,Efpz,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)


!!!!!!!! COMPUTING MAX AND MIN VALUES !!!!!!!!!!!!
!!!!! FINDING MIN/MAX vorticity
      CALL laplak2(ps,C1)               ! make - W_2D
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      tmp=maxval(R1) !max vorticity
      tmp1=minval(R1) !min vorticity
      call MPI_REDUCE(tmp,maxw,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(tmp1,minw,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

!!!!! FINDING MIN/MAX vorticity
      CALL derivk2(vz,C1,0) ! Copies vz
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      tmp=maxval(R1) !max vz
      tmp1=minval(R1) !min vz
      call MPI_REDUCE(tmp,maxv,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(tmp1,minv,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

!!!!!! SYNCHRONIZATION ORDER PARAMETERS !!!!
!!!!!!  B1 and P1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(n)**4
         DO i = ista,iend

          DO j = 1,n
            ! Calculate classic Kuramoto parameters
            ! Only counting half of the wave numbers, because the
            ! condition of u being real means that z1 should be zero.
            ! Filter out small scales, as well.
            IF (((ka(j).ge.0).or.(i.gt.1)).and.(ka2(j,i).le.(70**2))) THEN
            cnt_tmp = cnt_tmp+1
            z1_c = z1_c + exp(im*atan2(aimag(ps(j,i)),real(ps(j,i))))
            z2_c = z2_c + exp(2*im*atan2(aimag(ps(j,i)),real(ps(j,i)))) 
            ENDIF
            
            ! This is just setting c1 to zero for use below.  
            c1(j,i) = 0.0d0 
          END DO
         END DO
      
     ! Sum over cores
      call MPI_REDUCE(z1_c,z1_c_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(z2_c,z2_c_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(cnt_tmp,cnt,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      !!if (myrank.eq.0) print*, 'DBG:   count = ',cnt, 'n^2/2',dble(n)**2/2, 1/(1.0d0/dble(n)/dble(n/2))

      ! Normalize, by total number of modes included 
      z1_c_sum = z1_c_sum/cnt
      z2_c_sum = z2_c_sum/cnt

      ! Take abs value and angle
      z1 = abs(z1_c_sum)
      zp1 = atan2(aimag(z1_c_sum),real(z1_c_sum))
      z2 = abs(z2_c_sum)
      zp2 = atan2(aimag(z2_c_sum),real(z2_c_sum))


      ! Now for energy based order parameters
         CALL laplak2(ps,c1)              ! make - W_z
         CALL poisson(ps,c1,c1)       ! -[psi,W_z]? or + ?
    
      IF (myrank.eq.0) THEN
          j = 6
          i = ista+5
          !!print*, '---- DBG ----,  kx = ',ka(j),'ky = ',ka(i)
          !!print*, '---- DBG ----,  abs(ps(j,i)) = ',abs(ps(j,i))
          NL = abs(ps(j,i))*c1(j,i)*tmp
          ! Order params
          
          b1 = abs(NL)
          p1 = atan2(aimag(NL),real(NL))
          pk = atan2(aimag(ps(j,i)),real(ps(j,i)))
      ENDIF   

!!!!!!  B2 and P2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = 0.0d0
               C2(j,i) = 0.0d0 
            END DO
         END DO

        ! Convert to chi
          DO i = ista,iend
             DO j = 1,n
                IF (abs(ps(j,i)).gt.tiny) THEN
                C2(j,i) = ps(j,i)*ps(j,i)/abs(ps(j,i))
                ELSE
                C2(j,i) = 0.0d0
                ENDIF
             END DO
          END DO

         CALL laplak2(C2,C1)          ! make lap(chi)
         CALL poisson(C2,C1,C1)       ! -[chi,lap(chi)]? or + ?

      IF (myrank.eq.0) THEN
          j = 6
          i = ista+5
          !!print*, '---- DBG ----,  kx = ',ka(j),'ky = ',ka(i)
          !!print*, '---- DBG ----,  abs(ps(j,i)) = ',abs(ps(j,i))
          NL = abs(c2(j,i))*c1(j,i)*tmp
          ! Order params
          
          b2 = abs(NL)
          p2 = atan2(aimag(NL),real(NL))
      ENDIF   
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy_bal.txt',position='append')
         WRITE(1,20) time,Ep1,Epv,injp1,injz,nu*Dsp,nuv*Dspv,hnu*Hds,hnuv*Hdsv,2*omega*coup,Efp,Efpz
   20    FORMAT(E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3, E23.14E3, E23.14E3)
         CLOSE(1)
         OPEN(1,file='enstrophy_bal.txt',position='append')
         WRITE(1,22) time,Ep2,injp2,nu*Dspw,hnu*Hdsw
   22    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3 ) 
         CLOSE(1)
         OPEN(1,file='helicity_bal.txt',position='append')
         WRITE(1,23) time,-Eph,-injh1,-injh2,-nu*Dsph1,-nuv*Dsph2,-hnu*Hdsh1,-hnuv*Hdsh2
   23    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3 ,E23.14E3,E23.14E3,E23.14E3)
         CLOSE(1)
         OPEN(1,file='max_min.txt',position='append')
         WRITE(1,24) time,maxw,minw,maxv,minv
   24    FORMAT( E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
         CLOSE(1)
         OPEN(1,file='sync.txt',position='append')
         WRITE(1,25) time,b1,p1,b2,p2,pk,z1,zp1,z2,zp2
   25    FORMAT(E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
         CLOSE(1)
         OPEN(1,file='mean.txt',position='append')
         WRITE(1,26) time,KE_ps_mean,KE_vz_mean,Inj_ps_mean,Inj_vz_mean,Diss_ps_mean, Diss_vz_mean, HDiss_ps_mean,HDiss_vz_mean, NL_ps_mean,NL_vz_mean
   26    FORMAT(E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3,E23.14E3)
         CLOSE(1)
      ENDIF
      
      RETURN
      END SUBROUTINE cond_check

!*****************************************************************
      SUBROUTINE spectrum(ps,vz,ext,odir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     b  : vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ezk,Ezktot1,Ektot1
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: ps,vz
      DOUBLE PRECISION        :: tmp,two,tmp1,tmp2,tmp3
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir

      tmp = 1.0d0/dble(n)**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the VELOCITY energy spectrum!!!!!!!!!!!!!!
      DO i = 1,n/2+1
         Ek(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+two*ka2(j,i)*abs(ps(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Ek,Ektot1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!! FOR DEBUGGING: spectrum of forcing !!!!!

!      DO i = 1,n/2+1
!         Fk(i) = 0.0d0
!      END DO
!      DO i = ista,iend
!         two=2.0d0
!         IF (i.eq.1) two=1.0d0
!         DO j = 1,n
!            kmn = int(sqrt(ka2(j,i))+.5d0)
!            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!               Fk(kmn) = Fk(kmn)+two*ka2(j,i)*abs(fp(j,i))**2*tmp
!            ENDIF
!         END DO
!      END DO
!      CALL MPI_REDUCE(Fk,Fktot1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)

!!!!!!! spectrum of vz !!!!!

      DO i = 1,n/2+1
         Ezk(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ezk(kmn) = Ezk(kmn)+two*abs(vz(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Ezk,Ezktot1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(odir) // '/spectrum.' // ext // '.txt')
!         WRITE(1,30) Q,0.0d0,0.5d0*tmp1
         do i=1,n/2+1
         WRITE(1,30) Ektot1(i),Ezktot1(i)
         enddo
         CLOSE(1)
      ENDIF
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  30    FORMAT( E24.15E3,E24.15E3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE spectrum

!*****************************************************************
      SUBROUTINE transfers(ps,vz,ext,odir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     b  : vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)        :: Fl
      DOUBLE PRECISION, DIMENSION(n/2+1)        :: Fl0,Fl1,Fl2
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: ps,vz
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: c1,c2
      DOUBLE PRECISION        :: tmp,two,tmp1,fx0,fx1,fx2
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir

      tmp = 1.0d0/dble(n)**4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = 0.0d0
            END DO
         END DO

         DO i = ista,iend
            DO j = 1,n
               C2(j,i) = 0.0d0
            END DO
         END DO


         CALL laplak2(ps,c1)               ! make - W_z
         CALL poisson(ps,c1,c2)            ! - curl(u x W_z) = [psi,W_z]
         CALL poisson(ps,vz,c1)            ! [psi,vz]

!!!!!!!!!!!!!!!!!!!!!!!!    Enstrophy flux 
      DO i = 1,n/2+1
         Fl(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            Fl(kmn) = Fl(kmn)+two*ka2(j,i)*dble(ps(j,i)*conjg(c2(j,i)))*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Fl,Fl0,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!    Energy flux: u_2D
      DO i = 1,n/2+1
         Fl(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            Fl(kmn) = Fl(kmn)+two*dble(ps(j,i)*conjg(c2(j,i)))*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Fl,Fl1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!    Energy flux: v_z
      DO i = 1,n/2+1
         Fl(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            Fl(kmn) = Fl(kmn)-two*dble(vz(j,i)*conjg(c1(j,i)))*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Fl,Fl2,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (myrank.eq.0) THEN
         fx0=0.0d0 
         fx1=0.0d0
         fx2=0.0d0 
        OPEN(1,file=trim(odir) // '/transfer.' // ext // '.txt')
         do i=1,n/2+1
         WRITE(1,40) Fl0(i), Fl1(i), Fl2(i)
         enddo
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/fluxes.' // ext // '.txt')
         do i=1,n/2+1
         fx0=Fl0(i)+fx0
         fx1=Fl1(i)+fx1
         fx2=Fl2(i)+fx2
         WRITE(1,40) fx0, fx1, fx2
         enddo
         CLOSE(1)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  40    FORMAT( E24.15E3, E24.15E3, E24.15E3 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE transfers

!*****************************************************************
      SUBROUTINE CFL_condition(cfl,ps,vz,nu,nn,nuv,nnv,omega,dt)
!-----------------------------------------------------------------
!        Parameters
!     cfl :cfl factor

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: ps
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: vz
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2,c3
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j,nn,nnv
      DOUBLE PRECISION        :: tmp,dt,cfl,nu,nuv,omega
      DOUBLE PRECISION        :: tmp1,tmp2,kcut,nrm,maxv

      kcut=(dble(n)/3.0d0)
       nrm=(dble(n))**2

!!!!! FINDING MAX(|nabla psi|) 
      CALL derivk2(ps,c1,1)
      CALL derivk2(ps,c2,2)    
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
         END DO
      END DO


      tmp=maxval(r3) !max energy density
      call MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
!!!!! FINDING MAX(VX)
      CALL derivk2(vz,c3,0) 
      CALL fftp2d_complex_to_real(plancr,c3,r1,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r2(i,j) = r1(i,j)*r1(i,j)
         END DO
      END DO
      tmp2 = maxval(r2) !max vz^2
      call MPI_REDUCE(tmp2,maxv,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!!!!!

      kcut=(dble(n)/3.0d0)
      tmp=sqrt(tmp1)/nrm  ! max horizontal velocity 
      tmp2 = sqrt(maxv)/nrm ! max vertical velocity
      dt = cfl/(kcut*tmp+nu*(kcut**(2*nn))+nuv*(kcut**(2*nnv))+2*omega)

 
      RETURN
      END SUBROUTINE CFL_condition

!*****************************************************************
      SUBROUTINE test_sub(time,ps,nu,nn,dt)
!-----------------------------------------------------------------
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: ps
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j,nn
      DOUBLE PRECISION        :: tmp,nu,dt,div,cfl,time
      DOUBLE PRECISION        :: tmp1,tmp2,kcut,nrm

      kcut=(dble(n)/3.0d0)
      nrm=(dble(n))**2
      tmp = 1/nrm

      CALL derivk2(ps,c1,1)     !D_x psi = -vy
      CALL derivk2(ps,c2,2)     !D_y psi =  vx

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Divergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL derivk2(c2,c3,1)     !c3 = D_x vx
      CALL derivk2(-c1,c4,2)     !c2 = D_y vy
      CALL fftp2d_complex_to_real(plancr,c3,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c4,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            tmp1 = (r1(i,j)+r2(i,j))*tmp
         END DO
      END DO
       
      CALL MPI_REDUCE(tmp1,div,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       CFL     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
                r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
         END DO
      END DO

      tmp=maxval(r3)
      call MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      kcut=(dble(n)/3.0d0)
      tmp=sqrt(tmp1)/nrm
      cfl = dt*(kcut*tmp+nu*(kcut**(2*nn)))

!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='test_sub.txt',position='append')
         WRITE(1,20) time, div, cfl
   20 FORMAT(E23.14E3,E23.14E3,E23.14E3 )
         CLOSE(1)

      ENDIF

      RETURN
      END SUBROUTINE test_sub


!*****************************************************************
      SUBROUTINE const_inj(ps,kdn,kup,fp0,fp,kin,seed1)
!-----------------------------------------------------------------
!       This subroutine assures that we inject constant energy.
!       It is called when iflow == 2
!       kin == 0 for vz
!           == 1 for ps

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE var
      USE random
      IMPLICIT NONE
!                                               ps
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: ps,fp
      INTEGER :: i,j,seed1,kin
      DOUBLE PRECISION        :: tmp,kdn,kup,Efp,fp0
      DOUBLE PRECISION        :: tmp1,tmp2,tmp3,two,phase2d


!!!!!!!!!!! 2D Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the energy at FORCING scale!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(n)**4
      tmp1= 0.0d0
      tmp2= 0.0d0

      DO i = ista,iend
         IF (i.ne.1) THEN
         two = 2.0d0
         DO j = 1,n
            IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn )) THEN
               fp(j,i)=  ps(j,i)/(abs(ps(j,i))**2+1.0d0)
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(n/2+1,i) = 0.0d0
         DO j = 2,n/2
            IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn ) ) THEN
               fp(    j,i) =  ps(j,i)/(abs(ps(j,i))**2+1.0d0)
               fp(n-j+2,i) = conjg(fp(j,i))
            ELSE
               fp(    j,i) = 0.0d0
               fp(n-j+2,i) = 0.0d0
            ENDIF
         END DO
         ENDIF
      END DO
        
      CALL inerprod(ps,fp,kin,Efp) ! Finds the dot product: either |nabla psi \nabla fpsi| or |vz fz| 
!!!!!Rescaling of forcing!!!
      seed1=myrank
      DO i = ista,iend
         IF (i.ne.1) THEN
         DO j = 1,n
            IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn )) THEN
               tmp3    = randu(seed1)*sqrt(ka2(j,i)) 
               fp(j,i) = fp(j,i)*fp0/Efp + im*tmp3*ps(j,i)
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(n/2+1,i) = 0.0d0
         DO j = 2,n/2
            IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn )) THEN
               tmp3    = randu(seed1)*sqrt(ka2(j,i))
               fp(    j,i) = fp(j,i)*fp0/Efp + im*tmp3*ps(j,i)
               fp(n-j+2,i) = conjg(fp(j,i))
            ELSE
               fp(    j,i) = 0.0
               fp(n-j+2,i) = 0.0
            ENDIF
         END DO
         ENDIF
      END DO
      
      RETURN
      END SUBROUTINE const_inj

!*****************************************************************
      SUBROUTINE rand_force(kdn,kup,fp0,dt,seed,kin,fp,phase,kx,ky)
!-----------------------------------------------------------------
!       This subroutine creates random forcing.
!       It is called when iflow == 3.
!       kin == 0   for vz   
!           == 1   for ps
      USE var
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE random
        IMPLICIT NONE
!                                              
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::fp
      INTEGER :: i,j,seed,kin
      DOUBLE PRECISION        :: tmp,kdn,kup,fp0,energyfp,energyfp2
      DOUBLE PRECISION        :: dt,tmp1,tmp2,two,theta
      DOUBLE PRECISION, intent(out) :: phase
      INTEGER, intent(out)    :: kx,ky

!!!!!!!!!!! Like Chan et al. 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Choose random vector of length kup and random phase !!!
        IF (myrank.eq.0) THEN
                 theta = randu(seed)*pi
                 phase=randu(seed+1)*pi
! Aug 28                 phase = 0.2342987*pi
                 kx = floor(kup*cos(theta)+0.5)
                 ky = floor(kup*sin(theta)+0.5)
!                 print*,"DBG (kx,ky)", kx,ky
         ENDIF
         CALL MPI_BCAST(kx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(ky,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(phase,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp = sqrt(2.0d0*fp0)/(sqrt(dt)*kup**kin) ! kup**kin to distinguish between vz and ps
         tmp1 = tmp*cos(phase)
         tmp2 = tmp*sin(phase)
      DO i = ista,iend
         IF (i.ne.1) THEN
         DO j = 1,n
            IF (((ka(i).eq.ky).or.(ka(i).eq.(-ky))).and.(ka(j).eq.kx)) THEN
            fp(j,i) = (tmp1+im*tmp2)*dble(n)**2
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(n/2+1,i) = 0.0d0
         DO j = 2,n/2
            IF ((ka(i).eq.ky).and.(ka(j).eq.kx)) THEN
               fp(j,i) = (tmp1+im*tmp2)*dble(n)**2
               fp(n-j+2,i) = conjg(fp(j,i))
            ELSE
               fp(j,i) = 0.0d0
               fp(n-j+2,i) = 0.0d0
            ENDIF
         END DO
         ENDIF
      END DO
      

      RETURN
      END SUBROUTINE rand_force

!*****************************************************************
      SUBROUTINE sync_shell(ps,fp,nu,hnu,nn,mm,phasefp,kxfp,kyfp,ext,odir)
!-----------------------------------------------------------------
!
! Computes the synchronization between nodes. Averages over wavenumber
! shells. Uses FFT to make this computation. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     b  : vector potential
!     ext: the extension used when writting the file
!     odir: output directory
!
      USE kes
      USE ali
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1)        :: B,P,PKt,trans,hdiss,diss,en,force
      DOUBLE PRECISION, DIMENSION(n/2+1)        :: B1,P1,B2,P2,cosk,sink,T,D,H,E,F,Z1_a,Z2_a,Zp1,Zp2
      DOUBLE PRECISION, DIMENSION(n/2+1)        :: NORM,PKt2,PKt3,Pkt4,cos2k,sin2k,phipsi_R,phipsi_I
      INTEGER,          DIMENSION(n/2+1)        :: CNT,CNT_SUM
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: ps, fp
      DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: c1,c2,vl,vs,c3,c4,c5,c6
      DOUBLE PRECISION        :: tmp,two,tmp1,nu,hnu,phik
      DOUBLE COMPLEX          :: NL,NLll,NLss,NLls,NLsl
      INTEGER     :: kmn,nn,mm,kmnf
      INTEGER     :: i,j
      INTEGER, intent(in) :: kxfp,kyfp
      DOUBLE PRECISION, intent(in) :: phasefp
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir

      ! For list making
      double precision, dimension (:), allocatable :: kxs_kf,kys_kf,phi_kf,psi_kf
      double precision, dimension (:), allocatable :: kxs_k6,kys_k6,phi_k6,psi_k6
      double precision, dimension (:), allocatable :: kxs_k25,kys_k25,phi_k25,psi_k25
      double precision, dimension (:), allocatable :: psi_ll_k6,psi_sl_k6,psi_ls_k6,psi_ss_k6
      double precision, dimension (:), allocatable :: B_ll_k6,B_sl_k6,B_ls_k6,B_ss_k6
      INTEGER       :: ic,id,iu,iif,ii6,ii25
      CHARACTER     :: c,y,u
      CHARACTER*3   :: node

      ic = 48+int(myrank/100)
      id = 48+int(myrank/10)-int(myrank/100)*10
      iu = 48+int(myrank)-int(myrank/10)*10
      c = char(ic)
      y = char(id)
      u = char(iu)
      node = c // y // u

      tmp = 1.0d0/dble(n)**4


!!!!!!!!!!!!!!!!!!!!!!!!    
      DO i = 1,n/2+1
         CNT(i) = 0
         B(i) = 0.0d0
         P(i) = 0.0d0
         PKt(i) = 0.0d0
         PKt2(i) = 0.0d0
         PKt3(i) = 0.0d0
         PKt4(i) = 0.0d0
         trans(i) = 0.0d0
         diss(i) = 0.0d0
         hdiss(i) = 0.0d0
         en(i) = 0.0d0
         force(i) = 0.0d0
         Z1_a(i) = 0.0d0
         Z2_a(i) = 0.0d0
         Zp1(i) = 0.0d0
         Zp2(i) = 0.0d0
         NORM(i) = 0.0d0
      END DO

!!!!!!  B2 and P2, Count
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Convert to chi
          DO i = ista,iend
             DO j = 1,n
                IF (abs(ps(j,i)).gt.tiny) THEN
                c2(j,i) = ps(j,i)*ps(j,i)/abs(ps(j,i))
                ELSE
                c2(j,i) = 0.0d0
                ENDIF
             END DO
          END DO

         CALL laplak2(c2,c1)               ! make lap(chi)
         CALL poisson(c2,c1,c1)            ! [chi,lap(chi)]

!!!!!!!!!!!!!!!!!!!!!!!!    B_2 and Psi_2 
      DO i = 1,n/2+1
         B(i) = 0.0d0
         P(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
           kmn = int(sqrt(ka2(j,i))+.5d0)
           IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            NL = -abs(c2(j,i))*c1(j,i)
            ! Order params. Only integrate over half shell, so that
            ! it's not purely a real value.                 
            IF ((ka(j).ge.0).or.(i.gt.1)) THEN
            CNT(kmn) = CNT(kmn) + 1
            B(kmn) = B(kmn) +  abs(NL)*tmp
            P(kmn) = P(kmn) + atan2(aimag(NL),real(NL)) ! to be Shell averaged
            ENDIF
           ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(CNT,CNT_SUM,n/2+1,MPI_INTEGER,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(B,B2,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(P,P2,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ! Divide by count, so that it's a shell average:
      DO i = 1,n/2+1
        P2(i) = MODULO(P2(i)/dble(CNT_SUM(i))+pi,2*pi) - pi
      END Do

!!!!!!  B1 and P1 and more
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = ista,iend
            DO j = 1,n
               c1(j,i) = 0.0d0 
               c2(j,i) = 0.0d0 
               ! Apply filtering
               kmn = int(sqrt(ka2(j,i))+.5d0)
               IF (kmn.gt.6) THEN 
               vs(j,i) = ps(j,i)
               vl(j,i) = 0.0d0 + im*0.0d0
               ELSE
               vl(j,i) = ps(j,i)
               vs(j,i) = 0.0d0 + im*0.0d0
               END IF
            END DO
         END DO

         ! Regular 
         CALL laplak2(ps,c1)               ! make - W_z
         CALL poisson(ps,c1,c2)            ! -curl(u_2D x w_2D) = u.grad(u)
         ! <,<
         CALL laplak2(vl,c1)               ! make - W_z
         CALL poisson(vl,c1,c3)            ! -curl(u_2D x w_2D)
         ! <,>
         CALL laplak2(vs,c1)               ! make - W_z
         CALL poisson(vl,c1,c4)            ! -curl(u_2D x w_2D)
         ! >,<
         CALL laplak2(vl,c1)               ! make - W_z
         CALL poisson(vs,c1,c5)            ! -curl(u_2D x w_2D)
         ! >,>
         CALL laplak2(vs,c1)               ! make - W_z
         CALL poisson(vs,c1,c6)            ! -curl(u_2D x w_2D)

!!!! For lists. Notice that each core will save a separate list.
! First find forcing scale
      kmnf = int(sqrt(dble(kxfp)**2+dble(kyfp)**2)+.5d0)
      allocate(kxs_kf(CNT(kmnf))) 
      allocate(kys_kf(CNT(kmnf)))
      allocate(psi_kf(CNT(kmnf)))
      allocate(phi_kf(CNT(kmnf)))
      allocate(kxs_k6(CNT(6)))
      allocate(kys_k6(CNT(6)))
      allocate(psi_k6(CNT(6)))
      allocate(psi_ll_k6(CNT(6)))
      allocate(psi_ls_k6(CNT(6)))
      allocate(psi_sl_k6(CNT(6)))
      allocate(psi_ss_k6(CNT(6)))
      allocate(B_ll_k6(CNT(6)))
      allocate(B_ls_k6(CNT(6)))
      allocate(B_sl_k6(CNT(6)))
      allocate(B_ss_k6(CNT(6)))
      allocate(phi_k6(CNT(6)))
      allocate(kxs_k25(CNT(25)))
      allocate(kys_k25(CNT(25)))
      allocate(psi_k25(CNT(25)))
      allocate(phi_k25(CNT(25)))
      iif = 1
      ii6 = 1
      ii25 = 1


      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
                kmn = int(sqrt(ka2(j,i))+.5d0)
                IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                phik = atan2(aimag(ps(j,i)),real(ps(j,i)))
                NL = -abs(ps(j,i))*c2(j,i)
                ! Filtered terms
                NLll = -abs(ps(j,i))*c3(j,i)
                NLsl = -abs(ps(j,i))*c4(j,i)
                NLls = -abs(ps(j,i))*c5(j,i)
                NLss = -abs(ps(j,i))*c6(j,i)

                ! Order params. Only integrate over half shell, so that
                ! it's not purely a real value.                 
                IF ((ka(j).ge.0).or.(i.gt.1)) THEN
                B(kmn) = B(kmn) +  abs(NL)*tmp ! Shell integral
                P(kmn) = P(kmn) + atan2(aimag(NL),real(NL)) ! to be Shell averaged
                PKt(kmn) = PKt(kmn) +  cos(atan2(aimag(NL),real(NL))-phik) !to be Shell averaged
                PKt2(kmn) = PKt2(kmn) +  sin(atan2(aimag(NL),real(NL))-phik) !to be Shell averaged
                PKt3(kmn) = PKt3(kmn) +  cos(2.0d0*atan2(aimag(NL),real(NL))-2.0d0*phik) !to be Shell averaged
                PKt4(kmn) = PKt4(kmn) +  sin(2.0d0*atan2(aimag(NL),real(NL))-2.0d0*phik) !to be Shell averaged
                  !! Making lists to output.
                  !IF (((ka(i).eq.kyfp).or.(ka(i).eq.(-kyfp))).and.(ka(j).eq.kxfp)) THEN
                  !ENDIF 
                  IF (kmn.eq.kmnf) THEN 
                    kxs_kf(iif) = ka(j)
                    kys_kf(iif) = ka(i)
                    phi_kf(iif) = phik 
                    psi_kf(iif) = atan2(aimag(NL),real(NL))
                    iif = iif + 1
                  ELSEIF (kmn.eq.6) THEN 
                    kxs_k6(ii6) = ka(j)
                    kys_k6(ii6) = ka(i)
                    phi_k6(ii6) = phik 
                    psi_k6(ii6) = atan2(aimag(NL),real(NL))
                    psi_ll_k6(ii6) = atan2(aimag(NLll),real(NLll))
                    psi_sl_k6(ii6) = atan2(aimag(NLsl),real(NLsl))
                    psi_ls_k6(ii6) = atan2(aimag(NLls),real(NLls))
                    psi_ss_k6(ii6) = atan2(aimag(NLss),real(NLss))
                    B_ll_k6(ii6) = abs(NLll)
                    B_sl_k6(ii6) = abs(NLsl)
                    B_ls_k6(ii6) = abs(NLls)
                    B_ss_k6(ii6) = abs(NLss)
                    ii6 = ii6 + 1
                  ELSEIF (kmn.eq.25) THEN 
                    kxs_k25(ii25) = ka(j)
                    kys_k25(ii25) = ka(i)
                    phi_k25(ii25) = phik 
                    psi_k25(ii25) = atan2(aimag(NL),real(NL))
                    ii25 = ii25 + 1
                  ENDIF
                ENDIF

                ! Energy budget
                trans(kmn) = trans(kmn) + two*real(exp(-im*phik)*NL)*tmp ! = B_1 cos(phi_k-Psi)
                diss(kmn) = diss(kmn) - two*nu*ka2(j,i)**(1+nn)*abs(ps(j,i))**2*tmp 
                hdiss(kmn) = hdiss(kmn) - two*hnu*ka2(j,i)**(1-mm)*abs(ps(j,i))**2*tmp 
                en(kmn) = en(kmn) + two*ka2(j,i)*abs(ps(j,i))**2*tmp ! e(k) = k^2 * abs(psi)^2
                force(kmn) = force(kmn) + 0.5*two*real(ka2(j,i)*abs(ps(j,i))*exp(-im*phik)*fp(j,i))*tmp ! e(k) = k^2 * abs(psi)^2
                ! The 0.5 comes from random forcing.
                ENDIF   
         END DO
      END DO
      CALL MPI_REDUCE(B,B1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(P,P1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(PKt,cosk,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(PKt2,sink,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(PKt3,cos2k,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(PKt4,sin2k,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trans,T,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(diss,D,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(hdiss,H,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(en,E,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(force,F,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ! Divide by count, so that it's a shell average:
      DO i = 1,n/2+1
        P1(i) = MODULO(P1(i)/dble(CNT_SUM(i))+pi,2*pi) - pi
        cosk(i) = cosk(i)/dble(CNT_SUM(i)) ! cos(phi_k - Psi_1(k))
        sink(i) = sink(i)/dble(CNT_SUM(i)) ! sins(phi_k - Psi_1(k))
        cos2k(i) = cos2k(i)/dble(CNT_SUM(i)) ! cos(2phi_k - 2Psi_1(k))
        sin2k(i) = sin2k(i)/dble(CNT_SUM(i)) ! sin(2phi_k - 2Psi_1(k))
      END Do
      

      IF (size(phi_kf).gt.0) THEN
          OPEN(1,file=trim(odir)//'/ph_kf.'//node//'.'//ext//'.txt')
             do i=1,size(phi_kf)
             WRITE(1,45) kxs_kf(i),kys_kf(i), phi_kf(i),psi_kf(i)
             enddo
          CLOSE(1)
      ENDIF
      IF (size(phi_k6).gt.0) THEN
      OPEN(1,file=trim(odir)//'/ph_k6.'//node//'.'//ext//'.txt')
         do i=1,size(phi_k6)
         WRITE(1,47) kxs_k6(i),kys_k6(i),phi_k6(i),psi_k6(i),psi_ll_k6(i),psi_sl_k6(i),psi_ls_k6(i),psi_ss_k6(i), B_ll_k6(i),B_sl_k6(i),B_ls_k6(i),B_ss_k6(i)
         enddo
      CLOSE(1)
      ENDIF
      IF (size(phi_k25).gt.0) THEN
      OPEN(1,file=trim(odir)//'/ph_k25.'//node//'.'//ext//'.txt')
         do i=1,size(phi_k25)
         WRITE(1,45) kxs_k25(i),kys_k25(i),phi_k25(i),psi_k25(i)
         enddo
      CLOSE(1)
      ENDIF
  45    FORMAT( E24.15E3 ,  E24.15E3 , E24.15E3,  E24.15E3 )
  47    FORMAT( E24.15E3 ,  E24.15E3 , E24.15E3,  E24.15E3, E24.15E3, E24.15E3 , E24.15E3 , E24.15E3, E24.15E3, E24.15E3 , E24.15E3 , E24.15E3 )
      IF (myrank.eq.0) THEN
      OPEN(1,file=trim(odir)//'/ph_f.'//ext//'.txt')
         WRITE(1,46) dble(kxfp),dble(kyfp), phasefp
      CLOSE(1)
      ENDIF
  46    FORMAT( E24.15E3 ,  E24.15E3 , E24.15E3 )


!!!!!!  Z1 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Just keep the phase.
          DO i = ista,iend
             DO j = 1,n
                IF (abs(ps(j,i)).gt.tiny) THEN
                c2(j,i) = ps(j,i)/abs(ps(j,i))
                ELSE
                c2(j,i) = 0.0d0
                ENDIF
             END DO
          END DO

         CALL laplak2(c2,c1)               ! make lap(chi)
         CALL poisson(c2,c1,c1)            ! [chi,lap(chi)]? or + ?

!!!!!!!!!!!!!!!!!!!!!!!!  
      DO i = 1,n/2+1
         B(i) = 0.0d0
         P(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            NL = abs(c2(j,i))*c1(j,i)
            ! Order params. Only integrate over half shell, so that
            ! it's not purely a real value.                 
            IF ((ka(j).ge.0).or.(i.gt.1)) THEN
            B(kmn) = B(kmn) +  abs(NL)
            P(kmn) = P(kmn) + atan2(aimag(NL),real(NL)) ! to be Shell averaged
            ENDIF

            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(B,Z1_a,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(P,Zp1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ! Divide by count, so that it's a shell average:
      DO i = 1,n/2+1
        Z1_a(i) = Z1_a(i)/dble(CNT_SUM(i))
        Zp1(i) = MODULO(Zp1(i)/dble(CNT_SUM(i))+pi,2*pi) - pi
      END Do

!!!!!!  Z2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Just keep the phase.
          DO i = ista,iend
             DO j = 1,n
                IF (abs(ps(j,i)).gt.tiny) THEN
                c2(j,i) = ps(j,i)*ps(j,i)/abs(ps(j,i))**2
                ELSE
                c2(j,i) = 0.0d0
                ENDIF
             END DO
          END DO

         CALL laplak2(c2,c1)               ! make lap(chi)
         CALL poisson(c2,c1,c1)            ! [chi,lap(chi)]? or + ?

!!!!!!!!!!!!!!!!!!!!!!!!  
      DO i = 1,n/2+1
         B(i) = 0.0d0
         P(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            NL = abs(c2(j,i))*c1(j,i)
            ! Order params. Only integrate over half shell, so that
            ! it's not purely a real value.                 
            IF ((ka(j).ge.0).or.(i.gt.1)) THEN
            B(kmn) = B(kmn) +  abs(NL)
            P(kmn) = P(kmn) + atan2(aimag(NL),real(NL)) ! to be Shell averaged
            ENDIF

            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(B,Z2_a,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(P,Zp2,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ! Divide by count, so that it's a shell average:
      DO i = 1,n/2+1
        Z2_a(i) = Z2_a(i)/dble(CNT_SUM(i))
        Zp2(i) = MODULO(Zp2(i)/dble(CNT_SUM(i))+pi,2*pi) - pi
      END Do

!!!!!!  NORM 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = ista,iend
            DO j = 1,n
               c1(j,i) = 0.0d0
               c2(j,i) = 1.0d0 
            END DO
         END DO

         CALL laplak2(c2,c1)               ! make lap(chi)
         CALL poisson(c2,c1,c1)            ! [chi,lap(chi)]? or + ?

!!!!!!!!!!!!!!!!!!!!!!!!  
      DO i = 1,n/2+1
         B(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            NL = abs(c2(j,i))*c1(j,i)
            ! Order params. Only integrate over half shell, so that
            ! it's not purely a real value.                 
            IF ((ka(j).ge.0).or.(i.gt.1)) THEN
            B(kmn) = B(kmn) +  abs(NL)
            ENDIF

            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(B,NORM,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      ! Divide by count, so that it's a shell average:
      DO i = 1,n/2+1
        NORM(i) = NORM(i)/dble(CNT_SUM(i))
      END Do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (myrank.eq.0) THEN
        OPEN(1,file=trim(odir) // '/sync.' // ext // '.txt')
         do i=1,n/2+1
         WRITE(1,40) B1(i), P1(i), B2(i), P2(i), cosk(i),sink(i), T(i), D(i), H(i), E(i), F(i), Z1_a(i),Z2_a(i),Zp1(i),Zp2(i), NORM(i), cos2k(i),sin2k(i)
         enddo
         CLOSE(1)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  40    FORMAT( E24.15E3, E24.15E3,E24.15E3,E24.15E3,E24.15E3, E24.15E3, E24.15E3, E24.15E3 , E24.15E3 , E24.15E3, E24.15E3 , E24.15E3, E24.15E3, E24.15E3, E24.15E3, E24.15E3, E24.15E3, E24.15E3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE sync_shell

