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
!     a:  input matrix
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
            IF ((ka2(j,i).ge.kmax).and.(ka2(j,i).le.tiny)) THEN
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
            IF ((ka2(j,i).ge.kmax).and.(ka2(j,i).le.tiny)) THEN
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
      SUBROUTINE cond_check(ps,vz,fp,fz,time,nn,nu,mm,hnu,nnv,nuv,mmv,hnuv,kup,omega)
!-----------------------------------------------------------------
!
! Condition check for the conservation of energy etc in HD 2D
!
      USE kes
      USE grid
      USE mpivars
      USE fft

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::ps,vz,fp,fz
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) ::C1
      DOUBLE PRECISION    :: Ep1,Ep2,Epv,Eph
      DOUBLE PRECISION    :: Ek1,Ek2
      DOUBLE PRECISION    :: Dsp,Hds,Dspw,Hdsw
      DOUBLE PRECISION    :: Dspv,Hdsv,Dsph1,Hdsh1,Dsph2,Hdsh2
      DOUBLE PRECISION    :: injp1,injp2,injz,injh1,injh2
      DOUBLE PRECISION    :: omega,coup
      DOUBLE PRECISION    :: nu,hnu,nuv,hnuv
      DOUBLE PRECISION    :: Efk, Efp, kup, kmn, Efpz
      DOUBLE PRECISION    :: tmq,tmp,tmp0,tmp1,tmp2,tmp3,tmp4,two,time
      INTEGER :: i,j,nn,mm,nnv,mmv


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
      SUBROUTINE rand_force(kdn,kup,fp0,dt,seed,kin,fp)
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
      DOUBLE PRECISION        :: dt,tmp1,tmp2,two,phase,kx,ky,theta

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
