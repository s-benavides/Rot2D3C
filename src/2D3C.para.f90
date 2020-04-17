!=================================================================
      PROGRAM 2D3C
!=================================================================
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!=================================================================

      USE mpivars
      USE fft
      USE ali
      USE var
      USE kes
      USE grid
      USE random
      IMPLICIT NONE

!
! Integration parameters
!     ord  : order of the Runge-Kutta method used

      INTEGER, PARAMETER :: ord = 4
      INTEGER :: ini
      INTEGER :: step
      INTEGER :: tstep
      INTEGER :: cstep
      INTEGER :: sstep

!
! streamfunction, vector potential, z component 
! of the fields and external force matrixes


      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: ps   ! 2D streamfunction
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: vz   ! v_x 
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fp
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fz
 
!
! Temporal data storing matrixes

      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C1
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C2
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C3
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C4
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C5
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C6

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: R1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: R2
!
! Some auxiliary matrixes

      DOUBLE PRECISION :: ener
      DOUBLE PRECISION :: enerv,enst
      DOUBLE PRECISION :: eneru,jnst
      DOUBLE PRECISION :: dt,dt_new,CFL,dt_corr,timecorr
      DOUBLE PRECISION :: kup,kdn,kr
      DOUBLE PRECISION :: kmup,kmdn
      DOUBLE PRECISION :: prm1,prm2
      DOUBLE PRECISION :: dump,tmp
      DOUBLE PRECISION :: tmp1,tmp2,tmp3,tmp4,tmp5
      DOUBLE PRECISION :: fp0,u0
      DOUBLE PRECISION :: fz0,v0
      DOUBLE PRECISION :: time,omega
      DOUBLE PRECISION :: nu,hnu
      DOUBLE PRECISION :: nuv,hnuv
      DOUBLE PRECISION :: phase1,phase2
      DOUBLE PRECISION :: phase3,phase4
      DOUBLE PRECISION :: cphi,dphi,fphi

      INTEGER :: stat
      INTEGER :: t,o,nn,mm,nnv,mmv
      INTEGER :: i,j,ir,jr
      INTEGER :: ki,kj
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju,jt
      INTEGER :: timet,timec,times,timec2
      INTEGER :: seed,iflow, seed1

      CHARACTER     :: c,d,u,th
      CHARACTER*3   :: node,ext
      CHARACTER*4   :: ext4
      CHARACTER*100 :: odir
      CHARACTER*100 :: idir

!
! Initializes the MPI library

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      ic = 48+int(myrank/100)
      id = 48+int(myrank/10)-int(myrank/100)*10
      iu = 48+int(myrank)-int(myrank/10)*10
      c = char(ic)
      d = char(id)
      u = char(iu)
      node = c // d // u

!
! Allocates memory for distributed blocks

      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,jsta,jend)

      ALLOCATE( R1(n,jsta:jend) )
      ALLOCATE( R2(n,jsta:jend) )
      ALLOCATE( C1(n,ista:iend) )
      ALLOCATE( C2(n,ista:iend) )
      ALLOCATE( C3(n,ista:iend) )
      ALLOCATE( C4(n,ista:iend) )
      ALLOCATE( C5(n,ista:iend) )
      ALLOCATE( C6(n,ista:iend) )
      ALLOCATE( ps(n,ista:iend) )
      ALLOCATE( vz(n,ista:iend) )
      ALLOCATE( fp(n,ista:iend) )
      ALLOCATE( fz(n,ista:iend) )
      ALLOCATE( ka(n), ka2(n,ista:iend) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Reads from the external file 'status.txt'
! the status of a previous run (if any)
!     stat: last output of a previous run

      IF (myrank.eq.0) THEN
         OPEN(1,file='status.inp',status='unknown')
         READ(1,*) stat
         READ(1,*) time
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown')
         READ(1,*) CFL                      ! 1
         READ(1,*) step                     ! 2
         READ(1,*) cstep                    ! 3
         READ(1,*) sstep                    ! 4
         READ(1,*) tstep                    ! 5
         READ(1,*) fp0                      ! 6
         READ(1,*) fz0                      ! 7
         READ(1,*) u0                       ! 8
         READ(1,*) v0                       ! 9
         READ(1,*) kdn                      ! 10
         READ(1,*) kup                      ! 11
         READ(1,*) nu                       ! 12
         READ(1,*) hnu                      ! 13
         READ(1,*) nn                       ! 14
         READ(1,*) mm                       ! 15
         READ(1,*) nuv                      ! 16
         READ(1,*) hnuv                     ! 17
         READ(1,*) nnv                      ! 18
         READ(1,*) mmv                      ! 19
         READ(1,*) omega                    ! 20
         READ(1,*) seed                     ! 21
         READ(1,*) iflow                    ! 22
         READ(1,*) prm1                     ! 23
         READ(1,*) prm2                     ! 24
         READ(1,*) dt_corr                  ! 25 
         READ(1,'(a100)') idir              ! binary input directory
         READ(1,'(a100)') odir              ! output directory
         CLOSE(1)
!         step = step
!         tstep = tstep
!         sstep = sstep
!         cstep = cstep
      ENDIF
      CALL MPI_BCAST(  CFL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 1
      CALL MPI_BCAST( step,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 2
      CALL MPI_BCAST(cstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 3
      CALL MPI_BCAST(sstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 4
      CALL MPI_BCAST(tstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 5
      CALL MPI_BCAST(  fp0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 6
      CALL MPI_BCAST(  fz0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 7
      CALL MPI_BCAST(   u0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 8
      CALL MPI_BCAST(   v0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 9
      CALL MPI_BCAST(  kdn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 10
      CALL MPI_BCAST(  kup,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 11
      CALL MPI_BCAST(   nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 12
      CALL MPI_BCAST(  hnu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 13
      CALL MPI_BCAST(   nn,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 14
      CALL MPI_BCAST(   mm,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 15
      CALL MPI_BCAST(  nuv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 16
      CALL MPI_BCAST( hnuv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 17
      CALL MPI_BCAST(  nnv,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 18
      CALL MPI_BCAST(  mmv,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 19
      CALL MPI_BCAST(omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 20
      CALL MPI_BCAST( seed,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 21
      CALL MPI_BCAST(iflow,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 22
      CALL MPI_BCAST( prm1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 23
      CALL MPI_BCAST( prm2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 24
      CALL MPI_BCAST(dt_corr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !25
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!
! Some numerical constants

      ic = 48
      id = 48
      iu = 48
      jt = 48
      jc = 48
      jd = 48
      ju = 48

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (dble(n)/3.d0)**2
      tiny =  0.000001d0

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         ka(i) = dble(i-1)
         ka(i+n/2) = dble(i-n/2-1)
      END DO
      DO i = ista,iend
         DO j = 1,n
            ka2(j,i) = ka(i)**2+ka(j)**2
         END DO
      END DO

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE 
! in long runs

      CALL fftp2d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp2d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)

!
! Sets the initial conditions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (stat.eq.0) THEN
         ini = 1
         timet = tstep
         timec = cstep
         timec2 = cstep
         times = sstep
!         timecorr = dt_corr

!STREAM FUNCTION R1 & PHI R2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO j = jsta,jend
            DO i = 1,n
               R1(i,j) = 0.0d0
               R2(i,j) = 0.0d0
            END DO
         END DO
         DO ir = 1,kup+1
         DO jr = 1,kup+1
            kr = sqrt(dble(ir*ir+jr*jr))
         phase1=randu(seed)*2.0d0 *pi
         phase2=randu(seed)*2.0d0 *pi
         phase3=randu(seed)*2.0d0 *pi
         phase4=randu(seed)*2.0d0 *pi
         DO j = jsta,jend
            DO i = 1,n
               R1(i,j) = R1(i,j) &
                       + cos(2.0d0*(ir-1)*pi*(dble(i)-1)/dble(n)+phase1) &
                       * cos(2.0d0*(jr-1)*pi*(dble(j)-1)/dble(n)+phase2)/kr
               R2(i,j) = R2(i,j) &
                       + sin(2.0d0*(ir-1)*pi*(dble(i)-1)/dble(n)+phase3) &
                       * cos(2.0d0*(jr-1)*pi*(dble(j)-1)/dble(n)+phase4)    ! No kr because we're forcing vz
            END DO
         END DO
         END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
         CALL fftp2d_real_to_complex(planrc,R2,vz,MPI_COMM_WORLD)
         CALL energy(ps,ener,1)
         CALL energy(vz,enerv,0)
         CALL MPI_BCAST(ener,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         CALL MPI_BCAST(enerv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         tmp=u0/sqrt(ener)
         tmp1=v0/sqrt(enerv)
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  ps(j,i) = tmp*ps(j,i)
                  vz(j,i) = tmp1*vz(j,i)
               ELSE
                  ps(j,i) = 0.0d0
                  vz(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE ! stat
         print*,'READING...',stat
         ini = int((stat-1)*tstep)
         dump = dble(ini)/dble(sstep)+1
         times = 0
         timet = 0
         timec = 0
         timec2 = 0
         timecorr = 0 
        
         jt = 48+int(dump/1000)
         jc = 48+int(dump/100)-int(dump/1000)*10
         jd = 48+int(dump/10)-int(dump/100)*10
         ju = 48+int(dump)-int(dump/10)*10

         ic = 48+int(float(stat)/100)
         id = 48+int(float(stat)/10)-int(float(stat)/100)*10
         iu = 48+int(stat)-int(float(stat)/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)

         OPEN(1,file=trim(idir) // '/ps.' // node // '.' &
                           // c // d // u //'.out',form='unformatted')
         READ(1) R1
         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
         OPEN(1,file=trim(idir) // '/vz.' // node // '.' &
                           // c // d // u //'.out',form='unformatted')
         READ(1) R2
         CALL fftp2d_real_to_complex(planrc,R2,vz,MPI_COMM_WORLD)

!        CALL energy(ps,ener,1)
!        CALL energy(phi,enerphi,0)
!          IF (myrank.eq.0) THEN
!           print*, "DBG:",ener,enerphi
!          ENDIF
!

      ENDIF ! stat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         FORCING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (iflow.eq.1) THEN
         kdn = kup
         DO j = jsta,jend
            DO i = 1,n
               R1(i,j) = sin(2*kup*pi*(dble(i)-1)/dble(n)) &
                       + sin(2*kup*pi*(dble(j)-1)/dble(n))
               R2(i,j) = cos(2*kup*pi*(dble(i)-1)/dble(n)) &  ! Following 3D ABC forcing
                       + sin(2*kup*pi*(dble(j)-1)/dble(n))
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,R1,fp,MPI_COMM_WORLD)
         CALL energy(fp,eneru,1)
         CALL MPI_BCAST(eneru,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         tmp=fp0/sqrt(eneru)
         CALL fftp2d_real_to_complex(planrc,R2,fz,MPI_COMM_WORLD)
         CALL energy(fz,enerv,0)
         CALL MPI_BCAST(enerv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp1=fz0/sqrt(enerv)
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
                  fp(j,i) = tmp*fp(j,i)
                  fz(j,i) = tmp1*fz(j,i)
               ELSE
                  fp(j,i) = 0.0d0
                  fz(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
         ELSEIF (iflow.eq.2) THEN
         seed1 = seed+1
         CALL const_inj(ps,kdn,kup,fp0,fp,1,seed1)
         CALL const_inj(vz,kdn,kup,fz0,fz,0,seed)  
         ELSEIF (iflow.eq.3) THEN
         CALL CFL_condition(CFL,ps,vz,nu,nn,nuv,nnv,omega,dt)
         CALL rand_force(kdn,kup,fp0,dt,seed,1,fp)
         CALL rand_force(kdn,kup,fz0,dt,seed1,0,fz)
         ENDIF ! iflow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*, "DBG pre RK:",ener,enerphi
!          ENDIF
!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################

 RK : DO t = ini,step
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*,"DBG top RK pre CFL",ener,enerphi
!         ENDIF
       CALL CFL_condition(CFL,ps,vz,nu,nn,nuv,nnv,omega,dt)
! Every 'cstep' steps, generates external files 
! to check consistency and convergence. See the 
! cond_check subroutine for details.
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*,"DBG top RK",ener,enerphi
!         ENDIF
          IF (timec.eq.cstep) THEN   
              timec = 0
              CALL cond_check(ps,vz,fp,fz,time,nn,nu,mm,hnu,nnv,nuv,mmv,hnuv,kup,omega)
          ENDIF

          IF (iflow.eq.3) THEN
           CALL rand_force(kdn,kup,fp0,dt,seed,1,fp)
           CALL rand_force(kdn,kup,fz0,dt,seed1,0,fz)
          ENDIF ! iflow3 
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*, "DBG post iflow:",ener,enerphi
!          ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Every 'sstep' steps, generates external files 
! with the power spectrum

         IF (times.eq.sstep) THEN
            times = 0
            ju = ju+1
            IF (ju.eq.58) THEN
               ju = 48
               jd = jd+1
            ENDIF
            IF (jd.eq.58) THEN
               jd = 48
               jc = jc+1
            ENDIF
            IF (jc.eq.58) THEN
               jc = 48
               jt = jt+1
            ENDIF
            th= char(jt)
            c = char(jc)
            d = char(jd)
            u = char(ju)
            ext4 = th // c // d // u 
           CALL spectrum(ps,fp,phi,ext4,odir)
           CALL transfers(ps,ext4,odir)

           IF (myrank.eq.0) THEN
            OPEN(1,file='time_spec.txt',position='append')
            WRITE(1,13) ext4,time
   13       FORMAT( A4,    F12.6)
            CLOSE(1)
           ENDIF
         ENDIF

! Every 'tstep' steps, stores the results of the integration

         IF (timet.eq.tstep) THEN
            timet = 0
            iu = iu+1
            IF (iu.eq.58) THEN
               iu = 48
               id = id+1
            ENDIF
            IF (id.eq.58) THEN
               id = 48
               ic = ic+1
            ENDIF
            c = char(ic)
            d = char(id)
            u = char(iu)
            ext = c // d // u
            tmp=1.0d0/dble(n)**2
            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = ps(j,i)*tmp
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(odir) // '/ps.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = phi(j,i)*tmp
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(odir) // '/phi.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
            DO i = ista,iend
               DO j = 1,n
                  C1(j,i) = ps(j,i)*ka2(j,i)*tmp
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(odir) // '/ww.' // node // '.' &
                 // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
           IF (myrank.eq.0) THEN
           OPEN(1,file='time_field.txt',position='append')
           WRITE(1,12) c//d//u,time
!   10      FORMAT( I10,   A30)
!   11      FORMAT( F12.6, A30)
   12      FORMAT( A3,    F12.6) 
           CLOSE(1)
           ENDIF      

         ENDIF

         timet = timet+1
         times = times+1
         timec = timec+1
         timec2 = timec2+1
!         timecorr = timecorr+dt
         time = time+dt
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*, "DBG post sstep:",ener,enerphi
!          ENDIF


! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

         DO i = ista,iend
            DO j = 1,n
               C1(j,i) = ps(j,i)
               C2(j,i) = phi(j,i)
            END DO
         END DO


! Runge-Kutta step 2

         DO o = ord,1,-1
!         CALL energy(C1,ener,1)
!         CALL energy(C2,enerphi,0)
!           IF (myrank.eq.0) THEN
!            print*,"DBG top ord",ener,enerphi
!           endif
!!!!!!!!!!!!!! iflow2!!!!! Change the forcing to keep constant energy
       IF (iflow.eq.2) THEN
               seed1 = seed+1
               CALL const_inj(C1,kdn,kup,fp0,fp,seed1)
       ENDIF ! iflow2

         CALL laplak2(C1,C4)               ! make - W_2D
         CALL poisson(C1,C4,C5)            ! -curl(u_2D x w_2D)
         CALL pmult(C2,C2,C3)               ! makes phi^2
         CALL pmult(C3,C2,C6)               ! makes phi^3
        ! Normalizing w term so that |w|^2 = 1.
         CALL energy(C1,tmp4,2) ! |w|^2
         CALL MPI_BCAST(tmp4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL pmult(C4/sqrt(tmp4),C2,C3)               ! makes -phi*w

         DO i = ista,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).gt.tiny)) THEN
            tmp1 = dt/dble(o) 
            tmp2 = 1.0d0/ka2(j,i) 

            !  ps
            tmp3 = (1.0d0 +(nu*ka2(j,i)**nn + hnu*tmp2**mm )*tmp1)
            C1(j,i) =  ps(j,i)+((-C5(j,i))*tmp2+fp(j,i))*tmp1 
            C1(j,i) =  C1(j,i)/tmp3

            ! phi
            C2(j,i) = phi(j,i) + (mu*C2(j,i)-cphi*C6(j,i)-fphi*C3(j,i) &
                                        -dphi*ka2(j,i)*C2(j,i))*tmp1
            ELSE  
            C1(j,i) = 0.0d0
            C2(j,i) = 0.0d0
            ENDIF
            END DO
         END DO
!         CALL energy(C1,ener,1)
!         CALL energy(C2,enerphi,0)
!           IF (myrank.eq.0) THEN
!                print*,"DBG end ord",ener,enerphi
!           ENDIF
         END DO  ! ord

! Runge-Kutta step 3
! Copies the result from the auxiliary matrixes into ps, az

         CALL derivk2(C1,ps,0)
         CALL derivk2(C2,phi,0)

!          IF (timec2.eq.cstep) THEN
!              timec2 = 0
!              CALL test_sub(time,ps,nu,nn,dt)
!          ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END DO RK
!##############  END OF MAIN LOOP ###################

!
! End of Runge-Kutta

      CALL MPI_FINALIZE(ierr)
      CALL fftp2d_destroy_plan(plancr)
      CALL fftp2d_destroy_plan(planrc)
      DEALLOCATE( R1 )
      DEALLOCATE( R2 )
      DEALLOCATE( ps )
      DEALLOCATE( vz )
      DEALLOCATE( fp )
      DEALLOCATE( fz )
      DEALLOCATE( C1 )
      DEALLOCATE( C2 )
      DEALLOCATE( C3 )
      DEALLOCATE( C4 )
      DEALLOCATE( C5 )
      DEALLOCATE( C6 )
      
      DEALLOCATE( ka,ka2 )

      END PROGRAM 2D3C
