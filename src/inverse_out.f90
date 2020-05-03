!=================================================================
      PROGRAM inverse_out
!=================================================================
!=================================================================
      IMPLICIT NONE

!     Program to invert fields, as well as x-axis. Used to test if 
!     certain solutions' mirror versions are also solutions.

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: Vold ! a slice of the field
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: Vnew ! the total NxN field


! Some auxiliary parameters

      REAL :: rfile,rnode

      INTEGER :: Nnodes,inode
      INTEGER :: n,ns
      INTEGER :: iname
      INTEGER :: ifile
      INTEGER :: j ,k
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju

      CHARACTER     :: c,d,u
      CHARACTER*3   :: node,ext
      CHARACTER*100 :: ldir
        
      ifile = 93                       ! File number
      Nnodes = 16                      ! Number of nodes used (ie number of slices )
      n     = 512                     ! Grid size
      ns    = n/Nnodes                 ! Slice thickness  

 
      ALLOCATE( Vold(n,ns) )
      ALLOCATE( Vnew(n,ns) )

        print*,"FILES N=",n," Nodes=",Nnodes," slice=",ns


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       DO iname=1,3 ! LOOP OVER HEADERS
       IF (iname.eq.1) ldir  = "./ps."
       IF (iname.eq.2) ldir  = "./vz."
       IF (iname.eq.3) ldir  = "./ww."

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       DO inode=0,Nnodes-1  ! LOOP OVER ALL OLD SLICES %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        Make the file name of slice # inode
         rfile = float(ifile)
         rnode = float(inode)

         jc = 48+int(rfile/100)
         jd = 48+int(rfile/10)-int(rfile/100)*10
         ju = 48+int(rfile)-int(rfile/10)*10

         ic = 48+int(rnode/100)
         id = 48+int(rnode/10 )-int(rnode/100)*10
         iu = 48+int(rnode)-int(rnode/10)*10

         c = char(ic)
         d = char(id)
         u = char(iu)
         node = c // d // u

         c = char(jc)
         d = char(jd)
         u = char(ju)
         ext = c // d // u

         print*," " 
         print*,"READING SLICE ",trim(ldir)//node//'.'//ext//'.out' 
         OPEN(1,file=trim(ldir) // node // '.'// ext // '.out', &
                                                form='unformatted')
         READ(1) Vold
         CLOSE(1)

        ! Take the negative value, and also flip the x-axis
        
         DO k=1,ns
            DO j=1,n
               Vnew(j,k)=-1.0d0*Vold(n-j,k) 
            END DO
         END DO

         print*,"writing SLICE ",trim(ldir) //node//'.001.out'
         OPEN(1,file=trim(ldir) // node //'.001.out',form='unformatted')
         WRITE(1) Vnew
         CLOSE(1)

       END DO !inode 
       END DO !iname


      DEALLOCATE( Vold)
      DEALLOCATE( Vnew)

      END PROGRAM inverse_out
