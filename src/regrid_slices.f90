!=================================================================
      PROGRAM rd_slices
!=================================================================
!=================================================================
      IMPLICIT NONE

!

!
!     fields 

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: Vold ! a slice of the field
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: Vnew ! the total NxN field


! Some auxiliary parameters

      REAL :: rfile,rnode

      INTEGER :: Nnodes,Mnodes,inode
      INTEGER :: n,ns
      INTEGER :: m,ms
      INTEGER :: is
      INTEGER :: iname
      INTEGER :: ifile
      INTEGER :: i ,j ,k
      INTEGER :: i2,j2,k2
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju
      INTEGER :: kc,kd,ku

      CHARACTER     :: c,d,u
      CHARACTER*3   :: node,ext
      CHARACTER*100 :: ldir
        
      !ldir  = "thin2D_vx."                ! directory and file head (here we load the x component of the velocity field)
      ifile = 1                       ! File number
      Nnodes= 64                       ! Number of nodes used (ie number of slices )
      Mnodes= 128                       ! Number of nodes to be used (ie number of slices )
      n     = 2048                       ! Old Grid size
      m     = 4096                       ! New grid size   m>n
      ns    = n/Nnodes                  ! Old Slice thickness  
      ms    = m/Mnodes                  ! New Slice thickness ms=<ns

 
      ALLOCATE( Vold(n,ns) )
      ALLOCATE( Vnew(m,ms) )

        print*,"OLD FILES N=",n," Nodes=",Nnodes," slice=",ns
        print*,"NEW FILES N=",m," Nodes=",Mnodes," slice=",ms
        print*,"ONE OLD SLICE GIVES ",Mnodes/Nnodes, " SLICES" 
        print*,"ONE POINT GIVES     ",m/n," POINTS"           


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       DO iname=1,3 ! LOOP OVER HEADERS
       IF (iname.eq.1) ldir  = "thin2D_ps."
       IF (iname.eq.2) ldir  = "thin2D_vx."
       IF (iname.eq.3) ldir  = "thin2D_vy."

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

         print*," " !NOTE: I am changing the '.out' to '.out2' in case input file # = output file #.
         print*,"READING SLICE ",trim(ldir)//node//'.'//ext//'.out2'!HINT:~$ rename '.out' '.out2' *.out 
         OPEN(1,file=trim(ldir) // node // '.' &
                                // ext // '.out2',form='unformatted')
         READ(1) Vold
         CLOSE(1)
        
!         print*, " DBG Vold ", SHAPE(Vold)
         DO is = 1,Mnodes/Nnodes
!         print*,"DBG making slice #",is,is+inode*ns/ms,ns/ms,ns,ms
         DO k=1,ms
            k2=(k+(is-1)*ms+1)*n/m
!            print*,"k=",k," k2=",k2  
            DO j=1,m
               j2=(j+1)*n/m
!               print*,"j=",j," j2=",j2
               Vnew(j,k)=Vold(j2,k2) 
            END DO
         END DO
!         print*, " DGB Vnew", SHAPE(Vnew) 
          rnode = float(inode*Mnodes/Nnodes+is-1)
          ic = 48+int(rnode/100)
          id = 48+int(rnode/10 )-int(rnode/100)*10
          iu = 48+int(rnode)-int(rnode/10)*10
          c = char(ic)
          d = char(id)
          u = char(iu)
          node = c // d // u
         print*,"writing SLICE ",trim(ldir) //node//'.001.out'
         OPEN(1,file=trim(ldir) // node //'.001.out',form='unformatted')
         WRITE(1) Vnew
         CLOSE(1)
         END DO !is

       END DO !inode 
       END DO !iname


      DEALLOCATE( Vold)
      DEALLOCATE( Vnew)

      END PROGRAM rd_slices
