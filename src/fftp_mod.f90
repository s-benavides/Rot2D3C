!=================================================================
! MODULES for FFTP
! Parallel Fast Fourier Transform in 2D and 3D
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
!                             7 when using GNU compilers
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb               
      INCLUDE 'fftw_f77.i'
      INTEGER, PARAMETER  :: ikind = 8
      INTEGER, PARAMETER  :: csize = 16
      TYPE FFTPLAN
         INTEGER(kind=ikind) :: planr,planc
         INTEGER :: n
         INTEGER, DIMENSION (:), POINTER :: itype1, itype2
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE mpivars
      INCLUDE 'mpif.h'
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
