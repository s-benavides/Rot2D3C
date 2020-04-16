!*****************************************************************
      SUBROUTINE initial(a1,a2,a3,v1,v2,v3)
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

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a1,a2,a2
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: v1,v2,v3
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j





      END SUBROUTINE initial

