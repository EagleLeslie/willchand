SUBROUTINE brakzero(rho,intc,z1,z2,n,iz1,nz)

  IMPLICIT NONE

  ! Find zeros
  REAL*8, INTENT(IN) :: rho(n), z1, z2, intc
  INTEGER, INTENT(IN) :: n
  REAL*8, DIMENSION(:),INTENT(OUT) :: iz1(n)
  REAL*8, ALLOCATABLE,DIMENSION(:,:) :: zsurf
  INTEGER, POINTER :: count
  INTEGER, TARGET,INTENT (OUT) :: nz
  INTEGER :: i,j,z
  REAL*8 :: dens, x, dz, f1, f2,iz

  ALLOCATE(count)

  !dz = (z2 -z1)/n

  f1 = rho(1) - intc
  nz = 0
  count = 0
  DO z = 1,n
    f2 = rho(z) - intc
    !WRITE(*,*) f1,f2,iz,z-1,z,rho(z)
    IF (f1*f2.lt.0) THEN
      nz = nz + 1
      count = nz
      iz = (-f1)*((z-(z-1))/(f2-f1)) + (z-1)
      iz1(count) = iz
    END IF
    f1 = f2
  END DO

  RETURN

  NULLIFY(count)

END SUBROUTINE brakzero
