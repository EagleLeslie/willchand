REAL FUNCTION zline(zarr,rhoarr,z)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: z
  REAL, INTENT(IN) :: zarr(z),rhoarr(z)
  REAL :: intc, func,y,dy
  INTEGER :: jlo
  COMMON / values / intc

  !WRITE(*,*) 'poop 2'
  CALL polint(zarr(z),rhoarr(z),2,intc,y,dy)

  func = y - intc
  WRITE(*,*) z,zarr(z),y,dy,intc,func
  !CALL hunt(rhoz,z,intc,jlo)

END FUNCTION zline
