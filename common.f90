MODULE common
  IMPLICIT NONE

  ! common variables - Fortran 90 syntax
  REAL*8, PARAMETER :: pi=3.1415926535897932385
  REAL*8, PARAMETER :: mfe = 9.2732796E-26, mmg = 4.0359398E-26, mo = 2.6566962E-26 ! masses in kg
  REAL*8, PARAMETER :: amufe = 55.845, amumg = 24.305, amuo = 15.999, amumgo = 40.304 ! masses in amu
  REAL*8, PARAMETER :: avog = 6.02214076E23 !mol^1
  SAVE

  INTERFACE
    SUBROUTINE den_field(a,b,c,xx,yy,zz,xref,yref,zref,xmax,ymax,zmax,xi,phi)
      IMPLICIT NONE
      ! a,b,c :: lattice constants of system
      ! xx,yy,zz :: position of atom x,y,z components
      ! xmax,ymax,zmax :: limits of cell
      ! phi :: returned density field
      REAL*8, INTENT(IN) :: xi,xref,yref,zref,xx,yy,zz,a,b,c
      INTEGER, INTENT(IN) :: xmax,ymax,zmax
      REAL*8, INTENT(OUT) :: phi
      REAL*8 :: dx,dy,dz,rdiff
      !REAL, PARAMETER :: pi=3.1415926535897932385

    END SUBROUTINE den_field
  END INTERFACE

END MODULE common

SUBROUTINE den_field(a,b,c,xx,yy,zz,xref,yref,zref,xmax,ymax,zmax,xi,phi)
  IMPLICIT NONE
  ! a,b,c :: lattice constants of system
  ! xx,yy,zz :: position of atom x,y,z components
  ! xmax,ymax,zmax :: limits of cell
  ! phi :: returned density field
  REAL*8, INTENT(IN) :: xi,xref,yref,zref,xx,yy,zz,a,b,c
  INTEGER, INTENT(IN) :: xmax,ymax,zmax
  REAL*8, INTENT(OUT) :: phi
  REAL*8 :: dx,dy,dz,rdiff
  REAL*8, PARAMETER :: pi=3.1415926535897932385

  dx = (xref - xx) - ANINT(xref - xx)
  dy = (yref - yy) - ANINT(yref - yy)
  dz = (zref - zz) - ANINT(zref - zz)

  ! Willard and Chandler density calculation
  rdiff = SQRT((a*(xmax)*dx**2. + b*(ymax)*dy**2.+ c*(zmax)*dz**2.))

  phi = ((2*pi*(xi**2.))**(-3./2.)) * EXP(-(rdiff**2.)/(2*(xi**2.)))

  !WRITE(*,*) dx,dy,dz,rdiff,phi

END SUBROUTINE den_field
