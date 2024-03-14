SUBROUTINE proximity(ntimesteps,nions,natom,acell,xyz,xmax,ymax,zmax,nphi,zintf,agridx,agridy,agridz,aintf)

    IMPLICIT NONE

    INTEGER :: i,j,k,l,m
    INTEGER :: x,y,z,z2,t0,tf,tdiff
    INTEGER, INTENT(IN) :: ntimesteps, nions, xmax, ymax, zmax
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    REAL*8 :: a,b,c,xref,yref,zref,zref2,dx,dy,dz,dz2,dz3,mdx,mdy,mdz,delta
    REAL*8, DIMENSION(:,:,:), INTENT(IN) :: xyz(ntimesteps,SUM(natom),3), acell(ntimesteps,3,3)
    REAL*8, DIMENSION(:,:,:,:), INTENT(IN) :: zintf(ntimesteps,xmax,ymax,2)
    REAL*8, DIMENSION(:,:,:,:,:), INTENT(IN) :: nphi(ntimesteps,xmax,ymax,zmax,3)
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: agridx(SUM(natom),xmax), agridy(SUM(natom),xmax), agridz(SUM(natom),1)
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: aintf(ntimesteps,SUM(natom))

    ! allocate arrays for 
    ! agridx, agridy, agridz (:,:)
    ! aintf(:,:,:)

    ! Initialize data
    a = acell(1,1,1)
    b = acell(1,2,2)
    c = acell(1,3,3)
    delta = 1.5

    t0 = 1
    tf = ntimesteps
    tdiff = tf - t0 + 1

    DO i = t0,tf
        DO j = 1,SUM(natom(:))
            DO x = 1,xmax
                ! convert grid points from angstroms to fractal units
                xref = REAL(x)/REAL((xmax))
                dx = ABS(xref - xyz(i,j,1) - ANINT(xref - xyz(i,j,1)))
                agridx(j,x) = dx
                DO y = 1,ymax
                    ! convert grid points from angstroms to fractal units
                    yref = REAL(y)/REAL((ymax))
                    dy = ABS(yref - xyz(i,j,2) - ANINT(yref - xyz(i,j,2)))
                    agridy(j,y) = dy

                    ! convert grid points from angstroms to fractal units
                    zref = (zintf(i,x,y,1))/REAL((zmax))
                    zref2 = (zintf(i,x,y,2))/REAL((zmax))

                    dz = ABS(zref - xyz(i,j,3) - ANINT(zref - xyz(i,j,3)))
                    dz2 = ABS(zref2 - xyz(i,j,3)- ANINT(zref2 - xyz(i,j,3)))

                    IF (dz.lt.dz2) agridz(j,1) = zintf(i,x,y,1) !- (0.25*c)
                    IF (dz2.lt.dz) agridz(j,1) = zintf(i,x,y,2) !- (0.25*c)
                    IF (agridz(j,1).eq.0) agridz(j,1) = NINT(c)

                END DO
            END DO

            k = MINLOC(agridx(j,:),DIM=1)
            l = MINLOC(agridy(j,:),DIM=1)
            m = agridz(j,1)

            dx = ((REAL(k)/xmax) - xyz(i,j,1)) - ANINT((REAL(k)/xmax) - xyz(i,j,1))
            dy = ((REAL(l)/ymax) - xyz(i,j,2)) - ANINT((REAL(l)/ymax) - xyz(i,j,2))
            dz3 = ((REAL(m)/zmax) - xyz(i,j,3)) - ANINT((REAL(m)/zmax) - xyz(i,j,3))

            mdx = dx/SQRT(dx**2. + dy**2. + dz3**2.)
            mdy = dy/SQRT(dx**2. + dy**2. + dz3**2.)
            mdz = dz3/SQRT(dx**2. + dy**2. + dz3**2.)

            aintf(i,j) = (a*dx*nphi(i,k,l,m,1)) + (b*dy*nphi(i,k,l,m,2)) + (c*dz3*nphi(i,k,l,m,3))
        END DO
        WRITE(*,'(i7,a4,i7,a)',advance='no') i, 'of',tf, ' '//CHAR(13)
    END DO


END SUBROUTINE proximity