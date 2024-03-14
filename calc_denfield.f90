SUBROUTINE calc_denfield(ntimesteps, nions, natom, acell, xyz, xmax, ymax, zmax, phi, mnorm, nphi)

    USE common

    IMPLICIT NONE

    INTEGER :: i,j
    INTEGER :: x,y,z,z2,t0,tf,tdiff
    INTEGER, INTENT(IN) :: ntimesteps, nions, xmax, ymax, zmax
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    REAL*8, DIMENSION(:,:,:), INTENT(IN) :: xyz(ntimesteps,SUM(natom),3), acell(ntimesteps,3,3)
    REAL*8 :: a,b,c,xi,denfi,xref,yref,zref,xgrad,ygrad,zgrad,rdiff,dx,dy,dz,totrho
    REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: phi(ntimesteps,xmax,ymax,zmax), mnorm(xmax,ymax,zmax,3)
    REAL*8, DIMENSION(:,:,:,:,:), INTENT(OUT) :: nphi(ntimesteps,xmax,ymax,zmax,3)

    OPEN(UNIT=100,FILE='DENFIELD',STATUS='UNKNOWN')
    OPEN(UNIT=101,FILE='NPHI',STATUS='UNKNOWN')

    ! Initialize data
    a = acell(1,1,1)
    b = acell(1,2,2)
    c = acell(1,3,3)
    xi = 2.4

    WRITE(*,*) 'xmax: ', xmax, ';  ymax: ', ymax, ';  zmax: ', zmax
    WRITE(*,*) 'Coarse graining length = ', xi

    t0 = 1
    tf = ntimesteps
    tdiff = tf - t0 + 1

    ! Calculate density field
    ! Equation 2 and equation 3
    DO i = t0, tf
        DO z = 1, zmax
            DO x = 1,xmax
                DO y = 1,ymax
                    totrho = 0
                    xgrad = 0
                    ygrad = 0
                    zgrad = 0
                    DO j = 1,natom(1)
                        ! convert grid points from angstroms to fractal units
                        xref = REAL(x)/REAL((xmax))
                        yref = REAL(y)/REAL((ymax))
                        zref = REAL(z)/REAL((zmax))
                        !WRITE(*,*) zref

                        CALL den_field(a,b,c,xyz(i,j,1),xyz(i,j,2),xyz(i,j,3),xref,yref,zref,xmax,ymax,zmax,xi,denfi)
                        totrho = totrho + denfi
                        phi(i,x,y,z) = totrho

                        dx = xref - xyz(i,j,1) - ANINT(xref - xyz(i,j,1))
                        dy = yref - xyz(i,j,2) - ANINT(yref - xyz(i,j,2))
                        dz = zref - xyz(i,j,3) - ANINT(zref - xyz(i,j,3))

                        rdiff = SQRT((a*(xmax)*dx**2. + b*(ymax)*dy**2.+ c*(zmax)*dz**2.))

                        ! Calculate normal wrt grid
                        xgrad = xgrad + (-dx*a* EXP(-(rdiff**2.)/(2.*(xi**2.)) ) / ( (8**(1./2.))*(pi**(3./2.)) &
                        & * (xi**2.)**(5./2.) ) )

                        ygrad = ygrad + (-dy*b* EXP(-(rdiff**2.)/(2.*(xi**2.)) ) / ( (8**(1./2.))*(pi**(3./2.)) &
                        & * (xi**2.)**(5./2.) ) )

                        zgrad = zgrad + (-dz*c* EXP(-(rdiff**2.)/(2.*(xi**2.)) ) / ( (8**(1./2.))*(pi**(3./2.)) &
                        & * (xi**2.)**(5./2.) ) )

                    END DO
                nphi(i,x,y,z,1) = xgrad/SQRT(xgrad**2. + ygrad**2. + zgrad**2.)
                nphi(i,x,y,z,2) = ygrad/SQRT(xgrad**2. + ygrad**2. + zgrad**2.)
                nphi(i,x,y,z,3) = zgrad/SQRT(xgrad**2. + ygrad**2. + zgrad**2.)

                z2 = z - (0.25*c)
                IF (z2.gt.1) z2 = z2 - c
                IF (z2.lt.0) z2 = z2 + c

                WRITE(100,*) i,x,y,z2, phi(i,x,y,z) ! Averaged density field
                WRITE(101,*) i,x,y,z, nphi(i,x,y,z,:) ! Write out normalizaed gradient rho
                IF (mod(i,100) .eq. 0) write(*,'(a19,a9,i7,a)',advance='NO') &
                & 'Computing rho(r,t)', ' of time',i,' '//CHAR(13)
                mnorm(x,y,z,1) = (SUM(nphi(:,x,y,z,1),DIM=1)/tdiff)
                mnorm(x,y,z,2) = (SUM(nphi(:,x,y,z,2),DIM=1)/tdiff)
                mnorm(x,y,z,3) = (SUM(nphi(:,x,y,z,3),DIM=1)/tdiff)
                END DO
            END DO
        END DO
    END DO

    ! Write out mean normalized gradient rho
    ! DO x = 1,xmax
    !     DO y = 1,ymax
    !         DO z = 1,zmax
    !             WRITE(23,*) x,y,z,mnorm(x,y,z,:)
    !         END DO
    !     END DO
    ! END DO

    CLOSE(100)
    CLOSE(101)

END SUBROUTINE calc_denfield