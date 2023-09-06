PROGRAM willchand
  USE common

  IMPLICIT NONE

  INTEGER :: i,j,k,l,m
  INTEGER :: totatoms, nions
  INTEGER :: x, y, z,z2,nz,tdiff
  INTEGER, PARAMETER :: mint=3
  INTEGER, ALLOCATABLE :: ntimesteps, xmax, ymax, zmax, nbins, tau,nsurf,t0,tend
  REAL*8 :: xi, r, totrho, volume, rdiff,rdiff2,dx,dy,dz,dz2,dz3,xgrad,ygrad,zgrad,zgrad2,sumdelt,rmax,n,delta
  REAL*8 :: a, b, c, d0, d1, d2, alat, xref, yref, zref,zref2,intc,denfi,iz,fz,test,boo,fecount1,fecount2
  REAL*8 :: counter,mnx,mny,mnz,mgcount1,ocount1,mgap,oap,feap,az1,az2,mid,mdx,mdy,mdz,readin
  REAL*8:: x1,x2,zz,dyy,sz,dzz,comx, comy, comz, mass_sum, drift,molfe,molmg,molo,totmol
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: xyz, aintf
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: phi, zsurf,zintf,mintf,mnorm
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: nphi
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ibins,acount,molpct,wtpct,tester,tester2,tester4
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhoz,com,dr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atype,natom
  INTEGER, ALLOCATABLE,DIMENSION(:,:) :: tester3
  LOGICAL :: file_exists
  !CHARACTER(len=20) :: typat1, typat2, typat3
  CHARACTER(len=5), ALLOCATABLE, DIMENSION(:) :: aname
  !REAL, PARAMETER :: pi=3.1415926535897932385

  OPEN(UNIT=10,FILE='XDATCAR',STATUS='OLD',ERR=200)
  OPEN(UNIT=100,FILE='DENFIELD',STATUS='UNKNOWN')
  OPEN(UNIT=800, FILE='ratom')
  OPEN(UNIT=777,FILE='ADATCAR',STATUS='UNKNOWN')
  OPEN(UNIT=1111,FILE='intf1',STATUS='UNKNOWN')
  OPEN(UNIT=1112,FILE='intf2',STATUS='UNKNOWN')
  OPEN(UNIT=2222,FILE='atomcar',STATUS='UNKNOWN')
  OPEN(UNIT=2223,FILE='acount',STATUS='UNKNOWN')
  OPEN(UNIT=2224,FILE='wtper',STATUS='UNKNOWN')
  OPEN(UNIT=2225,FILE='molpct',STATUS='UNKNOWN')

  ! Get number of timesteps from XDATCAR
  CALL EXECUTE_COMMAND_LINE("grep 'Direct' XDATCAR | tail -1| awk '{print $3}' > ntime")
  OPEN(1,FILE='ntime')
  READ(1,*,END=100) i
  100 CONTINUE
  ntimesteps = i
  CLOSE(1,STATUS='DELETE')

  WRITE(*,*) 'Total timesteps: ', ntimesteps

  CALL EXECUTE_COMMAND_LINE("sed '7q;d' XDATCAR > blah")
  CALL EXECUTE_COMMAND_LINE("awk '{print NF}' blah | sort -nu | tail -n 1 > nions")
  OPEN(1,FILE='nions')
  READ(1,*,END=101) i
  101 CONTINUE
  nions = i
  CALL EXECUTE_COMMAND_LINE("rm blah nions")
  WRITE(*,*) 'Total atom types: ', nions

  ALLOCATE(aname(nions))
  ALLOCATE(natom(nions))

  ! Get cell and natom from XDATCAR
  READ(10,*)
  READ(10,*) alat
  READ(10,*) a, d1, d2
  READ(10,*) d0, b, d2
  READ(10,*) d0, d1, c
  READ(10,*) aname(:)
  READ(10,*) natom(:)
  READ(10,*)

  totatoms = SUM(natom(:))
  WRITE(*,*) 'Total atoms: ', totatoms
  WRITE(*,*) a, d1, d2
  WRITE(*,*) d0, b, d2
  WRITE(*,*) d0, d1, c
  WRITE(*,*) aname(:)
  WRITE(*,*) natom(:)

  volume = a * b * c
  WRITE(*,*) 'Volume: ', volume

  ! xi = coarse-graining length; must choose a value
  !WRITE(*,*) 'Choose a xi value: '
  !READ(*,*) xi
  xi = 2.4

  WRITE(*,*) 'z-axis length:', c
  rmax = c
  nbins = 25
  delta = rmax/nbins

  ALLOCATE(xyz(ntimesteps,totatoms,3))
  ALLOCATE(aintf(ntimesteps,totatoms,3))
  ALLOCATE(com(ntimesteps))
  ALLOCATE(atype(nions))
  ALLOCATE(dr(nbins))
  ALLOCATE(ibins(ntimesteps,nbins))
  ALLOCATE(acount(ntimesteps,3))

  DO i = 1,nions
    IF (i.eq.1) atype(i) = natom(i)
    IF (i.eq.2) atype(i) = natom(i) + natom(i-1)
    IF (i.eq.3) atype(i) = natom(i) + natom(i-1) + natom(i-2)
    IF (i.eq.4) atype(i) = natom(i) + natom(i-1) + natom(i-2) + natom(i-3)
  END DO

  mass_sum = (natom(1)*mfe) !+ (natom2*mmg) + (natom3*mo)

  ! Read in atom positions
  ! Shift everything by half box length and calculate center of mass
  DO i = 1, ntimesteps
!    comz = 0.0
    DO j = 1, totatoms
      READ(10,*) xyz(i,j,:)
      WRITE(800,*) xyz(i,j,1)*a, xyz(i,j,2)*b, xyz(i,j,3)*c
      xyz(i,j,3) = xyz(i,j,3) + 0.25
      IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
      IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
!      IF (j.le.atype(1)) THEN
!        comz = comz + (mfe * xyz(i,j,3))
!      END IF
    END DO
!    com(i) = comz/mass_sum
    !WRITE(*,*) i, com(i)
    READ(10,*,END=300)
  END DO
  300 CONTINUE

  ! Initialize data
  xmax = INT(CEILING(a))
  ymax = INT(CEILING(b))
  zmax = INT(CEILING(c))
  x1 = 1
  x2 = zmax

  WRITE(*,*) 'xmax: ', xmax, ';  ymax: ', ymax, ';  zmax: ', zmax

  ALLOCATE(phi(ntimesteps,xmax,ymax,zmax))
  ALLOCATE(nphi(ntimesteps,xmax,ymax,zmax,3))
  ALLOCATE(zsurf(ntimesteps,xmax,ymax,zmax))
  ALLOCATE(zintf(ntimesteps,xmax,ymax,2))
  ALLOCATE(mintf(ntimesteps,xmax,ymax,2))
  ALLOCATE(mnorm(xmax,ymax,zmax,3))
  ALLOCATE(rhoz(zmax))

  !ALLOCATE(denfi(zmax))

  ! interface constant
  WRITE(*,*) natom(:)
  INQUIRE(FILE='INTFC',EXIST=file_exists)
  WRITE(*,*) file_exists
  IF(file_exists) THEN
    CALL EXECUTE_COMMAND_LINE("awk 'NR==1 {print $1; exit}' INTFC > tt")
    OPEN(3,FILE='tt')
    READ(3,*,END=333) readin
    333 CONTINUE
    intc = readin
    WRITE(*,*) "File INTFC read into program."
    WRITE(*,*) 'Interface constant: ', intc
    CLOSE(3,STATUS='DELETE')
  ELSE
    WRITE(*,*) "File INTFC does not exist. Using default interface constant value: &
     Mean of metal and oxide phase divided by volume "

     mid = (natom(1) + natom(2))/2
     WRITE(*,*) 'Mid point', mid
     intc = mid/volume
     WRITE(*,*) 'Interface constant: ', intc
     WRITE(*,*) 'Interface constant is about 1/2 bulk density (Willard and Chandler, 2010).'
  END IF

  !ntimesteps = 100

  ! Write out atom positions to a file for jupyter notebook interface figure
  ! subtract center of mass drift
!  DO i = 1,ntimesteps
!    WRITE(*,*) 'Time:', i, 'initial center of mass:', com(1),&
!    & 'current center of mass:', com(i), 'drift:', (com(i) - com(1))
!    WRITE(*,*) 'Time:', i, 'initial center of mass:', (com(1))*c,&
!    & 'shifted current center of mass:', (com(i))*c, 'drift:', (com(i) - com(1))*c
!    DO j = 1,totatoms
!      !WRITE(*,*) i, xyz(i,j,3), com(1), com(i), com(i) - com(1), com(1) - com(i)
!      IF (com(i).gt.com(1)) THEN
!        xyz(i,j,3) = xyz(i,j,3) - (com(i)-com(1))
!      ELSE IF (com(i).lt.com(1)) THEN
!        xyz(i,j,3) = xyz(i,j,3) + (com(i) - com(1))
!      END IF
!      WRITE(777,*) xyz(i,j,:)
!      WRITE(800,*) xyz(i,j,1)*a, xyz(i,j,2)*b, xyz(i,j,3)*c
!    END DO
!  END DO

!  DO i = 1, ntimesteps
!    comz = 0.0
!    DO j = 1, natom1
!      comz = comz + (mfe * xyz(i,j,3))
!    END DO
!    com(i) = comz/mass_sum
!    WRITE(*,*) i,com(i)
!  END DO

  WRITE(*,*)

  !ntimesteps = 20000

  t0 = 1 !ntimesteps - 50
  tend = ntimesteps
  tdiff = tend - t0 + 1

  ALLOCATE(molpct(ntimesteps,3))
  ALLOCATE(wtpct(ntimesteps,3))

  ! Calculate density field
  ! Equation 2 and equation 3
  DO i = t0,tend
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
            !rdiff = SQRT((dx**2. + dy**2.+ dz**2.))

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
          !WRITE(*,*) dx,dy,dz,xgrad,ygrad,zgrad
          !WRITE(*,*) i,x,y,z,nphi(i,x,y,z,1),nphi(i,x,y,z,2),nphi(i,x,y,z,3)

          z2 = z - (0.25*c)
          IF (z2.gt.1) z2 = z2 - c
          IF (z2.lt.0) z2 = z2 + c

          WRITE(100,*) i,x,y,z2, phi(i,x,y,z), SUM(phi(:,x,y,z),DIM=1)/tdiff!ntimesteps
          WRITE(*,'(a19,a3,i3,a3,i3,a3,i3,a9,i7,a)',advance='NO') 'Computing rho(r,t)', &
          & 'x:',x,'y:',y,'z:',z,' of time',i,' '//CHAR(13)
          mnorm(x,y,z,1) = (SUM(nphi(:,x,y,z,1),DIM=1)/tdiff)!ntimesteps)
          mnorm(x,y,z,2) = (SUM(nphi(:,x,y,z,2),DIM=1)/tdiff)!ntimesteps)
          mnorm(x,y,z,3) = (SUM(nphi(:,x,y,z,3),DIM=1)/tdiff)!ntimesteps)
        !  WRITE(*,*) mnorm(x,y,z,:)
        END DO
      END DO
    END DO
  END DO

  OPEN(UNIT=99,FILE='adens')
  DO i = tend,tend
    DO z = 1,zmax
      DO x = 1,xmax
        DO y = 1,ymax
          WRITE(99,*) i,x,y,z,phi(i,x,y,z), SUM(phi(:,x,y,z),DIM=1)/tdiff
        END DO
      END DO
    END DO
  END DO
  CLOSE(99)

  WRITE(*,*)
  WRITE(*,*) 'poop'

  ! Find interfaces
  ! Equation 4
  DO i = t0,tend
    DO x = 1,xmax
      DO y = 1,ymax
      !  DO z = 1,zmax
          CALL brakzero(phi(i,x,y,:),intc,0.0,REAL(zmax),zmax,rhoz,nz)
          !WRITE(*,*) nz
          DO z = 1,nz
              zsurf(i,x,y,z) = rhoz(z)
              zintf(i,x,y,1) = rhoz(1)
              zintf(i,x,y,2) = rhoz(2)
              IF (nz.gt.3) zintf(i,x,y,2) = 0
          END DO
          WRITE(1111,*) i,x,y,zintf(i,x,y,1) - (0.25*c)
          WRITE(1112,*) i,x,y,zintf(i,x,y,2) - (0.25*c)
          !WRITE(*,*) x,y,zintf(i,x,y,1) - (0.25*c), zintf(i,x,y,2) - (0.25*c),nz
          !zintf(i,x,y,1) = zintf(i,x,y,1) - (0.25*c)
          !zintf(i,x,y,2) = zintf(i,x,y,2) - (0.25*c)
          mintf(i,x,y,1) = SUM(zintf(:,x,y,1),DIM=1)/tdiff!ntimesteps
          mintf(i,x,y,2) = SUM(zintf(:,x,y,2),DIM=1)/tdiff!ntimesteps
      !  END DO
      END DO
    END DO
  END DO

!  OPEN(UNIT=1113,FILE='mintf',STATUS='UNKNOWN')
!  DO x = 1,xmax
!    DO y = 1,ymax
!      WRITE(1113,*) x,y,mintf(tend,x,y,1), mintf(tend,x,y,2)
!    END DO
!  END DO
!  CLOSE(1113)

  WRITE(*,*) 'poop 2'

  ! Calculate normal dotted with distance of atom wrt interface
  ! Part 1 of equation 5 or equation 6

  ALLOCATE(tester(totatoms,xmax))
  ALLOCATE(tester2(totatoms,ymax))
  ALLOCATE(tester3(totatoms,1))
  ALLOCATE(tester4(totatoms,1))

  DO i = t0,tend
    DO j = 1,totatoms
      !xyz(i,j,3) = xyz(i,j,3) - 0.25
      !IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
      !IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
      !WRITE(*,*) xyz(i,j,3)*c
      DO x = 1,xmax
        xref = REAL(x)/REAL((xmax))
        dx = ABS(xref - xyz(i,j,1) - ANINT(xref - xyz(i,j,1)))
        tester(j,x) = dx
        !WRITE(*,*) j,x,xref,xyz(i,j,1),dx,MINVAL(tester(j,:),DIM=1),MINLOC(tester(j,:),DIM=1)
      !  DO z = 1,zmax
          DO y = 1,ymax

            yref = REAL(y)/REAL((ymax))
            dy = ABS(yref - xyz(i,j,2) - ANINT(yref - xyz(i,j,2)))
            tester2(j,y) = dy
            ! toggle between zintf(i,x,y,1) and minf(tend,x,y,1) for instantaneous vs mean interface
            zref = (zintf(i,x,y,1))/REAL((zmax))
            zref2 = (zintf(i,x,y,2))/REAL((zmax))

            dz = ABS(zref - xyz(i,j,3) - ANINT(zref - xyz(i,j,3)))
            dz2 = ABS(zref2 - xyz(i,j,3)- ANINT(zref2 - xyz(i,j,3)))

            !tester4(j,1) = dz
            !tester3(j,1) = zintf(i,x,y,2) - (0.25*c)
            IF (dz.lt.dz2) tester4(j,1) = dz
            IF (dz2.lt.dz) tester4(j,1) = dz2

            IF (dz.lt.dz2) tester3(j,1) = zintf(i,x,y,1) !- (0.25*c)
            IF (dz2.lt.dz) tester3(j,1) = zintf(i,x,y,2) !- (0.25*c)
            IF (tester3(j,1).eq.0) tester3(j,1) = NINT(c)

            !WRITE(*,*) dz*c, dz2*c, zintf(i,x,y,1),zintf(i,x,y,2), tester3(j,1), tester4(j,1)

          END DO
      !  END DO
      END DO

      !dx = MINVAL(tester(j,:),DIM=1)
      !dy = MINVAL(tester2(j,:),DIM=1)
      !dz3 = tester4(j,2)
      k = MINLOC(tester(j,:),DIM=1)
      l = MINLOC(tester2(j,:),DIM=1)
      m = tester3(j,1)


      dx = ((REAL(k)/xmax) - xyz(i,j,1)) - ANINT((REAL(k)/xmax) - xyz(i,j,1))
      dy = ((REAL(l)/ymax) - xyz(i,j,2)) - ANINT((REAL(l)/ymax) - xyz(i,j,2))
      dz3 = ((REAL(m)/zmax) - xyz(i,j,3)) - ANINT((REAL(m)/zmax) - xyz(i,j,3))

      mdx = dx/SQRT(dx**2. + dy**2. + dz3**2.)
      mdy = dy/SQRT(dx**2. + dy**2. + dz3**2.)
      mdz = dz3/SQRT(dx**2. + dy**2. + dz3**2.)

      !WRITE(*,*) dz, dz2, tester4(j,2),dz3

      aintf(i,j,1) = (a*dx*nphi(i,k,l,m,1)) + (b*dy*nphi(i,k,l,m,2)) + (c*dz3*nphi(i,k,l,m,3))
      !aintf(i,j,1) = (a*dx*mnorm(k,l,m,1)) + (b*dy*mnorm(k,l,m,2)) + (c*dz*mnorm(k,l,m,3))
      !WRITE(*,*) j, k,l,m,dx,dy,dz3
      !WRITE(*,*) j,aintf(i,j,1),m, mdz,(xyz(i,j,3))*c,nphi(i,k,l,m,:)
    END DO
    WRITE(*,'(i7,a4,i7,a)',advance='no') i, 'of',tend, ' '//CHAR(13)
  END DO

  WRITE(*,*)
  WRITE(*,*) 'poop 3'
  WRITE(*,*)

  130  FORMAT(f18.10, 5X, f18.10,5X, f18.10,5X, f18.10,5X, f18.10)
  ! Write to ADATCAR file
  WRITE(777,*)
  WRITE(777,*) alat
  WRITE(777,*) a, d1, d2
  WRITE(777,*) d0, b, d2
  WRITE(777,*) d0, d1, c
  WRITE(777,*) aname(:)
  WRITE(777,*) natom(:)

  ! Calulcate proximity of ith atom to the surface
  ! Equation 5 or equation 6
  DO i = t0,tend
    fecount1 = 0
    mgcount1 = 0
    ocount1 = 0
    WRITE(777,*) "Direct configuration=", i
    DO j = 1,totatoms
      WRITE(2222,*) i, j, aintf(i,j,1)
      xyz(i,j,3) = xyz(i,j,3) - 0.25
      IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
      IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
      WRITE(777,130) xyz(i,j,1),xyz(i,j,2),xyz(i,j,3)-0.25, aintf(i,j,1)
      IF (ABS(aintf(i,j,1)).ge.1.2) THEN
        IF (j.lt.atype(1).and.aintf(i,j,1).lt.0) THEN
          fecount1 = fecount1 + 1
        END IF

        IF (j.gt.atype(1).and.j.le.atype(2).and.aintf(i,j,1).lt.0) THEN
          mgcount1 = mgcount1 + 1
        END IF

        IF (j.gt.atype(2).and.aintf(i,j,1).lt.0) THEN
          ocount1 = ocount1 + 1
        END IF
      END IF
    END DO
    acount(i,1) = fecount1
    acount(i,2) = mgcount1
    acount(i,3) = ocount1
    WRITE(2223,*) i, acount(i,1), acount(i,2), acount(i,3), acount(i,1)/volume, acount(i,2)/volume, acount(i,3)/volume
    WRITE(*,'(a15,i7,a5,a5,i7,a)',advance='NO') 'Computing a(', i, ')', &
    & ' of ',tend,' '//CHAR(13)
    !WRITE(*,*) 'Fe count: ', fecount1,'Fe count density: ', fecount1/volume, acount(i,1)
    !WRITE(*,*) 'Mg count: ', mgcount1,'Mg count density: ', mgcount1/volume, acount(i,2)
    !WRITE(*,*) 'O count: ', ocount1,'O count density: ', ocount1/volume, acount(i,3)
    mgap = (mgcount1)*100/(mgcount1+fecount1+ocount1)
    oap = (ocount1)*100/(mgcount1+fecount1+ocount1)
    feap = (fecount1)*100/(mgcount1+fecount1+ocount1)

    molfe = acount(i,1)
    molmg = acount(i,2)
    molo = acount(i,3)
    totmol = molfe + (0.5*(molmg + molo))

    molpct(i,1) = (molfe/totmol)*100
    molpct(i,2) = (0.5 * (molmg + molo))/totmol * 100
    ! molpct(i,3) = (molo/totmol)*100
    wtpct(i,1) = (feap*amufe)/( (0.5*(mgap*amumg) + (oap*amuo)) + (feap*amufe)) * 100
    wtpct(i,2) = ( (0.5*(mgap*amumg) + (oap*amuo)) )/( (0.5*(mgap*amumg) + (oap*amuo)) + (feap*amufe)) * 100
    !wtpct(i,3) = (oap*amuo)/((mgap*amumg) + (feap*amufe) + (oap*amuo)) * 100

    WRITE(2224,*) i, (feap*amufe)/( (0.5*(mgap*amumg) + (oap*amuo)) + (feap*amufe)) * 100, &
    & ( (0.5*(mgap*amumg) + (oap*amuo)) )/( (0.5*(mgap*amumg) + (oap*amuo)) + (feap*amufe)) * 100
    ! & (mgap*amumg)/((mgap*amumg) + (feap*amufe) + (oap*amuo)) * 100, &
    ! & (oap*amuo)/((mgap*amumg) + (feap*amufe) + (oap*amuo)) * 100
    WRITE(2225,*) i, (molfe/totmol)*100, (0.5 * (molmg + molo))/totmol * 100 !, (molmg/totmol)*100, (molo/totmol)*100
  END DO

  WRITE(*,*) SUM(wtpct(:,1),DIM=1)/tdiff,SUM(wtpct(:,2),DIM=1)/tdiff!,SUM(wtpct(:,3),DIM=1)/tdiff
  WRITE(*,*) SUM(molpct(:,1),DIM=1)/tdiff,SUM(molpct(:,2),DIM=1)/tdiff!,SUM(molpct(:,3),DIM=1)/tdiff

  WRITE(*,*)


  ! Calculate mean density profile wrt instantaneous or mean interface
  ! Equation 7
  DO i = t0, tend
    n = 0.0
    DO k = 1,nbins
      r = (k*delta) - (delta/2)
      dr(k) = r
      counter = 0
      n = 0
      DO j = 1,natom(1)
      !  WRITE(*,*) aintf(i,j,2)
        IF (aintf(i,j,1).gt.dr(k-1).and.aintf(i,j,1).lt.dr(k)) THEN
          counter = counter + 1
          n = n + (aintf(i,j,1) - r)
        END IF
      !  IF(abs(aintf(i,j,2)).gt.dr(k-1).and.abs(aintf(i,j,2)).lt.dr(k)) THEN
      !  counter = counter + 1
      !  END IF
      !ibins(i,k) = counter/(a**2.)/(natom1/volume)
      END DO
      ibins(i,k) = counter/(a**2.)/(natom(1)/volume)
      !WRITE(*,*) i, r, counter, counter/(a**2.)/(natom1/volume), counter/(a**2.)
      WRITE(2,*) i, r, SUM(ibins(:,k),DIM=1)/tdiff
    END DO
  END DO

  WRITE(*,*)
  WRITE(*,*) 'poop 4'

  DEALLOCATE(aname,natom,xyz,aintf,atype,dr,acount,zsurf,zintf,mintf,mnorm,phi)
  DEALLOCATE(nphi,rhoz,com,ibins,molpct,wtpct)
  DEALLOCATE(tester,tester2,tester3,tester4)
  !DEALLOCATE(denfi)

  200 CONTINUE
  CLOSE(10)
  CLOSE(100)
  CLOSE(200)
  CLOSE(300)
  CLOSE(800)
  CLOSE(2222)
  CLOSE(2223)
  CLOSE(2224)

END PROGRAM
