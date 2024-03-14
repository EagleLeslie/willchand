PROGRAM main

  IMPLICIT NONE

  INTEGER :: i,j,k,l
  INTEGER :: ntimesteps, nions, totatoms, xmax, ymax, zmax, dummy
  INTEGER, ALLOCATABLE, DIMENSION(:) :: natom
  REAL*8 :: start, finish, density, scale
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volume, rhoz
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: agridx, agridy, agridz, aintf
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: acell, xyz
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: phi, zintf, mnorm
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: nphi
  CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: atype
  LOGICAL :: sim_type, file_exists

  ! Time how long program runs
  CALL CPU_TIME(start)

  ! Find parameters: number of timesteps, number of ions in system, and total atoms
  CALL find_params(ntimesteps, nions, totatoms)

  ! Allocate arrays according to parameters
  ALLOCATE(natom(nions))
  ALLOCATE(atype(nions))
  ALLOCATE(volume(ntimesteps))
  ALLOCATE(acell(ntimesteps,3,3))
  ALLOCATE(xyz(ntimesteps,totatoms,3))

  ! Read in XDATCAR
  CALL read_xdat(scale,ntimesteps,nions,totatoms,acell,volume,natom,xyz,atype, sim_type)

  xmax = INT(CEILING(acell(1,1,1)))
  ymax = INT(CEILING(acell(1,2,2)))
  zmax = INT(CEILING(acell(1,3,3)))

  ALLOCATE(phi(ntimesteps,xmax,ymax,zmax))
  ALLOCATE(nphi(ntimesteps,xmax,ymax,zmax,3))
  ALLOCATE(mnorm(xmax,ymax,zmax,3))
  ALLOCATE(rhoz(zmax))
  ALLOCATE(zintf(ntimesteps,xmax,ymax,2))

  ALLOCATE(agridx(totatoms,xmax))
  ALLOCATE(agridy(totatoms,ymax))
  ALLOCATE(agridz(totatoms,1))
  ALLOCATE(aintf(ntimesteps,totatoms))

  ! True == NPT or NPH simulation
  ! False == NVT simulation
  PRINT*, sim_type

  ! Shift everything by half box length so that the interfaces aren't on box edges
  ! OPEN(UNIT=800, FILE='ratom',STATUS='UNKNOWN')
  DO i = 1, ntimesteps
      DO j = 1, SUM(natom)
          ! WRITE(800,*) xyz(i,j,1)*acell(1,1,1), xyz(i,j,2)*acell(1,2,2), xyz(i,j,3)*acell(1,3,3)
          xyz(i,j,3) = xyz(i,j,3) + 0.25
          IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
          IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
      END DO
  END DO
  ! CLOSE(800)

  ! Check to see if the Density field has already been calculated
  INQUIRE(FILE='DENFIELD',EXIST=file_exists)
  WRITE(*,*) file_exists
  IF(file_exists) THEN
    WRITE(*,*) "File DENFIELD exists. Reading in densities to phi array."
    OPEN(UNIT=15,FILE='DENFIELD',STATUS='OLD',ERR=501)
    OPEN(UNIT=16,FILE='NPHI',STATUS='OLD',ERR=502)
    DO i = 1,ntimesteps
      DO j = 1,zmax
        DO k = 1,xmax
          DO l = 1,ymax
            READ(15,*,END=500) dummy, dummy, dummy, dummy, density
            READ(16,*,END=503) dummy, dummy, dummy, dummy, nphi(i,k,l,j,:)
            phi(i,k,l,j) = density
            ! WRITE(*,*) phi(i,j,k,l)
          END DO
        END DO
      END DO
      IF (mod(i,200) .eq. 0) write(*,'(a19,a7,i7,a6,i7,a)',advance='NO') &
      & 'Reading in phi(r,t)', ' time',i,' of ', ntimesteps,' '//CHAR(13)
    END DO
    500 CONTINUE
    501 CONTINUE
    502 CONTINUE
    503 CONTINUE
    CLOSE(15)
    CLOSE(16)

    WRITE(*,*) "Now calculating the interfaces."
    ! Calculate the interfaces
    CALL interface(ntimesteps, nions, natom, volume(1), acell(1,3,3), xmax, ymax, zmax, rhoz, phi, zintf)
    ! Calculate atom proximity to interfaces
    CALL proximity(ntimesteps, nions, natom, acell, xyz, xmax, ymax, zmax, nphi, zintf, agridx, agridy, agridz, aintf)

  ELSE
    WRITE(*,*) "File DENFIELD does not exist. Moving on to calculate the density field right now. "
    ! Calculate density field
    CALL calc_denfield(ntimesteps, nions, natom, acell, xyz, xmax, ymax, zmax, phi, mnorm, nphi)
    ! Calculate the interfaces
    CALL interface(ntimesteps, nions, natom, volume(1), acell(1,3,3), xmax, ymax, zmax, rhoz, phi, zintf)
    ! Calculate atom proximity to interfaces
    CALL proximity(ntimesteps, nions, natom, acell, xyz, xmax, ymax, zmax, nphi, zintf, agridx, agridy, agridz, aintf)
  END IF

  OPEN(UNIT=777,FILE='ADATCAR',STATUS='UNKNOWN')

  130  FORMAT(f18.10, 5X, f18.10,5X, f18.10,5X, f18.10,5X, f18.10)
  ! Write to ADATCAR file
  WRITE(777,*) "ADATCAR HEADER"
  WRITE(777,*) scale
  WRITE(777,*) acell(1,1,:)
  WRITE(777,*) acell(1,2,:)
  WRITE(777,*) acell(1,3,:)
  WRITE(777,*) atype(:)
  WRITE(777,*) natom(:)

  ! Calulcate proximity of ith atom to the surface
  ! Equation 5 or equation 6
  DO i = 1,ntimesteps
      WRITE(777,*) "Direct configuration=", i
      DO j = 1,SUM(natom(:))
          xyz(i,j,3) = xyz(i,j,3) - 0.25

          IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
          IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
          WRITE(777,130) xyz(i,j,1),xyz(i,j,2),xyz(i,j,3)-0.25, aintf(i,j)

      END DO
  END DO

  CLOSE(777)

  ! Deallocate arrays
  DEALLOCATE(natom, atype, volume, acell, xyz, phi, nphi, mnorm, rhoz, zintf)
  DEALLOCATE(agridx, agridy, agridz, aintf)

  CALL CPU_TIME(finish)
  PRINT*, "Time = ", finish-start, "seconds ( = ", (finish-start)/60, "minutes.)"

END PROGRAM main