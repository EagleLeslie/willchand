SUBROUTINE interface(ntimesteps, nions, natom, volume, c, xmax, ymax, zmax, rhoz, phi, zintf)

    IMPLICIT NONE

    INTEGER :: i,j,k,l,x,y,z,nz
    INTEGER, INTENT(IN) :: ntimesteps, nions, xmax, ymax, zmax
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    REAL*8, INTENT(IN) :: volume,c
    REAL*8, DIMENSION(:), INTENT(IN) :: rhoz(zmax)
    REAL*8, DIMENSION(:,:,:,:), INTENT(IN) :: phi(ntimesteps,xmax,ymax,zmax)
    REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: zintf(ntimesteps,xmax,ymax,2)
    REAL*8 :: intc,mid,readin
    LOGICAL :: file_exists

    OPEN(UNIT=1111,FILE='intf1',STATUS='UNKNOWN')
    OPEN(UNIT=1112,FILE='intf2',STATUS='UNKNOWN')
    
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

    ! Find interfaces
    ! Equation 4
    DO i = 1,ntimesteps
        DO x = 1,xmax
            DO y = 1,ymax
                CALL brakzero(phi(i,x,y,:),intc,0.0,REAL(zmax),zmax,rhoz,nz)
                ! WRITE(*,*) phi(i,x,y,:)
                DO z = 1,nz
                    zintf(i,x,y,1) = rhoz(1)
                    zintf(i,x,y,2) = rhoz(2)
                    IF (nz.gt.3) zintf(i,x,y,2) = 0
                END DO
            WRITE(1111,*) i,x,y,zintf(i,x,y,1) - (0.25*c)
            WRITE(1112,*) i,x,y,zintf(i,x,y,2) - (0.25*c)
            END DO
        END DO
    END DO

    CLOSE(1111)
    CLOSE(1112)

END SUBROUTINE interface