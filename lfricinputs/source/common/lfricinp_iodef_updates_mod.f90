! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_iodef_updates_mod

! LFRic Inputs modules
USE lfricinp_unit_handler_mod, ONLY: get_free_unit

! MPI module
USE mpi,  ONLY: mpi_comm_rank, mpi_barrier

! LFRic modules
USE constants_mod,   ONLY: i_def
USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,         &
                           log_scratch_space
IMPLICIT NONE
PRIVATE
PUBLIC :: lfricinp_iodef_set_calendar

CHARACTER(LEN=9), PARAMETER :: iodef_file_name = 'iodef.xml'

CONTAINS

SUBROUTINE lfricinp_iodef_set_calendar(comm, calendar, start_date, time_origin)

IMPLICIT NONE

INTEGER(KIND=i_def), INTENT(IN) :: comm
CHARACTER(LEN=*),    INTENT(IN) :: calendar, start_date, time_origin

INTEGER(KIND=i_def)     :: local_rank, ierr
INTEGER, PARAMETER      :: max_num_lines = 1000
INTEGER                 :: num_lines, i, idx
INTEGER                 :: unit_number
INTEGER                 :: stat
CHARACTER(LEN=512)      :: message = 'IOMSG unset'
CHARACTER(LEN=512)      :: iodef_calendar_definition
CHARACTER(LEN=512)      :: line(max_num_lines), dummy_line

CALL mpi_comm_rank(comm, local_rank, ierr)

IF (local_rank == 0) THEN

  iodef_calendar_definition =  '<calendar type = "'                              &
                               // TRIM(calendar) //                              &
                               '" start_date = "'                                &
                               // TRIM(start_date) //                            &
                               '" time_origin = "'                               &
                               // TRIM(time_origin) //                           &
                               '"/>' 

  CALL get_free_unit(unit_number)
  OPEN(UNIT=unit_number, FILE=iodef_file_name, IOSTAT=stat, IOMSG=message)
  IF (stat /= 0) CALL log_event(message, LOG_LEVEL_ERROR)
  num_lines = 0
  DO 
    READ(unit_number,'(A)', IOSTAT=stat) dummy_line
    IF (stat /= 0) EXIT
    IF (INDEX(dummy_line, 'calendar type') /= 0) THEN
      idx = INDEX(dummy_line, '<')
      dummy_line(idx:) = ADJUSTL(iodef_calendar_definition) 
    END IF
    num_lines = num_lines + 1
    line(num_lines) = dummy_line
  END DO
  CLOSE(UNIT=unit_number)

  CALL log_event('Set calendar definition in iodef file', LOG_LEVEL_INFO)
  CALL get_free_unit(unit_number)
  OPEN(UNIT=unit_number, FILE=iodef_file_name, IOSTAT=stat, IOMSG=message)                 
  IF (stat /= 0) CALL log_event(message, LOG_LEVEL_ERROR)
  DO i = 1, num_lines
    WRITE(unit_number,'(A)') TRIM(line(i))
    FLUSH(unit_number)
  END DO
  CLOSE(UNIT=unit_number)

END IF

CALL mpi_barrier(comm, ierr)

END SUBROUTINE lfricinp_iodef_set_calendar

END MODULE lfricinp_iodef_updates_mod
