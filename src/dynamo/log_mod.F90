!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief A simple logging facility.
!>
!> @deprecated This is a minimal implementation which will be replaced when
!>             the project adopts the ESMF.
!>
module log_mod

  ! This is implemented as a simple Fortran90 module rather than a singleton
  ! class as was originally intended. While it may be possible to imlpement
  ! singletons in Fortran it will fight you every step of the way and it
  ! wasn't worth the effort.

  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

  implicit none

  private
  public log_set_info_stream, log_set_alert_stream, log_set_level, log_event

  !> Named logging level.
  !>
  !> Any integer can be used for the logging level but this name represents
  !> a break between level. Generally you will want to use these names.
  !>
  !> @{
  integer, public, parameter :: LOG_LEVEL_ERROR   = 200
  integer, public, parameter :: LOG_LEVEL_WARNING = 150
  integer, public, parameter :: LOG_LEVEL_INFO    = 100
  integer, public, parameter :: LOG_LEVEL_DEBUG   =  50
  integer, public, parameter :: LOG_LEVEL_TRACE   =   0
  !> @}

  !> Space in which to marshal log messages.
  !>
  !> Although any string can be passed to log_event() this space is provided to
  !> prevent a proliferation of work spaces all over the code. It also means
  !> that should 160 characters be found to be insufficient it need only be
  !> changed in one place.
  !>
  character( 160 ), public :: log_scratch_space

  integer, private, parameter :: EXIT_CODE_ON_ERROR = 1

  integer, private :: log_level  = LOG_LEVEL_INFO
  integer, private :: info_unit  = output_unit
  integer, private :: alert_unit = error_unit

contains

  !> Set where information goes.
  !>
  !> If this routine is never called then information will default to standard
  !> out.
  !>
  !> @param unit The unit to send information to
  !>
  subroutine log_set_info_stream(unit)

    implicit none

    integer, intent( in ) :: unit

    info_unit = unit

  end subroutine log_set_info_stream

  !> Set where alerts go.
  !>
  !> If this routine is never called then alerts will default to standard
  !> error.
  !>
  !> @param unit The unit to send alerts to
  !>
  subroutine log_set_alert_stream(unit)

    implicit none

    integer, intent( in ) :: unit

    alert_unit = unit

  end subroutine log_set_alert_stream

  !> Set the level this logger responds to.
  !>
  !> Events ranked lower than the logging level will be accepted and dropped
  !> on the floor.
  !>
  !> @param level The new logging level to adopt.
  !>
  subroutine log_set_level(level)

    implicit none

    integer, intent( in ) :: level

    log_level = level

  end subroutine log_set_level

  !> Log an event to the terminal.
  !>
  !> The event description is sent to the terminal along with timestamp and
  !> level information. For the most serious events (a severity level equal to
  !> or greater than LOG_LEVEL_ERROR), execution of the code will be aborted.
  !>
  !> @param message A description of the event.
  !> @param level   The severity of the event. Defaults to cInfoLevel.
  !>
  subroutine log_event(message, level)

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

    implicit none

    character (*), intent( in ) :: message
    integer,       intent( in ) :: level

    integer        :: unit
    character (5)  :: tag
    character (8)  :: date_string
    character (10) :: time_string
    character (5)  :: zone_string

    if (level >= log_level) then
      call date_and_time( date=date_string, time=time_string, zone=zone_string)

      select case (level)
        case ( : LOG_LEVEL_DEBUG - 1)
          unit = info_unit
          tag  = 'TRACE'
        case (LOG_LEVEL_DEBUG : LOG_LEVEL_INFO - 1 )
          unit = info_unit
          tag  = 'DEBUG'
        case ( LOG_LEVEL_INFO : LOG_LEVEL_WARNING - 1 )
          unit = info_unit
          tag  = 'INFO'
        case ( LOG_LEVEL_WARNING : LOG_LEVEL_ERROR - 1)
          unit = alert_unit
          tag  = 'WARN'
        case ( LOG_LEVEL_ERROR : )
          unit = alert_unit
          tag  = 'ERROR'
      end select

      write (unit, '(A,A,A,A,A,A,A)') date_string, time_string, zone_string, &
                                      ':', tag, ': ', trim( message )

      ! If the severity level of the event is serious enough, stop the code.
      if ( level >= LOG_LEVEL_ERROR )then
        error stop EXIT_CODE_ON_ERROR
      end if

    end if

  end subroutine log_event

end module log_mod
