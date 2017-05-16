!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

program log_mod_error_test

  use ESMF,    only : ESMF_Initialize, ESMF_Finalize
  use log_mod, only : log_event, LOG_LEVEL_ERROR

  call ESMF_Initialize()
  call log_event( 'An error was logged.', LOG_LEVEL_ERROR )
  call ESMF_Finalize()

end program log_mod_error_test
