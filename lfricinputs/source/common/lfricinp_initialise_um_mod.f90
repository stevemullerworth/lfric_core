! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This module has an underlying dependency on LFRic's logging module

MODULE lfricinp_initialise_um_mod

! Intrinsic modules
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

! Shumlib modules
USE f_shum_file_mod, ONLY: shum_file_type
USE f_shum_field_mod, ONLY: shum_field_type

! UM2LFRic modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_um_parameters_mod, ONLY: fnamelen, um_integer64

! LFRic modules
USE log_mod,     ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,             &
                       log_scratch_space

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_initialise_um, lfricinp_finalise_um, um_input_file

TYPE(shum_file_type), SAVE :: um_input_file

CONTAINS

! DEPENDS ON: c_shum_byteswap.o
! This is required to force fcm-make to compile the C code; whilst the built-in
! dependency analyser successfully works out that it needs to compile the
! Fortran side of the byte-swapping code, it requires an explicit statement
! to force it to compile the C part of the byte-swapping code. This is
! currently the approved way of linking Fortran and C in fcm-make.

!-------------------------------------------------------------------------------

SUBROUTINE lfricinp_initialise_um(fname)

  IMPLICIT NONE

  CHARACTER(LEN=fnamelen), INTENT(IN) :: fname
  CHARACTER(LEN=*), PARAMETER :: routinename='lfricinp_initialise_um'

  ! Load the UM file
  CALL log_event('Loading file '//TRIM(fname), LOG_LEVEL_INFO)
  CALL shumlib(routinename//'::open_file', um_input_file%open_file(fname),&
                print_on_success=.TRUE._C_BOOL)

  CALL shumlib(routinename//'::open_file', um_input_file%read_header(),  &
                print_on_success=.TRUE._C_BOOL)

END SUBROUTINE lfricinp_initialise_um

SUBROUTINE lfricinp_finalise_um()
  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER :: routinename='lfricinp_finalise_um'

  ! Close file
  CALL shumlib(routinename//'::close_file', um_input_file%close_file() )

END SUBROUTINE lfricinp_finalise_um

END MODULE lfricinp_initialise_um_mod
