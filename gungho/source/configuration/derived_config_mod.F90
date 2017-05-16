!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief module containing some variables derived from namelist input
!!        that cant yet be computed in namelist code
module derived_config_mod

use constants_mod,     only: i_def

  implicit none

  private
  integer(i_def), public, protected :: si_bundle_size
  integer(i_def), public, protected :: bundle_size
  public :: set_derived_config

contains
  subroutine set_derived_config()
    use formulation_config_mod, only: eliminate_p

    implicit none
    bundle_size = 3_i_def
    if ( eliminate_p ) then
      si_bundle_size = 3_i_def
    else
      si_bundle_size = 4_i_def
    end if

  end subroutine set_derived_config

end module derived_config_mod

