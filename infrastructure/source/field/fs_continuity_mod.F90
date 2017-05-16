!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief Define enumerator variables that describe the different types of continuity.
!>
!> @details Enumerator variables that describe the different types of continuity
!>          that can be used to construct function spaces

module fs_continuity_mod

use constants_mod, only: i_def

implicit none

private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
integer(i_def), public, parameter :: W0      = 100
integer(i_def), public, parameter :: W1      = 101
integer(i_def), public, parameter :: W2      = 102
integer(i_def), public, parameter :: W3      = 103
integer(i_def), public, parameter :: Wtheta  = 104
integer(i_def), public, parameter :: W2V     = 105
integer(i_def), public, parameter :: W2H     = 106
integer(i_def), public, parameter :: Wchi    = 107

character(len=7), public, parameter :: fs_name(100:107) &
          = [character(len=7) :: 'W0', 'W1', 'W2', 'W3', 'Wtheta', 'W2V', 'W2H', 'Wchi']

end module fs_continuity_mod
