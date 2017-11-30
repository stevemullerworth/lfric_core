!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page gung_ho GungHo Program
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate gungho functionality.

!> @details Calls various init subroutines to create a mesh, function spaces
!> and then prognostic fields on those function spaces.
!> Controls the main timestepping loop, initialises and calls
!> timestepping algorithms and also handles output of diagnostics at
!> specified intervals as well as checkpoint restart files.

program gungho

  use cli_mod,           only : get_initial_filename
  use gungho_driver_mod, only : initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program gungho
