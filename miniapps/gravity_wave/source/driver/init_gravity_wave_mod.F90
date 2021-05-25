!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief init functionality for gravity wave simulations

!> @details Handles init of prognostic and coordinate fields

module init_gravity_wave_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use gravity_wave_alg_mod,           only : gravity_wave_alg_init
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO

  implicit none

  contains

  subroutine init_gravity_wave( mesh_id, wind, buoyancy, pressure )
    implicit none

    integer(i_def), intent(in) :: mesh_id

    ! Prognostic fields
    type( field_type ), intent(inout) :: wind
    type( field_type ), intent(inout) :: buoyancy
    type( field_type ), intent(inout) :: pressure

    call log_event( 'gravity_wave: Initialising miniapp ...', LOG_LEVEL_INFO )

    ! Get a reference out of the depository collection and put it into
    ! the collection held in state

    call gravity_wave_alg_init(mesh_id, wind, pressure, buoyancy)

    call log_event( 'gravity_wave: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_gravity_wave

end module init_gravity_wave_mod
