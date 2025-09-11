!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief Module container for query functions related to science.
module sci_query_mod

  use constants_mod,   only: i_def
  use global_mesh_mod, only: global_mesh_type
  use log_mod,         only: log_event, log_scratch_space, &
                             log_level_error

  implicit none

  private
  public  :: valid_for_global_model


contains

!> @brief  Queries whether a global mesh is valid for a global model domain.
!>         A negative result does not imply that the mesh is valid for a
!>         regional model.
!> @param[in]  global_mesh
!> @return     answwer
function valid_for_global_model( global_mesh ) result ( answer )

   implicit none

   type(global_mesh_type), intent(in) :: global_mesh
   logical :: answer

   if ( global_mesh%is_topology_periodic() .and. &
        global_mesh%is_coord_sys_ll()      .and. &
        global_mesh%is_geometry_spherical() ) then

      ! Note: if these conditions where satisfied from the
      !       planar mesh generator then the model would be
      !       a torus, though as this is not supported, this
      !       is only .true. for a sphere.
      answer = .true.
  else
      answer = .false.
  end if

end function valid_for_global_model

end module sci_query_mod
