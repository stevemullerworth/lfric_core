!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_set_target_data_point_mod

    use constants_mod, only : i_def, r_def, imdi
    use field_mod, only : field_type, field_proxy_type

    implicit none
    public

contains



    !------------------------------------------------------------------------------
    !> invoke_set_target_data_point_kernel: explicitly set a single dof within a target field

    subroutine invoke_set_target_data_point_kernel(field, target_cell, target_layer, value)

        use local_mesh_mod, only : local_mesh_type
        use mesh_mod, only : mesh_type

        implicit none

        type(field_type), intent(inout) :: field
        real(kind = r_def), intent(in) :: value
        integer(kind = i_def), intent(in) :: target_cell, target_layer
        type(mesh_type), pointer :: mesh => null()
        type(local_mesh_type), pointer :: local_mesh => null()
        integer, pointer :: map_w3(:, :) => null()

        integer(i_def) :: local_id
        type(field_proxy_type) :: field_proxy

        field_proxy = field%get_proxy()

        !        field_proxy%data(target_variable)=value

        mesh => field_proxy%vspace%get_mesh()
        local_mesh => mesh%get_local_mesh()
        map_w3 => field_proxy%vspace%get_whole_dofmap()

        ! find out what the local mesh id is for the global ID provided (this could be working on a partition of the
        ! global mesh.
        local_id = local_mesh%get_lid_from_gid(target_cell)
        ! only update it if it is in our grid
        if(local_id /= IMDI) then

            field_proxy%data(map_w3(1, local_id) + target_layer - 1) = value
            call field_proxy%set_dirty()
        end if

    end subroutine invoke_set_target_data_point_kernel

end module psykal_lite_set_target_data_point_mod
