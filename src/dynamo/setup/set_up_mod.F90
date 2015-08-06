!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Sets up data required for the the function spaces

!> @details This code generates a mesh and determines the basis functions and
!> dofmaps. This will be replaced with code that reads this in from a mesh
!> generation and paritioning pre-processor stage.  

! There are no tests for this code as this will be replaced.

module set_up_mod

  use constants_mod,              only : i_def, r_def, str_def, PI, QUAD
  use function_space_mod,         only : function_space_type
  use reference_element_mod,      only : reference_cube, &
                                         reference_element, nfaces, nedges, nverts
  use num_dof_mod,                only : num_dof_init
  use basis_function_mod,         only : get_basis, &
              w0_nodal_coords, w1_nodal_coords, w2_nodal_coords, w3_nodal_coords

  use dofmap_mod,                 only : get_dofmap, get_orientation, &
              w0_dofmap, w1_dofmap, w2_dofmap, w3_dofmap
  implicit none
  
contains 

!> @brief Generates a mesh and determines the basis functions and dofmaps
!> @details This will be replaced with code that reads the information in
!> @param[out] mesh Mesh object to run model on
  subroutine set_up(mesh)

    use log_mod,         only : log_event, LOG_LEVEL_INFO
    use slush_mod,       only : element_order,                                &
                                l_spherical, w_unique_dofs, w_dof_entity,     &
                                dx, dy, num_cells_x, num_cells_y,             &
                                xproc, yproc, local_rank, total_ranks,        &
                                l_fplane, f_lat

    use mesh_mod,        only : mesh_type
    use partition_mod,   only : partition_type,                 &
                                partitioner_interface,          &
                                partitioner_cubedsphere_serial, &
                                partitioner_cubedsphere,        &
                                partitioner_biperiodic

    use global_mesh_mod, only : global_mesh_type

    implicit none

    character(len = str_def) :: filename
    type (global_mesh_type)  :: global_mesh
    type (partition_type)    :: partition

    type (mesh_type), intent(out) :: mesh

    procedure (partitioner_interface), pointer :: partitioner_ptr => null ()


    real(r_def)    :: dz
    integer(i_def) :: nlayers
    !Get the processor decomposition
    !Code is not set up to run in parallel - so hardcode for now
    xproc = 1
    yproc = 1
!> @todo Eventually xproc and yproc will be inputted into Dynamo (and not hard-coded).
!>       When this happens their values will need to be checked to make sure they are
!>       sensible  - e.g. that they are consistent with the values of num_cells_x
!>       and num_cells_y 

    ! hard-coded these numbers are
    l_fplane      = .true.
    num_cells_x   = 12
    num_cells_y   = 12
    nlayers       = 5
    element_order = 0
    l_spherical   = .true.
! Horizontal spacings for cartesian grid    
    dx = 6000.0_r_def 
    dy = 2000.0_r_def
! Vertical spacing for all grids    
    dz = 2000.0_r_def

    filename = 'ugrid_quads_2d.nc' 
    call log_event( "set_up: generating/reading the mesh", LOG_LEVEL_INFO )

    reference_element = QUAD

    ! Setup reference cube  
    call reference_cube()

    ! Generate the global mesh and choose a partitioning strategy by setting
    ! a function pointer to point at the appropriate partitioning routine
    if ( l_spherical ) then
      global_mesh=global_mesh_type( filename )
      partitioner_ptr => partitioner_cubedsphere_serial
    else
      global_mesh=global_mesh_type( num_cells_x, &
                                    num_cells_y, &
                                    dx,& 
                                    dy )
      partitioner_ptr => partitioner_biperiodic
      if ( l_fplane ) f_lat = PI/4.0_r_def
    end if

    ! Generate the partition object
    partition=partition_type( global_mesh, &
                              partitioner_ptr, &
                              xproc, &
                              yproc, &
                              1, &
                              local_rank, &
                              total_ranks)

    mesh = mesh_type(partition, global_mesh, nlayers, dz)

! -----------------------------------------------------------
! Initialise FE elements on the mesh constructed above
! really another pre-processor step
! ----------------------------------------------------------

    ! initialise numbers of dofs    
    call num_dof_init( mesh, element_order, w_unique_dofs, w_dof_entity )
         
    call log_event( "set_up: computing basis functions", LOG_LEVEL_INFO )

    ! read the values of the basis functions. 
    call get_basis( k=element_order,             &
                    w_unique_dofs=w_unique_dofs, &
                    w_dof_entity=w_dof_entity )  

    call log_event( "set_up: computing the dof_map", LOG_LEVEL_INFO )

    ! compute the dof maps for each function space
    call get_dofmap( mesh=mesh,                  &
                     w_dof_entity=w_dof_entity,  &
                     w_unique_dofs=w_unique_dofs)
    
    ! compute cell local orientations for vector spaces
    call get_orientation( mesh, w_unique_dofs, w_dof_entity )

    return

  end subroutine set_up

end module set_up_mod
