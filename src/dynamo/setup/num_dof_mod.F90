!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Computes the global and local number of dofs for the 4 element spaces
!>        V0..V3
!>
module num_dof_mod

  use mesh_generator_mod, only : nface_g,nedge_g,nvert_g

contains 

  !> Compute the local and global number of dofs.
  !>
  !> @param ncells Number of cells.
  !> @param nlayers Number of vertical cells.
  !> @param k Order of RT space ( = 0 for lowest order )
  !>
  subroutine num_dof_init( ncells, nlayers, k, v_unique_dofs, v_dof_entity )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    integer, intent( in ) :: ncells
    integer, intent( in ) :: nlayers
    integer, intent( in ) :: k

    integer, intent( out ) :: v_unique_dofs(4,2) ! there are 4 vspaces
    integer, intent( out ) :: v_dof_entity(4,0:3)

    ! numbers of dofs in each space for an element 
    integer :: nv0,nv0_cell,nv0_face,nv0_edge,nv0_vert,            &
               nv1,nv1_cell,nv1_face,nv1_edge,                     &
               nv2,nv2_cell,nv2_face,                              &
               nv3,nv3_cell

    ! global numbers of unique dofs
    integer :: nv0_g, nv1_g, nv2_g, nv3_g

    integer :: ndof_entity_v0(0:3), ndof_entity_v1(0:3),           &
               ndof_entity_v2(0:3), ndof_entity_v3(0:3)

    ! local values
    nv0 = (k+2)*(k+2)*(k+2)
    nv0_cell = k*k*k
    nv0_face = k*k
    nv0_edge = k
    nv0_vert = 1

    nv1 = 3*(k+2)*(k+2)*(k+1)
    nv1_cell = 3*k*k*(k+1)
    nv1_face = 2  *k*(k+1)
    nv1_edge =       (k+1)

    nv2  = 3*(k+2)*(k+1)*(k+1)
    nv2_cell = 3*k*(k+1)*(k+1)
    nv2_face =     (k+1)*(k+1)

    nv3 = (k+1)*(k+1)*(k+1)
    nv3_cell = nv3

    ! global numbers of dofs per function space
    nv3_g = ncells*nlayers*nv3_cell
    nv2_g = ncells*nlayers*nv2_cell + nface_g*nv2_face
    nv1_g = ncells*nlayers*nv1_cell + nface_g*nv1_face + nedge_g*nv1_edge
    nv0_g = ncells*nlayers*nv0_cell + nface_g*nv0_face + nedge_g*nv0_edge + nvert_g*nv0_vert

    ! populate the returned arrays
    v_unique_dofs(1,1) = nv0_g
    v_unique_dofs(2,1) = nv1_g
    v_unique_dofs(3,1) = nv2_g
    v_unique_dofs(4,1) = nv3_g

    v_unique_dofs(1,2) = nv0
    v_unique_dofs(2,2) = nv1
    v_unique_dofs(3,2) = nv2
    v_unique_dofs(4,2) = nv3

    ! Number of dofs per mesh entity for each space
    ndof_entity_v0(:) = (/ nv0_vert, nv0_edge, nv0_face, nv0_cell /)
    ndof_entity_v1(:) = (/ 0       , nv1_edge, nv1_face, nv1_cell /)
    ndof_entity_v2(:) = (/ 0       , 0       , nv2_face, nv2_cell /)
    ndof_entity_v3(:) = (/ 0       , 0       , 0       , nv3_cell /)

    !populate the returned arrays
    v_dof_entity(1,:) = ndof_entity_v0(:)
    v_dof_entity(2,:) = ndof_entity_v1(:)
    v_dof_entity(3,:) = ndof_entity_v2(:)
    v_dof_entity(4,:) = ndof_entity_v3(:)

    ! diagnostic output
    write( log_scratch_space, '(A, I0, A, I0)' ) &
        'ncells = ', ncells, ', nlayers = ', nlayers
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( '   space     |   V0   |   V1   |   V2   |   V3   |', &
                    LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'global dof    ', nv0_g, '   ', nv1_g, '   ', nv2_g, '   ', nv3_g
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'local dof     ', nv0, '   ', nv1, '   ', nv2, '   ', nv3
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof in volume ', nv0_cell, '   ', nv1_cell, '   ', nv2_cell, &
        '   ', nv3
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on face   ', nv0_face, '   ', nv1_face, '   ', nv2_face, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on edge   ', nv0_edge, '   ', nv1_edge, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on vert   ', nv0_vert, '   ', 0, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine num_dof_init

end module num_dof_mod
