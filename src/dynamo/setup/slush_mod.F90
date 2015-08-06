!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines variables which are either temporary or have not yet been
!> given a home.

!> @details When developing dynamo code, developers should this as a scratch
!> place for variables. This module is intending to be regularly reviewed to
!> move out or delete temporary variables

module slush_mod

  use constants_mod,  only : r_def

  implicit none

  !> For a biperiodic mesh, the number of horizontal cells in the x-direction of
  !> the global domain
  !> For a cubed-sphere mesh it is the number of cells across a face.
  integer :: num_cells_x
  !> Number of horizontal cells in the y-direction of the global domain
  !> (not used for cubed-sphere meshes)
  integer :: num_cells_y
  !> Order of the function space
  integer :: element_order
  !> Flag for whether mesh is on a sphere or not
  logical :: l_spherical
  !> Flag for whether a plane is with constant f (omega)
  logical :: l_fplane

  !> Number of unique dofs in a particular function space (4,:),
  !> either globally (:,1) or per cell (:,2)
  integer :: w_unique_dofs(4,2)
  !> Number of dofs in a particular function space (4,:) per entity (:,0:3)
  integer :: w_dof_entity(4,0:3)

  !> Grid spacing in the x-direction
  real(kind=r_def)  :: dx
  !> Grid spacing in the y-direction
  real(kind=r_def)  :: dy

  real(kind=r_def)  :: f_lat            ! Latitude for f-plane tests

  !> Total number of MPI ranks available
  integer :: total_ranks
  !> The local rank number
  integer :: local_rank

  !> Number of processors along x-direction
  integer :: xproc
  !> Number of processors along y-direction
  integer :: yproc
  
  !> Number of cells that are wholly owned by the partition
  !> (i.e. all dofs in these cells are wholly owned by the partition)
  integer :: num_core
  !> Number of cells that are owned by the partition, but may have dofs that
  !> are also owned by halo cells
  integer :: num_owned
  !> Number of cells in the halo for this partition
  integer :: num_halo

end module slush_mod

