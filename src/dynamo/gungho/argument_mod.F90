!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief The argument type to hold kernel metadata required by the psy layer.

!> @details Metadata for the kernels. For each field passed to a kernel 
!> the psy layer needs to know how this field is to be accessed.
!> read, write etc, to which function space it belongs and what
!> stencil it operates over. These are the three integers of
!> the arg_type and the values are then one of the parameters
!> defined in this module. 
!> field metadata also has three logicals controlling whether the psy layer
!> needs to pass the basis function, the differential basis function,
!> and the guassian quadrature type.
!> Another metadatum which describes the kernel, not the fields
!> is what the kernel will iterate over. Usually cells, sometimes
!> all the dofs.

module argument_mod
  implicit none

! access descriptors
  integer, public, parameter :: gh_read  = 1 
  integer, public, parameter :: gh_write = 2
  integer, public, parameter :: gh_rw    = 3
  integer, public, parameter :: gh_inc   = 4
  integer, public, parameter :: gh_sum   = 5
  integer, public, parameter :: gh_min   = 6
  integer, public, parameter :: gh_max   = 7

! vspace labels
  integer, public, parameter :: v0 = 1
  integer, public, parameter :: v1 = 2
  integer, public, parameter :: v2 = 3
  integer, public, parameter :: v3 = 4
  integer, public, parameter :: any_space = 0

! stencil label
  integer, public, parameter :: fe = 1 

! kernel iterator
  integer, public, parameter :: cells     = 1
  integer, public, parameter :: all_dofs  = 2

  type, public :: arg_type
     integer :: arg_intent
     integer :: vspace
     integer :: stencil
     logical :: basis
     logical :: diff_basis
     logical :: nodal_coords
     logical :: gaussian_quad
  end type arg_type

end module argument_mod

