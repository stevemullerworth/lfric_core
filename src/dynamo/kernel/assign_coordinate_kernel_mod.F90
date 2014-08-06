!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the coordinate field from the grid node values

module assign_coordinate_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, any_space, fe, cells ! the enums                                 
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: assign_coordinate_kernel_type
  private
  type(arg_type) :: meta_args(3) = [  &
       arg_type(gh_write,any_space,fe,.false.,.false.,.false.,.false.),        &
       arg_type(gh_write,any_space,fe,.false.,.false.,.false.,.false.),        &
       arg_type(gh_write,any_space,fe,.false.,.false.,.false.,.false.)         &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass :: assign_coordinate_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface assign_coordinate_kernel_type
   module procedure assign_coordinate_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public assign_coordinate_code
contains

type(assign_coordinate_kernel_type) function assign_coordinate_kernel_constructor() result(self)
  return
end function assign_coordinate_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[out] chi_1 Real array the first component of the genralised physical coordinates 
!! @param[out] chi_2 Real array the second component of the genralised physical coordinates 
!! @param[out] chi_3 Real array the third component of the genralised physical coordinates 
!! @param[in] vertex_coords Real array. the cordinates of the vertices
!! @param[in] x_node Real array. the cordinates of the nodal points in the function space
subroutine assign_coordinate_code(nlayers,ndf,nverts,map,chi_1,chi_2,chi_3, &
                                  vertex_coords,chi_hat_node,chi_hat_vert)

  !Arguments
  integer, intent(in) :: nlayers, ndf, nverts
  integer, intent(in) :: map(ndf)  
  real(kind=r_def), intent(out) :: chi_1(*), chi_2(*), chi_3(*)
  real(kind=r_def), intent(in)  :: vertex_coords(3,nverts,nlayers)
  real(kind=r_def), intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)

  !Internal variables
  integer          :: k, df, dfk, vert
  
  real(kind=r_def) :: interp_weight, x, y, z
  
  ! compute the representation of the coordinate field
  do k = 0, nlayers-1
    do df = 1, ndf 
! compute interpolation weights
      x = 0.0_r_def
      y = 0.0_r_def
      z = 0.0_r_def
      do vert = 1,nverts
        interp_weight = (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                       *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                       *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
      
        x = x + interp_weight*vertex_coords(1,vert,k+1)
        y = y + interp_weight*vertex_coords(2,vert,k+1)
        z = z + interp_weight*vertex_coords(3,vert,k+1)
      end do
! For spherical domains we need to project x,y,z back onto spherical shells
! This still needs to be done
    
      dfk = map(df)+k 
      chi_1(dfk) = x
      chi_2(dfk) = y
      chi_3(dfk) = z
    end do
  end do
  
end subroutine assign_coordinate_code

end module assign_coordinate_kernel_mod
