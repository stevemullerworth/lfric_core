!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial u field for the gravity wave test

!> @detail The kernel computes initial u perturbation field for the Klemp & Skamarock 
!>         nonhydrostatic gravity wave test

module initial_u_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use argument_mod,            only : arg_type, &          ! the type
                                    gh_rw, v2, fe, cells ! the enums

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_u_kernel_type
  private
  type(arg_type) :: meta_args(1) = [ &
       arg_type(gh_rw,v2,fe,.false.,.false.,.false.,.false.) &
       ]
  integer :: iterates_over = cells

contains
  procedure, nopass :: initial_u_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_u_kernel_type
   module procedure initial_u_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_u_code
contains

type(initial_u_kernel_type) function initial_u_kernel_constructor() result(self)
  return
end function initial_u_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[inout] u Real array, the actual data
subroutine initial_u_code(nlayers,ndf,map,u)
  
  !Arguments
  integer, intent(in) :: nlayers, ndf
  integer, intent(in) :: map(ndf)
  real(kind=r_def), intent(inout) :: u(*)

  !Internal variables
  integer               :: df, k
   
  ! compute the pointwise theta profile
  do k = 0, nlayers-1
    do df = 1, ndf
       u(map(df) + k) = 0.0_r_def
    end do
  end do
  
end subroutine initial_u_code

end module initial_u_kernel_mod
