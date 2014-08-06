!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial rho field for the gravity wave test

!> @detail The kernel computes initial rho perturbation field for the Klemp & Skamarock 
!>         nonhydrostatic gravity wave test

module initial_rho_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only: r_def
use argument_mod,            only: arg_type, &          ! the type
                                   gh_rw, v3, fe, cells ! the enums

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_rho_kernel_type
  private
  type(arg_type) :: meta_args(1) = [ &
       arg_type(gh_rw,v3,fe,.false.,.false.,.false.,.false.) &
       ]
  integer :: iterates_over = cells

contains
  procedure, nopass :: initial_rho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_rho_kernel_type
   module procedure initial_rho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_rho_code
contains

type(initial_rho_kernel_type) function initial_rho_kernel_constructor() result(self)
  return
end function initial_rho_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[inout] rho Real array, the actual data
subroutine initial_rho_code(nlayers,ndf,map,rho)
  
  !Arguments
  integer, intent(in) :: nlayers, ndf
  integer, intent(in) :: map(ndf)
  real(kind=r_def), intent(inout) :: rho(*)

  !Internal variables
  integer               :: df, k
   
  ! compute the pointwise theta profile
  do k = 0, nlayers-1
    do df = 1, ndf
       rho(map(df) + k) = 0.0_r_def
    end do
  end do
  
end subroutine initial_rho_code

end module initial_rho_kernel_mod
