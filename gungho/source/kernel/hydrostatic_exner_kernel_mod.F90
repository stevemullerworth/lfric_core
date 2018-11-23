!-----------------------------------------------------------------------------
! Copyright (c) 2018,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes Exner pressure distribution via hydrostatic balance equation (lowest order only)

module hydrostatic_exner_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_READ, GH_WRITE,    &
                                       ANY_SPACE_9, GH_BASIS,          &
                                       GH_DIFF_BASIS,                  &
                                       CELLS, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use planet_config_mod,          only : gravity, cp, rd, p_zero
use idealised_config_mod,       only : test
use kernel_mod,                 only : kernel_type
use fs_continuity_mod,          only : WTHETA, W3
use formulation_config_mod,     only : init_exner_bt
use coord_transform_mod,        only : xyz2llr

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hydrostatic_exner_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,  GH_WRITE,  W3),                             &
       arg_type(GH_FIELD,  GH_READ,   WTHETA),                         &
       arg_type(GH_FIELD,  GH_READ,   WTHETA),                         &
       arg_type(GH_FIELD,  GH_READ,   W3),                             &
       arg_type(GH_FIELD*3,GH_READ, ANY_SPACE_9)                       &
       /)
  type(func_type) :: meta_funcs(1) = (/                             &
       func_type(ANY_SPACE_9, GH_BASIS)                             &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass ::hydrostatic_exner_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface hydrostatic_exner_kernel_type
   module procedure hydrostatic_exner_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public hydrostatic_exner_code
contains

type(hydrostatic_exner_kernel_type) &
   function hydrostatic_exner_kernel_constructor() result(self)
  return
end function hydrostatic_exner_kernel_constructor

!> @brief Computes hydrostatic Exner function
!! @param[in] nlayers Number of layers
!! @param[inout] exner Exner pressure field
!! @param[in] theta Potential temperature field
!! @param[in] height_wt Height coordinate in wtheta
!! @param[in] height_w3 Height coordinate in w3
!! @param[in] chi_1 X component of the chi coordinate field
!! @param[in] chi_2 Y component of the chi coordinate field
!! @param[in] chi_3 Z component of the chi coordinate field
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] ndf_wt Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wt Number unique of degrees of freedom  for wtheta
!! @param[in] map_wt Dofmap for the cell at the base of the column for wt
!! @param[in] ndf_chi Number of degrees of freedom per cell for wchi
!! @param[in] undf_chi Number unique of degrees of freedom  for wchi
!! @param[in] map_chi Dofmap for the cell at the base of the column for wchi
!! @param[in] basis_chi Dofmap for the cell at the base of the column for wchi
subroutine hydrostatic_exner_code(nlayers, exner, theta, height_wt, height_w3,    &
                                  chi_1, chi_2, chi_3, ndf_w3, undf_w3, map_w3,   &
                                  ndf_wt, undf_wt, map_wt, ndf_chi, undf_chi, map_chi, basis_chi)

  use analytic_pressure_profiles_mod, only : analytic_pressure  

  implicit none   

  !Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_w3, ndf_wt, ndf_chi, undf_w3, undf_wt, undf_chi
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(undf_w3), intent(inout) :: exner
  real(kind=r_def), dimension(undf_wt), intent(in) :: theta, height_wt
  real(kind=r_def), dimension(undf_w3), intent(in) :: height_w3
  real(kind=r_def), dimension(undf_chi),intent(in) :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(1,ndf_chi,ndf_wt),  intent(in)  :: basis_chi

  !Internal variables
  integer(kind=i_def)                    :: k, dfc, layers_offset, wt_dof
  real(kind=r_def)                       :: dz  
  real(kind=r_def),   dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: x(3)

  real(kind=r_def)            :: exner_start
   
  if(init_exner_bt) then
    layers_offset = 0
    wt_dof = 1
  else
    layers_offset = nlayers - 1
    wt_dof = 2
  end if

  do dfc = 1, ndf_chi
    chi_1_e(dfc) = chi_1( map_chi(dfc) + layers_offset )
    chi_2_e(dfc) = chi_2( map_chi(dfc) + layers_offset )
    chi_3_e(dfc) = chi_3( map_chi(dfc) + layers_offset )
  end do

  ! Horizontal coordinates of cell bottom or top
  x(:) = 0.0_r_def
  do dfc = 1, ndf_chi
    x(1) = x(1) + chi_1_e(dfc)*basis_chi(1,dfc,wt_dof)
    x(2) = x(2) + chi_2_e(dfc)*basis_chi(1,dfc,wt_dof)
    x(3) = x(3) + chi_3_e(dfc)*basis_chi(1,dfc,wt_dof)
  end do
  
  ! Exner at the model surface or top
  exner_start = analytic_pressure( (/x(1), x(2), x(3)/), test, 0.0_r_def)

  if(init_exner_bt) then

    ! Bottom-up initialization
    ! Exner at the bottom level
    dz = height_w3(map_w3(1))-height_wt(map_wt(1))
    exner(map_w3(1)) = exner_start - gravity * dz /(cp*theta(map_wt(1)))

    ! Exner on other levels
    do k = 1, nlayers-1
      dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
      exner(map_w3(1)+k) = exner(map_w3(1)+k-1) - gravity * dz/(cp*theta(map_wt(1)+k))
    end do

  else

    ! Top-down initialization
    ! Exner at the top level
    dz = height_wt(map_wt(1)+nlayers) - height_w3(map_w3(1)+nlayers-1)
    exner(map_w3(1)+nlayers-1) = exner_start + gravity * dz /(cp*theta(map_wt(1)+nlayers))

    ! Exner on other levels
    do k = nlayers-1, 1, -1
      dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
      exner(map_w3(1)+k-1) = exner(map_w3(1)+k) + gravity * dz/(cp*theta(map_wt(1)+k))
    end do

  end if

end subroutine hydrostatic_exner_code

end module hydrostatic_exner_kernel_mod
