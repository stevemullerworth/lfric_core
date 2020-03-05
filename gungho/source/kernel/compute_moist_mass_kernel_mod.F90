!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the cell integrated masses of a given moisture species.
!>
!> @details The kernel computes the cell integrated mass of a water species,
!> \f[ \int( \rho * m_X dV \f]
!>
module compute_moist_mass_kernel_mod

  use argument_mod,      only : arg_type, func_type,                  &
                                GH_FIELD, GH_WRITE, GH_READ,          &
                                ANY_SPACE_9,                          &
                                GH_BASIS, GH_DIFF_BASIS,              &
                                CELLS, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type
  use planet_config_mod, only : scaled_radius

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: compute_moist_mass_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                  &
         arg_type(GH_FIELD,   GH_WRITE, W3),                             &
         arg_type(GH_FIELD,   GH_READ,  Wtheta),                         &
         arg_type(GH_FIELD,   GH_READ,  W3),                             &
         arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                     &
         /)
    type(func_type) :: meta_funcs(3) = (/                                &
         func_type(W3,          GH_BASIS),                               &
         func_type(Wtheta,      GH_BASIS),                               &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
         /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_moist_mass_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_moist_mass_code
contains

!> @brief Compute the cell-integrated mass of a water species
!! @param[in] nlayers The number of layers
!! @param[out] water_mass The cell integrated mass of the water species
!! @param[in] mr_i The mixing ratio of the i-th moisture species
!! @param[in] rho The (dry) density
!! @param[in] chi_1 The physical x coordinate in chi
!! @param[in] chi_2 The physical y coordinate in chi
!! @param[in] chi_3 The physical z coordinate in chi
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis 4-dim array holding basis functions evaluated at
!!                     gaussian quadrature points
!! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
!! @param[in] undf_wtheta The number of unique degrees of freedom for wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for wtheta
!! @param[in] wtheta_basis 4-dim array holding basis functions evaluated at
!!                         quadrature points
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis 4-dim array holding differential of the basis
!!                           functions evaluated at gaussian quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_moist_mass_code(                                             &
                                    nlayers,                                    &
                                    water_mass,                                 &
                                    mr_i, rho,                                  &
                                    chi_1, chi_2, chi_3,                        &
                                    ndf_w3, undf_w3, map_w3, w3_basis,          &
                                    ndf_wtheta, undf_wtheta, map_wtheta,        &
                                    wtheta_basis,                               &
                                    ndf_chi, undf_chi, map_chi, chi_diff_basis, &
                                    nqp_h, nqp_v, wqp_h, wqp_v                  &
                                  )

  use coordinate_jacobian_mod,  only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_wtheta, ndf_chi
  integer(kind=i_def), intent(in) :: undf_w3, undf_wtheta, undf_chi
  integer(kind=i_def), dimension(ndf_w3), intent(in)     :: map_w3
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi), intent(in)    :: map_chi
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)     :: w3_basis
  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in)    :: chi_diff_basis
  real(kind=r_def), dimension(undf_w3),     intent(out)  :: water_mass
  real(kind=r_def), dimension(undf_w3),     intent(in)   :: rho
  real(kind=r_def), dimension(undf_wtheta), intent(in)   :: mr_i
  real(kind=r_def), dimension(undf_chi),    intent(in)   :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(nqp_h), intent(in)         :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)         :: wqp_v

  ! Internal variables
  integer(kind=i_def)                                    :: df, k, l
  integer(kind=i_def)                                    :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)                   :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)               :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)           :: jac
  real(kind=r_def), dimension(ndf_w3)                    :: rho_e, water_mass_e
  real(kind=r_def), dimension(ndf_wtheta)                :: mr_e
  real(kind=r_def)                                       :: mr_at_quad, rho_at_quad

  do k = 0, nlayers-1
    ! Extract coordinates for this element
    do df = 1, ndf_chi
      l = map_chi(df) + k
      chi_1_e(df) = chi_1(l)
      chi_2_e(df) = chi_2(l)
      chi_3_e(df) = chi_3(l)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)
    ! Loop through dofs, grabbing the values for this cell for reference element
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
      water_mass_e(df) = 0.0_r_def
    end do

    do df = 1, ndf_wtheta
      mr_e(df) = mr_i( map_wtheta(df) + k )
    end do

    ! Get values at quadrature points
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def
        mr_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad = rho_at_quad + rho_e(df) * w3_basis(1,df,qp1,qp2)
        end do
        do df = 1, ndf_wtheta
          mr_at_quad = mr_at_quad + mr_e(df) * wtheta_basis(1,df,qp1,qp2)
        end do
        ! Contribution to integral for water mass for this cell from quad point
        do df = 1, ndf_w3
          water_mass_e(df) = water_mass_e(df) + wqp_h(qp1) * wqp_v(qp2) &
            * rho_at_quad * mr_at_quad * dj(qp1,qp2)
        end do
      end do
    end do

    ! Return from reference element to physical element
    do df = 1, ndf_w3
      water_mass( map_w3(df) + k ) = water_mass_e(df)
    end do
  end do

end subroutine compute_moist_mass_code

end module compute_moist_mass_kernel_mod