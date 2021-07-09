!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Apply 3D viscosity mu * (d2dx2 + d2dy2 + d2dz2) to the components
!>        of the momentum equation for lowest order elements.
!>
module momentum_viscosity_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_INC,           &
                                GH_SCALAR, ANY_SPACE_9,    &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                STENCIL, CROSS, CELL_COLUMN
  use chi_transform_mod, only : chi2xyz
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: momentum_viscosity_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2, STENCIL(CROSS)),        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: momentum_viscosity_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: momentum_viscosity_code

contains

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] u_inc Diffusion increment for wind field
!! @param[in] u_n Input wind field
!! @param[in] map_w2_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2 Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] panel_id Field describing the IDs of the mesh panels
!! @param[in] viscosity_mu Viscosity constant
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  Number of unique degrees of freedom  for wind_space
!! @param[in] cell_map_w2 Array holding the dofmap for the cell at the base of the column for w2
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom  for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi
!! @param[in] ndf_pid Number of degrees of freedom per cell for panel ID
!! @param[in] undf_pid Number of unique degrees of freedom for panel ID
!! @param[in] map_pid Dofmap for the cell at the base of the column for panel_id
subroutine momentum_viscosity_code(nlayers,                               &
                                   u_inc, u_n,                            &
                                   map_w2_size, map_w2,                   &
                                   chi1, chi2, chi3,                      &
                                   panel_id, viscosity_mu,                &
                                   ndf_w2, undf_w2, cell_map_w2,          &
                                   ndf_chi, undf_chi, map_chi,            &
                                   ndf_pid, undf_pid, map_pid             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: map_w2_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_size), intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_chi),            intent(in)  :: map_chi
  integer(kind=i_def), dimension(ndf_pid),            intent(in)  :: map_pid
  integer(kind=i_def), dimension(ndf_w2),             intent(in)  :: cell_map_w2

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  real(kind=r_def), intent(in) :: viscosity_mu

  ! Internal variables
  integer(kind=i_def)                      :: k, km, kp, df, ipanel
  real(kind=r_def)                         :: d2dx, d2dy, d2dz
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2, idz2
  real(kind=r_def), dimension(ndf_chi)     :: chi_x_e, chi_y_e, chi_z_e

  !  ----------
  !  |    |   |
  !  |  w | i |
  !  -----x----
  !  |    |   |
  !  | sw | s |
  !  ----------
  !  y
  !  ^
  !  |_> x
  !

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Compute grid spacing
  do k = 0, nlayers - 1
    do df = 1,ndf_chi
      call chi2xyz(chi1(map_chi(df)+k), chi2(map_chi(df)+k), chi3(map_chi(df)+k), &
                   ipanel, chi_x_e(df), chi_y_e(df), chi_z_e(df))
    end do
    idx2(k) = 1.0_r_def/(maxval(chi_x_e) - minval(chi_x_e))**2
    idy2(k) = 1.0_r_def/(maxval(chi_y_e) - minval(chi_y_e))**2
    idz2(k) = 1.0_r_def/(maxval(chi_z_e) - minval(chi_z_e))**2
  end do

  ! Velocity diffusion
  k = 0
  km = k
  kp = k+1
  ! u
  d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
  u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! v
  d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
  u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! w
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  do k = 1, nlayers-2
    km = k - 1
    kp = k + 1
    ! u
    d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
    u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! v
    d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
    u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! w
    d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
    u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
  end do

  k = nlayers-1
  km = k - 1
  kp = k
  ! u
  d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
  u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! v
  d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
  u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! w
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  k = nlayers
  km = k - 1
  kp = k
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(nlayers-1)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(nlayers-1)
  d2dz = (u_n(map_w2(5,1) + kp) - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(nlayers-1)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! Enforce zero flux boundary conditions
  u_inc(cell_map_w2(5) )             = 0.0_r_def
  u_inc(cell_map_w2(6) + nlayers-1 ) = 0.0_r_def

end subroutine momentum_viscosity_code

end module momentum_viscosity_kernel_mod
