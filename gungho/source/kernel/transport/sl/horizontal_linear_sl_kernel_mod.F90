!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Calculates the fields in x and y at time n+1 using linear
!!          semi-Lagrangian transport.
!> @details This kernel using linear interpolation to solve the one-dimensional
!!          advection equation in both x and y.
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_linear_sl_kernel_mod

  use argument_mod,       only : arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 CELL_COLUMN, GH_WRITE,     &
                                 GH_READ, GH_SCALAR,        &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 STENCIL, CROSS, GH_INTEGER
  use constants_mod,      only : i_def, r_tran
  use fs_continuity_mod,  only : W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_linear_sl_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! field_out_x
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! field_out_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts_y
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     )                                         & ! extent_size
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_linear_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: horizontal_linear_sl_code

contains

  !> @brief Compute the advective field in x and y using linear interpolation.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] field_x        Field at time n+1 in x direction
  !> @param[in,out] field_y        Field at time n+1 in y direction
  !> @param[in]     field          Field to transport
  !> @param[in]     stencil_size_c Local length of field cross stencil
  !> @param[in]     stencil_map    Dofmap for the field stencil
  !> @param[in]     dep_pts_x      Departure points in x
  !> @param[in]     dep_pts_y      Departure points in y
  !> @param[in]     extent_size    Stencil extent needed for the LAM edge
  !> @param[in]     ndf_wf         Number of degrees of freedom for field per cell
  !> @param[in]     undf_wf        Number of unique degrees of freedom for field
  !> @param[in]     map_wf         Map for field
  !> @param[in]     ndf_w2h        Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h       Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h        Map for W2h

  subroutine horizontal_linear_sl_code( nlayers,        &
                                        field_x,        &
                                        field_y,        &
                                        field,          &
                                        stencil_size_c, &
                                        stencil_map,    &
                                        dep_pts_x,      &
                                        dep_pts_y,      &
                                        extent_size,    &
                                        ndf_wf,         &
                                        undf_wf,        &
                                        map_wf,         &
                                        ndf_w2h,        &
                                        undf_w2h,       &
                                        map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_size_c

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_wf),  intent(in) :: map_wf
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_wf,stencil_size_c), intent(in) :: stencil_map

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_wf),  intent(inout) :: field_x
    real(kind=r_tran), dimension(undf_wf),  intent(inout) :: field_y
    real(kind=r_tran), dimension(undf_wf),  intent(in)    :: field
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: dep_pts_x
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: dep_pts_y
    integer(kind=i_def), intent(in)                       :: extent_size

    ! Local fields
    real(kind=r_tran)   :: field_local(1:(stencil_size_c+ 1) / 2)
    real(kind=r_tran)   :: departure_dist, departure_dist_w3, departure_dist_wt

    ! Interpolation coefficients
    real(kind=r_tran)   :: frac_d, x0, x1, xx, q0, q1

    ! Indices
    integer(kind=i_def) :: ind_lo, ind_hi, nl, k_w2h
    integer(kind=i_def) :: k, km1, kp1, jj, int_d

    ! Stencils
    integer(kind=i_def) :: stencil_size, stencil_half, lam_edge_size

    ! Cross stencil has order e.g.
    !                           | 17 |
    !                           | 16 |
    !                           | 15 |
    !                           | 14 |
    !       |  5 |  4 |  3 |  2 |  1 | 10 | 11 | 12 | 13 | for extent 4
    !                           |  6 |
    !                           |  7 |
    !                           |  8 |
    !                           |  9 |
    !
    ! Local fields have order e.g.  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | for extent 4
    ! in x and y
    ! Advection calculated for centre cell, e.g. cell 1 for stencil, cell 5 for local

    ! nl = nlayers-1  for w3
    !    = nlayers    for wtheta
    nl = nlayers - 1 + (ndf_wf - 1)

    ! stencil is a cross stencil we need 1D stencil size from this
    stencil_size = (stencil_size_c + 1_i_def) / 2_i_def
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 4_i_def*extent_size+1_i_def

    if ( lam_edge_size > stencil_size_c) then

      ! At edge of LAM, so set output to zero
      do k = 0, nl
         field_x( map_wf(1) + k ) = 0.0_r_tran
         field_y( map_wf(1) + k ) = 0.0_r_tran
      end do

    else

      ! Loop over k levels
      do k = 0, nl

        k_w2h = min(k,nlayers-1)
        km1 = max(k-1,0)
        kp1 = min(k,nlayers-1)

        ! x direction departure distance at centre
        departure_dist_w3 = ( dep_pts_x( map_w2h(1) + k_w2h)+dep_pts_x( map_w2h(3) + k_w2h) )/2.0_r_tran
        departure_dist_wt = ( dep_pts_x( map_w2h(1) + km1)+dep_pts_x( map_w2h(3) + km1) + &
                              dep_pts_x( map_w2h(1) + kp1)+dep_pts_x( map_w2h(3) + kp1) )/4.0_r_tran
        departure_dist = (2 - ndf_wf) * departure_dist_w3 + (ndf_wf - 1) * departure_dist_wt

        ! Calculates number of cells of interest to move
        frac_d = departure_dist - int(departure_dist)
        int_d = int(departure_dist,i_def)

        ! Get local field values from cross stencil
        do jj = 1, stencil_half
          field_local(jj) = field(stencil_map(1,stencil_half+1-jj) + k)
        end do
        do jj = stencil_half+1, stencil_size
          field_local(jj) = field(stencil_map(1,stencil_half-1+jj) + k)
        end do

        ! Set up linear interpolation in correct cell
        x0 = 0.0_r_tran
        x1 = 1.0_r_tran
        if (departure_dist >= 0.0_r_tran) then
          ind_hi = stencil_half - int_d
          ind_lo = ind_hi - 1
          xx = x1 - frac_d
        else
          ind_hi = stencil_half - int_d + 1
          ind_lo = ind_hi - 1
          xx = x0 - frac_d
        end if

        q0 = field_local(ind_lo)
        q1 = field_local(ind_hi)

        field_x(map_wf(1)+k) = (xx-x1)/(x0-x1)*q0 + (xx-x0)/(x1-x0)*q1

        ! y direction departure distance at centre
        departure_dist_w3 = ( dep_pts_y( map_w2h(2) + k_w2h)+dep_pts_y( map_w2h(4) + k_w2h) )/2.0_r_tran
        departure_dist_wt = ( dep_pts_y( map_w2h(2) + km1)+dep_pts_y( map_w2h(4) + km1) + &
                              dep_pts_y( map_w2h(2) + kp1)+dep_pts_y( map_w2h(4) + kp1) )/4.0_r_tran
        departure_dist = (2 - ndf_wf) * departure_dist_w3 + (ndf_wf - 1) * departure_dist_wt

        ! Calculates number of cells of interest to move
        frac_d = departure_dist - int(departure_dist)
        int_d = int(departure_dist,i_def)

        ! Get local field values from cross stencil
        field_local(stencil_half) = field(stencil_map(1,1) + k)
        do jj = 1, stencil_half-1
          field_local(jj) = field(stencil_map(1,stencil_size+1-jj) + k)
        end do
        do jj = stencil_half+1, stencil_size
          field_local(jj) = field(stencil_map(1,stencil_size-1+jj) + k)
        end do

        ! Set up linear interpolation in correct cell
        x0 = 0.0_r_tran
        x1 = 1.0_r_tran
        if (departure_dist >= 0.0_r_tran) then
          ind_hi = stencil_half - int_d
          ind_lo = ind_hi - 1
          xx = x1 - frac_d
        else
          ind_hi = stencil_half - int_d + 1
          ind_lo = ind_hi - 1
          xx = x0 - frac_d
        end if
        q0 = field_local(ind_lo)
        q1 = field_local(ind_hi)

        field_y(map_wf(1)+k) = (xx-x1)/(x0-x1)*q0 + (xx-x0)/(x1-x0)*q1

      end do ! vertical levels k

    end if

  end subroutine horizontal_linear_sl_code

end module horizontal_linear_sl_kernel_mod
