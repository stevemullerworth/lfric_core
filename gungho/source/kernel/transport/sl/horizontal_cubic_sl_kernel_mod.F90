!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 using cubic
!!          semi-Lagrangian transport.
!> @details This kernel using cubic interpolation to solve the one-dimensional
!!          advection equation in both x and y, giving advective increments
!!          in both directions. This is the second part of the COSMIC splitting,
!!          so the x increment works on the field previously advected in the y-direction
!!          (and vice versa).
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_cubic_sl_kernel_mod

  use argument_mod,       only : arg_type,                   &
                                 GH_FIELD, GH_REAL,          &
                                 CELL_COLUMN, GH_WRITE,      &
                                 GH_READ, GH_SCALAR,         &
                                 STENCIL, CROSS, GH_INTEGER, &
                                 ANY_DISCONTINUOUS_SPACE_1,  &
                                 ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,      only : r_tran, i_def, r_def
  use transport_enumerated_types_mod, only : horizontal_monotone_strict
  use fs_continuity_mod,  only : W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! increment_x
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! increment_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field_y
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2 ),                & ! ix_start
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2 ),                & ! ix_end
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2 ),                & ! iy_start
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_2 ),                & ! iy_end
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts_y
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                        & ! monotone
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                        & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                                         & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: horizontal_cubic_sl_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] increment_x       Advective increment in x direction
  !> @param[in,out] increment_y       Advective increment in y direction
  !> @param[in]     field_x           Field from x direction
  !> @param[in]     stencil_size_x    Local length of field_x stencil
  !> @param[in]     stencil_map_x     Dofmap for the field_x stencil
  !> @param[in]     field_y           Field from y direction
  !> @param[in]     stencil_size_y    Local length of field_y stencil
  !> @param[in]     stencil_map_y     Dofmap for the field_y stencil
  !> @param[in]     ix_start          Start index in x for change in panel ID orientation
  !> @param[in]     ix_end            End index in x for change in panel ID orientation
  !> @param[in]     iy_start          Start index in y for change in panel ID orientation
  !> @param[in]     iy_end            End index in y for change in panel ID orientation
  !> @param[in]     dep_pts_x         Departure points in x
  !> @param[in]     dep_pts_y         Departure points in y
  !> @param[in]     monotone          Horizontal monotone option for cubic SL
  !> @param[in]     extent_size       Stencil extent needed for the LAM edge
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_wf            Number of degrees of freedom for field per cell
  !> @param[in]     undf_wf           Number of unique degrees of freedom for Wf
  !> @param[in]     map_wf            Map for Wf
  !> @param[in]     ndf_wp            Number of degrees of freedom for panel ID
  !!                                  index function space per cell
  !> @param[in]     undf_wp           Number of unique degrees of freedom for
  !!                                  panel ID index function space
  !> @param[in]     map_wp            Map for panel ID index function space
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Map for W2h

  subroutine horizontal_cubic_sl_code( nlayers,        &
                                       increment_x,    &
                                       increment_y,    &
                                       field_x,        &
                                       stencil_size_x, &
                                       stencil_map_x,  &
                                       field_y,        &
                                       stencil_size_y, &
                                       stencil_map_y,  &
                                       ix_start,       &
                                       ix_end,         &
                                       iy_start,       &
                                       iy_end,         &
                                       dep_pts_x,      &
                                       dep_pts_y,      &
                                       monotone,       &
                                       extent_size,    &
                                       dt,             &
                                       ndf_wf,         &
                                       undf_wf,        &
                                       map_wf,         &
                                       ndf_wp,         &
                                       undf_wp,        &
                                       map_wp,         &
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
    integer(kind=i_def), intent(in) :: undf_wp
    integer(kind=i_def), intent(in) :: ndf_wp
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_y
    integer(kind=i_def), intent(in) :: monotone
    integer(kind=i_def), intent(in) :: extent_size
    real(kind=r_tran),   intent(in) :: dt

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_wf),  intent(in) :: map_wf
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_wp),  intent(in) :: map_wp
    integer(kind=i_def), dimension(ndf_wf,stencil_size_x), intent(in) :: stencil_map_x
    integer(kind=i_def), dimension(ndf_wf,stencil_size_y), intent(in) :: stencil_map_y

    ! Arguments: Fields
    real(kind=r_tran),   dimension(undf_wf),  intent(inout) :: increment_x
    real(kind=r_tran),   dimension(undf_wf),  intent(inout) :: increment_y
    real(kind=r_tran),   dimension(undf_wf),  intent(in)    :: field_x
    real(kind=r_tran),   dimension(undf_wf),  intent(in)    :: field_y
    integer(kind=i_def), dimension(undf_wp),  intent(in)    :: ix_start
    integer(kind=i_def), dimension(undf_wp),  intent(in)    :: ix_end
    integer(kind=i_def), dimension(undf_wp),  intent(in)    :: iy_start
    integer(kind=i_def), dimension(undf_wp),  intent(in)    :: iy_end
    real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: dep_pts_x
    real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: dep_pts_y

    ! Variables for flux calculation
    real(kind=r_tran) :: departure_dist, departure_dist_w3, departure_dist_wt

    ! Local fields
    real(kind=r_tran)   :: field_local_for_x(1:(stencil_size_x+ 1) / 2)
    real(kind=r_tran)   :: field_local_for_y(1:(stencil_size_y+ 1) / 2)
    real(kind=r_tran)   :: field_x_local(1:(stencil_size_x+ 1) / 2)
    real(kind=r_tran)   :: field_y_local(1:(stencil_size_y+ 1) / 2)

    ! Cubic interpolation coefficients
    real(kind=r_tran)   :: frac_d, x0, x1, x2, x3, xx, q0, q1, q2, q3
    real(kind=r_tran)   :: field_out_x, field_out_y
    real(kind=r_tran)   :: qx_max, qx_min, qy_max, qy_min

    ! Indices
    integer(kind=i_def) :: ind_lo, ind_hi, k_w2h
    integer(kind=i_def) :: k, km1, kp1, jj, ijk, int_d, nl

    ! Stencils
    integer(kind=i_def) :: stencil_half, stencil_size, lam_edge_size

    ! Stencil has order e.g.
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
    ! Increments calculated for centre cell, e.g. cell 1 for stencil, cell 5 for local

    ! nl = nlayers-1  for w3
    !    = nlayers    for wtheta
    nl = nlayers - 1 + (ndf_wf - 1)

    ! Use stencil_size_y as each stencil size should be equal - as stencil_size_y
    ! is a cross stencil we need 1D stencil size from this
    stencil_size = (stencil_size_y + 1_i_def) / 2_i_def
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 4_i_def*extent_size+1_i_def

    if (lam_edge_size > stencil_size_x) then

      ! At edge of LAM, so set output to zero
      do k = 0, nl
        increment_x( map_wf(1) + k ) = 0.0_r_tran
        increment_y( map_wf(1) + k ) = 0.0_r_tran
      end do

    else

      ! Not at edge of LAM so compute increments

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

        ! Get local field values - this will depend on panel ID at cubed sphere edges
        do jj = 1, stencil_half
          field_y_local(jj) = field_y(stencil_map_y(1,stencil_half+1-jj) + k)
          field_x_local(jj) = field_x(stencil_map_x(1,stencil_half+1-jj) + k)
        end do
        do jj = stencil_half+1, stencil_size
          field_y_local(jj) = field_y(stencil_map_y(1,stencil_half-1+jj) + k)
          field_x_local(jj) = field_x(stencil_map_x(1,stencil_half-1+jj) + k)
        end do
        field_local_for_x(:) = field_y_local(:)
        do ijk = iy_start(map_wp(1)), iy_end(map_wp(1))
          field_local_for_x(ijk) = field_x_local(ijk)
        end do

        ! Set up cubic interpolation at correct index
        x0 = 0.0_r_tran
        x1 = 1.0_r_tran
        x2 = 2.0_r_tran
        x3 = 3.0_r_tran
        if (departure_dist >= 0.0_r_tran) then
          ind_hi = stencil_half - int_d
          ind_lo = ind_hi - 1
          xx = x2 - frac_d
        else
          ind_hi = stencil_half - int_d + 1
          ind_lo = ind_hi - 1
          xx = x1 - frac_d
        end if
        q0 = field_local_for_x(ind_lo-1)
        q1 = field_local_for_x(ind_lo)
        q2 = field_local_for_x(ind_hi)
        q3 = field_local_for_x(ind_hi+1)
        qx_max = max(q0,q1,q2,q3)
        qx_min = min(q0,q1,q2,q3)

        field_out_x = ((xx-x1)*(xx-x2)*(xx-x3))/((x0-x1)*(x0-x2)*(x0-x3))*q0 + &
                      ((xx-x0)*(xx-x2)*(xx-x3))/((x1-x0)*(x1-x2)*(x1-x3))*q1 + &
                      ((xx-x0)*(xx-x1)*(xx-x3))/((x2-x0)*(x2-x1)*(x2-x3))*q2 + &
                      ((xx-x0)*(xx-x1)*(xx-x2))/((x3-x0)*(x3-x1)*(x3-x2))*q3

        ! y direction departure distance at centre
        departure_dist_w3 = ( dep_pts_y( map_w2h(2) + k_w2h)+dep_pts_y( map_w2h(4) + k_w2h) )/2.0_r_tran
        departure_dist_wt = ( dep_pts_y( map_w2h(2) + km1)+dep_pts_y( map_w2h(4) + km1) + &
                              dep_pts_y( map_w2h(2) + kp1)+dep_pts_y( map_w2h(4) + kp1) )/4.0_r_tran
        departure_dist = (2 - ndf_wf) * departure_dist_w3 + (ndf_wf - 1) * departure_dist_wt

        ! Calculates number of cells of interest to move
        frac_d = departure_dist - int(departure_dist)
        int_d = int(departure_dist,i_def)

        ! Get local field values - this will depend on panel ID at cubed sphere edges
        field_y_local(stencil_half) = field_y(stencil_map_y(1,1)+k)
        field_x_local(stencil_half) = field_x(stencil_map_x(1,1)+k)
        do jj = 1, stencil_half-1
          field_y_local(jj) = field_y(stencil_map_y(1,stencil_size+1-jj) + k)
          field_x_local(jj) = field_x(stencil_map_x(1,stencil_size+1-jj) + k)
        end do
        do jj = stencil_half+1, stencil_size
          field_y_local(jj) = field_y(stencil_map_y(1,stencil_size-1+jj) + k)
          field_x_local(jj) = field_x(stencil_map_x(1,stencil_size-1+jj) + k)
        end do
        field_local_for_y(:) = field_x_local(:)
        do ijk = ix_start(map_wp(1)), ix_end(map_wp(1))
          field_local_for_y(ijk) = field_y_local(ijk)
        end do

        ! Set up cubic interpolation at correct index
        x0 = 0.0_r_tran
        x1 = 1.0_r_tran
        x2 = 2.0_r_tran
        x3 = 3.0_r_tran
        if (departure_dist >= 0.0_r_tran) then
          ind_hi = stencil_half - int_d
          ind_lo = ind_hi - 1
          xx = x2 - frac_d
        else
          ind_hi = stencil_half - int_d + 1
          ind_lo = ind_hi - 1
          xx = x1 - frac_d
        end if
        q0 = field_local_for_y(ind_lo-1)
        q1 = field_local_for_y(ind_lo)
        q2 = field_local_for_y(ind_hi)
        q3 = field_local_for_y(ind_hi+1)
        qy_max = max(q0,q1,q2,q3)
        qy_min = min(q0,q1,q2,q3)

        field_out_y = ((xx-x1)*(xx-x2)*(xx-x3))/((x0-x1)*(x0-x2)*(x0-x3))*q0 + &
                      ((xx-x0)*(xx-x2)*(xx-x3))/((x1-x0)*(x1-x2)*(x1-x3))*q1 + &
                      ((xx-x0)*(xx-x1)*(xx-x3))/((x2-x0)*(x2-x1)*(x2-x3))*q2 + &
                      ((xx-x0)*(xx-x1)*(xx-x2))/((x3-x0)*(x3-x1)*(x3-x2))*q3

        ! Monotone
        if (monotone == horizontal_monotone_strict) then
          if (field_out_x > qx_max) then
            field_out_x = qx_max
          else if (field_out_x < qx_min) then
            field_out_x = qx_min
          end if
          if (field_out_y > qy_max) then
            field_out_y = qy_max
          else if (field_out_y < qy_min) then
            field_out_y = qy_min
          end if
        end if

        ! Get increment
        increment_x(map_wf(1)+k) = (field_local_for_x(stencil_half) - field_out_x) / dt
        increment_y(map_wf(1)+k) = (field_local_for_y(stencil_half) - field_out_y) / dt

      end do ! vertical levels k

    end if

  end subroutine horizontal_cubic_sl_code

end module horizontal_cubic_sl_kernel_mod
