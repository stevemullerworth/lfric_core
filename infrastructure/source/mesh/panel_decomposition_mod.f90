!-----------------------------------------------------------------------------
! Copyright (c) 2025,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Provides strategies for decomposing rectangular panels into partitions
!>
module panel_decomposition_mod

  use constants_mod, only: i_def, l_def, r_def
  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  ! Parent for panel decomposition types
  type, public, abstract :: panel_decomposition_type
  contains
    procedure(get_partition_interface), deferred :: get_partition
  end type panel_decomposition_type

  ! Decomposition that accepts user specified number of xprocs and yprocs
  type, extends(panel_decomposition_type), public :: custom_decomposition_type
    integer(i_def) :: num_xprocs, num_yprocs
  contains
    procedure, public :: get_partition => get_custom_partition
  end type custom_decomposition_type
  ! Constructor
  interface custom_decomposition_type
    module procedure custom_decomposition_constructor
  end interface

  ! Decomposition that automatically determines number of xprocs and yprocs
  type, extends(panel_decomposition_type), public :: auto_decomposition_type
  contains
    procedure, public :: get_partition => get_auto_partition
  end type auto_decomposition_type

  ! Decomposition only in x direction
  type, extends(panel_decomposition_type), public :: row_decomposition_type
  contains
    procedure, public :: get_partition => get_row_partition
  end type row_decomposition_type

  ! Decomposition only in y direction
  type, extends(panel_decomposition_type), public :: column_decomposition_type
  contains
    procedure, public :: get_partition => get_column_partition
  end type column_decomposition_type


  ! Interface for routines that generate partition shape and location
  abstract interface

    subroutine get_partition_interface( self,             &
                                        relative_rank,    &
                                        panel_ranks,      &
                                        num_cells_x,      &
                                        num_cells_y,      &
                                        any_maps,         &
                                        partition_width,  &
                                        partition_height, &
                                        partition_x_pos,  &
                                        partition_y_pos )
      use constants_mod, only: i_def
      import :: panel_decomposition_type

      class(panel_decomposition_type), intent(in) :: self

      integer(i_def), intent(in)    :: relative_rank, &
                                       panel_ranks,   &
                                       num_cells_x,   &
                                       num_cells_y
      logical,        intent(in)    :: any_maps

      integer(i_def), intent(inout) :: partition_width,  &
                                       partition_height, &
                                       partition_x_pos,  &
                                       partition_y_pos

    end subroutine get_partition_interface

  end interface

contains

  !> @brief Partition the panel into a given number of x and y processes
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_custom_partition( self,             &
                                   relative_rank,    &
                                   panel_ranks,      &
                                   num_cells_x,      &
                                   num_cells_y,      &
                                   any_maps,         &
                                   partition_width,  &
                                   partition_height, &
                                   partition_x_pos,  &
                                   partition_y_pos )
    implicit none

    class(custom_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    num_xprocs = self%num_xprocs
    num_yprocs = self%num_yprocs

    if ( panel_ranks /= num_xprocs * num_yprocs ) then
      write( log_scratch_space, "(a,i0,a,i0,a,i0)" ) "Total ranks per panel ", panel_ranks, " must be the product of xprocs ", self%num_xprocs, " and yprocs ", self%num_yprocs
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if


    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

  end subroutine get_custom_partition


  !> @brief Constructor for custom_decomposition_type
  !> @param[in] xprocs The requested number of partitions in the x direction
  !> @param[in] yprocs The requested number of partitions in the y direction
  function custom_decomposition_constructor(xprocs, yprocs) result(self)
    implicit none

    type(custom_decomposition_type), target :: self
    integer(i_def) :: xprocs, yprocs

    self%num_xprocs = xprocs
    self%num_yprocs = yprocs

  end function custom_decomposition_constructor


  !> @brief Partition the panel into an automatically determined number of x and
  !         y processes
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_auto_partition( self,             &
                                 relative_rank,    &
                                 panel_ranks,      &
                                 num_cells_x,      &
                                 num_cells_y,      &
                                 any_maps,         &
                                 partition_width,  &
                                 partition_height, &
                                 partition_x_pos,  &
                                 partition_y_pos )
    implicit none

    class(auto_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    integer(i_def) :: start_factor
    integer(i_def) :: end_factor
    integer(i_def) :: fact_count
    logical(l_def) :: found_factors

    integer(i_def), parameter :: max_factor_iters = 10000

    ! For automatic partitioning, try to partition into the squarest
    ! possible partitions by finding the two factors of panel_ranks
    ! that are closest to sqrt(panel_ranks). If two factors can't
    ! be found after max_factor_iters attempts, they would provide
    ! partitions that are too un-square, so an error is produced.
    start_factor = nint(sqrt(real(panel_ranks, kind=r_def)), kind=i_def)
    end_factor = max(1,(start_factor-max_factor_iters))
    found_factors = .false.
    do fact_count = start_factor, end_factor, -1
      if (mod(panel_ranks, fact_count) == 0) then
        found_factors = .true.
        exit
      end if
    end do

    if (.not. found_factors) then
      call log_event( "Could not automatically partition domain.", &
                      LOG_LEVEL_ERROR )
    end if

    num_xprocs = fact_count
    num_yprocs = panel_ranks / fact_count

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

  end subroutine get_auto_partition


  !> @brief Partition the panel only in the x direction
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_row_partition( self,             &
                                relative_rank,    &
                                panel_ranks,      &
                                num_cells_x,      &
                                num_cells_y,      &
                                any_maps,         &
                                partition_width,  &
                                partition_height, &
                                partition_x_pos,  &
                                partition_y_pos )
    implicit none

    class(row_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    num_xprocs = panel_ranks
    num_yprocs = 1_i_def

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

  end subroutine get_row_partition


  !> @brief Partition the panel only in the y direction
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    any_maps         Whether there exist maps between meshes that
  !>                                must having aligning partitions
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine get_column_partition( self,             &
                                   relative_rank,    &
                                   panel_ranks,      &
                                   num_cells_x,      &
                                   num_cells_y,      &
                                   any_maps,         &
                                   partition_width,  &
                                   partition_height, &
                                   partition_x_pos,  &
                                   partition_y_pos )
    implicit none

    class(column_decomposition_type), intent(in) :: self
    integer(i_def), intent(in)    :: relative_rank,    &
                                     panel_ranks,      &
                                     num_cells_x,      &
                                     num_cells_y
    logical,        intent(in)    :: any_maps
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: num_xprocs, num_yprocs

    num_xprocs = 1_i_def
    num_yprocs = panel_ranks

    call xy_defensive_checks( num_cells_x, &
                              num_cells_y, &
                              num_xprocs,  &
                              num_yprocs,  &
                              panel_ranks, &
                              any_maps )

    call xy_decomposition( relative_rank,    &
                           num_cells_x,      &
                           num_cells_y,      &
                           num_xprocs,       &
                           num_yprocs,       &
                           partition_width,  &
                           partition_height, &
                           partition_x_pos,  &
                           partition_y_pos )

  end subroutine get_column_partition


  !> @brief Helper function for generating identical partitions arranged in a
  !         rectangular grid
  !> @param[in]    relative_rank    The number of this rank in the order of all
  !                                 ranks on the panel
  !> @param[in]    panel_ranks      The total number of ranks on the panel
  !> @param[in]    num_cells_x      The panel's size in the x direction
  !> @param[in]    num_cells_y      The panel's size in the y direction
  !> @param[in]    num_xprocs       The number of partitions in the x direction
  !> @param[in]    num_yprocs       The number of partitions in the y direction
  !> @param[inout] partition_width  The partition's size in the x direction
  !> @param[inout] partition_height The partition's size in the y direction
  !> @param[inout] partition_x_pos  The x index of the partition
  !> @param[inout] partition_y_pos  The y index of the partition
  subroutine xy_decomposition( relative_rank,    &
                               num_cells_x,      &
                               num_cells_y,      &
                               num_xprocs,       &
                               num_yprocs,       &
                               partition_width,  &
                               partition_height, &
                               partition_x_pos,  &
                               partition_y_pos )
    implicit none

    integer(i_def), intent(in)    :: relative_rank,    &
                                     num_cells_x,      &
                                     num_cells_y,      &
                                     num_xprocs,       &
                                     num_yprocs
    integer(i_def), intent(inout) :: partition_width,  &
                                     partition_height, &
                                     partition_x_pos,  &
                                     partition_y_pos

    integer(i_def) :: xproc, yproc

    xproc = mod(relative_rank - 1, num_xprocs)
    yproc = (relative_rank - 1) / num_xprocs

    partition_x_pos  = ( (xproc * num_cells_x) / num_xprocs ) + 1
    partition_width  = ( ((xproc+1)*num_cells_x) / num_xprocs ) - partition_x_pos + 1
    partition_y_pos  = ( (yproc * num_cells_y) / num_yprocs ) + 1
    partition_height = ( ((yproc+1)*num_cells_y) / num_yprocs ) - partition_y_pos + 1

  end subroutine xy_decomposition


  !> @brief Defensive checks common for all decomposition strategies into a
  !         rectangular grid
  !> @param[in] num_cells_x The panel's size in the x direction
  !> @param[in] num_cells_y The panel's size in the y direction
  !> @param[in] num_xprocs  The number of partitions in the x direction
  !> @param[in] num_yprocs  The number of partitions in the y direction
  subroutine xy_defensive_checks( num_cells_x, &
                                  num_cells_y, &
                                  num_xprocs,  &
                                  num_yprocs,  &
                                  panel_ranks, &
                                  any_maps )
    implicit none

    integer(i_def), intent(in) :: num_cells_x, &
                                  num_cells_y, &
                                  num_xprocs,  &
                                  num_yprocs,  &
                                  panel_ranks
    logical       , intent(in) :: any_maps

    if ( num_xprocs <=0 .or. num_yprocs <= 0 ) then
      call log_event("Number of x and y processes must be strictly positive.", LOG_LEVEL_ERROR)
    end if

    if ( num_cells_x < num_xprocs .or. num_cells_y < num_yprocs ) then
      write(log_scratch_space, "(a,i0,a,i0)") "Must have more cells than partitions in both x and y directions."
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Equal divisions are only required if there are maps between meshes
    if (any_maps) then
      if ( mod(num_cells_x, num_xprocs) /= 0 ) then
        write(log_scratch_space, "(a,i0,a,i0)") "Requested number of ranks in x direction ", num_xprocs, &
          " must divide panel x dimension ", num_cells_x
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if

      if ( mod(num_cells_y, num_yprocs) /= 0 ) then
        write(log_scratch_space, "(a,i0,a,i0)") "Requested number of ranks in y direction ", num_yprocs, &
          " must divide panel y dimension ", num_cells_y
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end if

    if ( num_xprocs * num_yprocs /= panel_ranks ) then
      write(log_scratch_space, "(a,i0,a,i0)") "Requested number of partitions ", num_xprocs * num_yprocs, &
        " must equal available number of ranks per panel ", panel_ranks
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end subroutine xy_defensive_checks

end module panel_decomposition_mod