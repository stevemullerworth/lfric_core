!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of IO within the
!! driver layer. This module also contains the model's io_context object and
!! associated getter routines
!>
module driver_io_mod

  use calendar_mod,            only: calendar_type
  use constants_mod,           only: i_native
  use driver_model_data_mod,   only: model_data_type
  use field_mod,               only: field_type
  use io_context_mod,          only: io_context_type
  use io_config_mod,           only: use_xios_io, subroutine_timers
  use log_mod,                 only: log_event, log_level_error, &
                                     log_level_trace
  use linked_list_mod,         only: linked_list_type
  use model_clock_mod,         only: model_clock_type
#ifdef USE_XIOS
  use lfric_xios_context_mod,  only: lfric_xios_context_type
#endif

  implicit none

  private
  public :: init_io, final_io,  &
            init_io_context,    &
            get_io_context,     &
            filelist_populator

  class(io_context_type), allocatable, target :: model_context

  abstract interface
    subroutine filelist_populator(files_list, model_data)
      import linked_list_type, model_data_type
      type(linked_list_type), intent(out) :: files_list
      class(model_data_type), optional, target, intent(in) :: model_data
    end subroutine filelist_populator
  end interface

contains

  !> @brief  Initialises the model I/O
  !>
  !> @param[in] id                A string identifier for the model
  !> @param[in] communicator      The ID for the model MPI communicator
  !> @param[in] chi               The model coordinate field
  !> @param[in] panel_id          Field containing the panel ID for each mesh
  !!                              vertex
  !> @param[in] model_clock       The model clock
  !> @param[in] calendar          The model calendar
  !> @param[in] populate_filelist Optional procedure for creating a list of
  !!                              file descriptions used by the model I/O
  !> @param[in] model_data        Optional Model data object
  !> @param[in] alt_coords        Optional array of coordinate fields
  !!                              for alternative meshes
  !> @param[in] alt_panel_ids     Optional panel ID fields for alternative meshes
  subroutine init_io( id, communicator,      &
                      chi, panel_id,         &
                      model_clock, calendar, &
                      populate_filelist,     &
                      model_data,            &
                      alt_coords, alt_panel_ids )

    implicit none

    character(*),                     intent(in)    :: id
    integer(i_native),                intent(in)    :: communicator
    class(field_type),                intent(in)    :: chi(:)
    class(field_type),                intent(in)    :: panel_id
    type(model_clock_type),           intent(inout) :: model_clock
    class(calendar_type),             intent(in)    :: calendar
    procedure(filelist_populator), &
                   pointer, optional, intent(in)    :: populate_filelist
    class(model_data_type), optional, intent(in)    :: model_data
    type(field_type),       optional, intent(in)    :: alt_coords(:,:)
    type(field_type),       optional, intent(in)    :: alt_panel_ids(:)

    ! Initialise model's main I/O context
    call init_io_context( model_context,         &
                          id, communicator,      &
                          chi, panel_id,         &
                          model_clock, calendar, &
                          populate_filelist,     &
                          model_data,            &
                          alt_coords,            &
                          alt_panel_ids )

  end subroutine init_io

  !> @brief  Initialises an I/O context based on user input
  !>
  !> @param[in] io_context        The I/O context to be set up
  !> @param[in] id                A string identifier for the model
  !> @param[in] communicator      The ID for the model MPI communicator
  !> @param[in] chi               The model coordinate field
  !> @param[in] panel_id          Field containing the panel ID for each mesh
  !!                              vertex
  !> @param[in] model_clock       The model clock
  !> @param[in] calendar          The model calendar
  !> @param[in] populate_filelist Optional procedure for creating a list of
  !!                              file descriptions used by the model I/O
  !> @param[in] model_data        Optional Model data object
  !> @param[in] alt_coords        Optional array of coordinate fields
  !!                              for alternative meshes
  !> @param[in] alt_panel_ids     Optional panel ID fields for alternative meshes
  subroutine init_io_context( io_context,            &
                              id, communicator,      &
                              chi, panel_id,         &
                              model_clock, calendar, &
                              populate_filelist,     &
                              model_data,            &
                              alt_coords,            &
                              alt_panel_ids )

    implicit none

    class(io_context_type), allocatable, intent(inout) :: io_context
    character(*),                        intent(in)    :: id
    integer(i_native),                   intent(in)    :: communicator
    class(field_type),                   intent(in)    :: chi(:)
    class(field_type),                   intent(in)    :: panel_id
    type(model_clock_type),              intent(inout) :: model_clock
    class(calendar_type),                intent(in)    :: calendar
    procedure(filelist_populator), &
                      pointer, optional, intent(in)    :: populate_filelist
    class(model_data_type),    optional, intent(in)    :: model_data
    type(field_type),          optional, intent(in)    :: alt_coords(:,:)
    type(field_type),          optional, intent(in)    :: alt_panel_ids(:)

    type(linked_list_type), pointer :: file_list
    integer(i_native) :: rc

    ! Allocate IO context type based on model configuration
    if ( use_xios_io ) then
#ifdef USE_XIOS
      allocate( lfric_xios_context_type::io_context, stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate LFRic-XIOS context object", &
                        log_level_error )
      end if

      ! Set up context
      select type(io_context)
      type is (lfric_xios_context_type)

        ! Populate list of I/O files if procedure passed through
        if (present(populate_filelist)) then
          file_list => io_context%get_filelist()
          call populate_filelist(file_list, model_data)
        end if
        call io_context%set_timer_flag(subroutine_timers)
        call io_context%initialise( id, communicator,      &
                                    chi, panel_id,         &
                                    model_clock, calendar, &
                                    alt_coords, alt_panel_ids )
      end select
#else
      call log_event( "Cannot use XIOS I/O: model has not been built with " // &
                      "enabled", log_level_error )
#endif
    else
      call log_event( "No I/O context to be set up for this run", &
                      log_level_trace )
      return
    end if

  end subroutine init_io_context

  !> @brief  Finalises the model I/O
  subroutine final_io()

    implicit none

    if (allocated(model_context)) deallocate(model_context)

  end subroutine final_io

  !> @brief  Returns the model io context.
  function get_io_context() result(context_ptr)

    implicit none

    class(io_context_type), pointer :: context_ptr

    context_ptr => model_context

  end function get_io_context

end module driver_io_mod
