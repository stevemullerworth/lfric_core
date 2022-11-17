!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the da_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module da_dev_driver_mod

  use checksum_alg_mod,         only: checksum_alg
  use cli_mod,                  only: get_initial_filename
  use configuration_mod,        only: final_configuration
  use constants_mod,            only: i_def, i_native, &
                                      PRECISION_REAL, r_def
  use clock_mod,                only: clock_type
  use driver_comm_mod,          only: init_comm, final_comm
  use driver_model_data_mod,    only: model_data_type
  use driver_log_mod,           only: init_logger, final_logger
  use driver_mesh_mod,          only: init_mesh, final_mesh
  use driver_fem_mod,           only: init_fem, final_fem
  use driver_io_mod,            only: init_io, final_io, filelist_populator, &
                                      get_io_context
  use io_context_mod,           only: io_context_type
  use field_mod,                only: field_type
  use da_dev_init_mod,          only: init_da_dev
  use log_mod,                  only: log_event,          &
                                      log_scratch_space,  &
                                      LOG_LEVEL_ALWAYS,   &
                                      LOG_LEVEL_INFO
  use mesh_mod,                 only: mesh_type
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: get_comm_size, &
                                      get_comm_rank
  use da_dev_mod,               only: load_configuration
  use da_dev_increment_alg_mod, only: da_dev_increment_alg
  use da_dev_init_files_mod,    only: init_da_dev_files
  use da_dev_config_mod,        only: pseudomodel, write_data, test_field
  use lfric_xios_read_mod,      only: read_state
  use lfric_xios_write_mod,     only: write_state
  use lfric_xios_context_mod,   only: lfric_xios_context_type

  implicit none

  private
  public initialise, run, finalise

  character(*), parameter :: program_name = "da_dev"

  type(model_data_type)  :: model_data
  type(model_clock_type) :: model_clock

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(mesh_type),  pointer              :: mesh      => null()
  type(mesh_type),  pointer              :: twod_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise()

    use step_calendar_mod,       only : step_calendar_type
    use time_config_mod,         only : timestep_start, timestep_end
    use timestepping_config_mod, only : dt, spinup_period

    implicit none

    character(:), allocatable              :: filename
    integer(i_native)                      :: model_communicator
    type(step_calendar_type)               :: calendar
    procedure(filelist_populator), pointer :: fl_populator => null()

    call init_comm( program_name, model_communicator )

    call get_initial_filename( filename )
    call load_configuration( filename, program_name )

    call init_logger(get_comm_size(), get_comm_rank(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh, &
                    twod_mesh = twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Initialise I/O context
    fl_populator => init_da_dev_files
    call init_io( program_name, model_communicator, chi, panel_id, populate_filelist=fl_populator )

    !Initialise model clock
    model_clock = model_clock_type( calendar%parse_instance(timestep_start), &
                                    calendar%parse_instance(timestep_end),   &
                                    dt, spinup_period )

    ! Initialise DA may eventually go here

    ! Create and initialise prognostic fields
    call init_da_dev( mesh, model_data )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()
    implicit none

    class(clock_type), pointer :: io_clock => null()
    class(io_context_type), pointer :: io_context => null()
    logical :: dummy

    io_context => get_io_context()
    select type (io_context)
    class is (lfric_xios_context_type)
      io_clock => io_context%get_clock()
      dummy = io_clock%tick()
    end select

    call model_clock%add_clock(io_clock)

    call log_event(program_name//": Run begins", LOG_LEVEL_INFO)

    do while( model_clock%tick() )
      call step()
    end do

    call log_event(program_name//": Run ends", LOG_LEVEL_INFO)
  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs a single time step.
  !>
  subroutine step()
    implicit none
    type(field_type), pointer  :: working_field => null()

    call model_data%depository%get_field( test_field, working_field )

    if ( pseudomodel ) then
      call read_state( model_data%depository, prefix="read_" )
    end if

    if( .not. pseudomodel ) then
      call da_dev_increment_alg( working_field )
    end if

    if ( write_data ) then
      call write_state( model_data%depository, prefix="write_" )
    end if

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    type(field_type), pointer :: working_field => null()

    call model_data%depository%get_field( test_field, working_field )

    !---------------------------------------------------------------------------
    ! Model finalise
    !---------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg( program_name, working_field, test_field )

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io()

    call final_fem()

    call final_mesh()

    call final_configuration()

    call final_logger( program_name )

    call final_comm()

  end subroutine finalise

end module da_dev_driver_mod
