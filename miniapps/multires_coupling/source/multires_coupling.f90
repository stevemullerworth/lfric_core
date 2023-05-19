!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp multires_coupling program

!> @brief Main program used for multires_coupling miniapp

!> @details Calls init, run and finalise routines from a driver module

program multires_coupling

  use cli_mod,                      only: get_initial_filename
  use driver_comm_mod,              only: init_comm, final_comm
  use driver_config_mod,            only: init_config, final_config
  use gungho_model_data_mod,        only: model_data_type
  use mpi_mod,                      only: global_mpi
  use multires_coupling_mod,        only: program_name, &
                                          multires_required_namelists
  use multires_coupling_driver_mod, only: initialise, run, finalise

  implicit none

  ! Model run working data set
  type (model_data_type) :: dynamics_mesh_model_data
  type (model_data_type) :: physics_mesh_model_data

  character(:), allocatable :: filename

  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, multires_required_namelists )
  deallocate( filename )

  call initialise( dynamics_mesh_model_data, &
                   physics_mesh_model_data,  &
                   global_mpi )

  call run( dynamics_mesh_model_data, physics_mesh_model_data )

  call finalise( dynamics_mesh_model_data, physics_mesh_model_data )
  call final_config()
  call final_comm()

end program multires_coupling
