!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @mainpage miniapp_skeleton
!> Barebones miniapp that can be taken and adapted for new science development

program miniapp_skeleton

  use constants_mod,                  only : i_def
  use cli_mod,                        only : get_initial_filename
  use miniapp_skeleton_mod,           only : load_configuration
  use init_gungho_mod,                only : init_gungho
  use init_miniapp_skeleton_mod,      only : init_miniapp_skeleton
  use ESMF
  use field_mod,                      only : field_type
  use miniapp_skeleton_alg_mod,       only : miniapp_skeleton_alg
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO
  use output_alg_mod,                 only : output_alg
  use checksum_alg_mod,               only : checksum_alg

  implicit none

  character(:), allocatable :: filename

  type(ESMF_VM)      :: vm
  integer            :: rc
  integer            :: total_ranks, local_rank
  integer            :: petCount, localPET

  integer            :: mesh_id

  ! prognostic fields
  type( field_type ) :: field_1
  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise ESMF and get the rank information from the virtual machine
  CALL ESMF_Initialize(vm=vm, defaultlogfilename="miniapp_skeleton.Log", &
                  logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  total_ranks = petCount
  local_rank  = localPET

  call log_event( 'skeleton miniapp running...', LOG_LEVEL_INFO )

  call get_initial_filename( filename )
  call load_configuration( filename )
  deallocate( filename )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------
  ! Create the mesh and function space collection
  call init_gungho(mesh_id, local_rank, total_ranks)

  ! Create and initialise prognostic fields
  call init_miniapp_skeleton(mesh_id, field_1)

  ! Call an algorithm
  call miniapp_skeleton_alg(field_1)
  

  ! Write out output file
  call log_event("skeleton miniapp: writing diagnostic output", LOG_LEVEL_INFO)
  call output_alg('skeleton_field', 0, field_1, mesh_id)
  
  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Write checksums to file
  call checksum_alg('miniapp_skeleton', field_1, 'skeleton_field_1')

  call log_event( 'Skeleton miniapp completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Close down ESMF
  call ESMF_Finalize(rc=rc)

end program miniapp_skeleton
