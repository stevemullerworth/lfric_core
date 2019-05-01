!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create empty fields for use by the gungho model
!> @details Creates the empty prognostic fields that will be later 
!>          initialise and used by the gungho model

module create_gungho_prognostics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, &
                                             write_diag_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use finite_element_config_mod,      only : element_order, &
                                             vorticity_in_w1
  use formulation_config_mod,         only : use_moisture
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use function_space_collection_mod , only : function_space_collection
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use mr_indices_mod,                 only : nummr, &
                                             mr_names
  use moist_dyn_mod,                  only : num_moist_factors
  use io_mod,                         only : xios_write_field_node, &
                                             xios_write_field_face, &
                                             checkpoint_read_xios,   &
                                             checkpoint_read_netcdf, &
                                             checkpoint_write_netcdf,&
                                             checkpoint_write_xios
  use io_config_mod,                  only : use_xios_io,     &
                                             write_diag,      &
                                             checkpoint_read, &
                                             checkpoint_write
  implicit none

  private
  public :: create_gungho_prognostics

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Create empty fields to be used as prognostics by the gungho model
  !> @param[in] mesh_id The identifier given to the current 3d mesh
  !> @param[in] chi. An array of fields that hold the finite-element
  !>                 representation of the mesh
  !> @param[inout] prognostic_fields A collection of the fields that make up the
  !>                                 prognostic variables in the model 
  !> @param[inout] mr An array of fields that hold the moisture mixing ratios
  !> @param[inout] moist_dyn An array of the moist dynamics fields
  subroutine create_gungho_prognostics(mesh_id, chi, prognostic_fields, mr, &
                                       moist_dyn, xi)
    implicit none

    integer(i_def), intent(in)                 :: mesh_id
    ! Prognostic fields
    type(field_collection_type), intent(inout) :: prognostic_fields
    type( field_type ), intent(inout), target  :: mr(nummr)
    ! Diagnostic fields
    type( field_type ), intent(inout)        :: xi
    type( field_type ), intent(inout)        :: moist_dyn(num_moist_factors)
    ! Coordinate fields
    type( field_type ), intent(in)           :: chi(:)

    type( field_type ), pointer              :: tmp_ptr

    integer(i_def)                           :: imr

    procedure(write_diag_interface), pointer :: tmp_write_diag_ptr
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr
    procedure(checkpoint_read_interface), pointer  :: tmp_checkpoint_read_ptr

    ! Temp fields to create prognostics
    type( field_type )                         :: u, rho, theta, exner

    call log_event( 'GungHo: Creating prognostics...', LOG_LEVEL_INFO )

    ! Create prognostic fields
    theta = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, Wtheta), &
                        name= "theta" )
    if ( vorticity_in_w1 ) then
      xi    = field_type( vector_space = &
                          function_space_collection%get_fs(mesh_id, element_order, W1), &
                          name = "xi" )
    else
      xi    = field_type( vector_space = &
                          function_space_collection%get_fs(mesh_id, element_order, W2), &
                          name = "xi" )
    end if
    u     = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W2), &
                        name = "u" )
    rho   = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W3), &
                        name = "rho" )
    exner = field_type( vector_space = &
                        function_space_collection%get_fs(mesh_id, element_order, W3), &
                        name = "exner" )

    ! The moisture mixing ratio fields (mr) and moist dynamics fields
    ! (moist_dyn) are always passed into the timestep algorithm, so are
    ! always created here (even if use_moisture is false).
    do imr = 1,nummr
      mr(imr) = field_type( vector_space = &
      function_space_collection%get_fs(mesh_id, element_order, theta%which_function_space()), &
                            name = trim(mr_names(imr)) )
    end do

    ! Auxilliary fields holding moisture-dependent factors for dynamics
    do imr = 1, num_moist_factors
      moist_dyn(imr) = field_type( vector_space = &
      function_space_collection%get_fs(mesh_id, element_order, theta%which_function_space()) )
    end do

    ! Set I/O behaviours for diagnostic output

    if (write_diag .and. use_xios_io) then

       ! Set diagnostic output handlers
 
       ! Face domain

       tmp_write_diag_ptr => xios_write_field_face

       ! Vector fields that are projected to scalar components
       call xi%set_write_diag_behaviour(tmp_write_diag_ptr)
       call u%set_write_diag_behaviour(tmp_write_diag_ptr)

       ! Scalar fields
       call rho%set_write_diag_behaviour(tmp_write_diag_ptr)
       call exner%set_write_diag_behaviour(tmp_write_diag_ptr)

       ! Theta is a special case as it can be on face (if function space is WTheta)
       ! or node (if function space is W0)
       if (theta%which_function_space() == Wtheta) then

         call theta%set_write_diag_behaviour(tmp_write_diag_ptr)

       else

        tmp_write_diag_ptr => xios_write_field_node

        call theta%set_write_diag_behaviour(tmp_write_diag_ptr)

       end if

       ! Moisture uses the same type of field write as Theta

       call theta%get_write_diag_behaviour(tmp_write_diag_ptr)

       do imr = 1,nummr
         call mr(imr)%set_write_diag_behaviour(tmp_write_diag_ptr)
       end do

    end if 

    if ( checkpoint_write .or. checkpoint_read) then

      if ( use_xios_io ) then

        ! Use XIOS for checkpoint / restart

        tmp_checkpoint_write_ptr => checkpoint_write_xios
        tmp_checkpoint_read_ptr => checkpoint_read_xios

        call log_event( 'GungHo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )

      else

        ! Use old checkpoint and restart methods

        tmp_checkpoint_write_ptr => checkpoint_write_netcdf
        tmp_checkpoint_read_ptr => checkpoint_read_netcdf

       call log_event( 'GungHo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )

      end if

      call u%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call rho%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call theta%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call exner%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call u%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call rho%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call theta%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call exner%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)

      do imr = 1,nummr
        call mr(imr)%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
        call mr(imr)%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      end do

    end if

    ! Create the prognostic field collection and populate it

    prognostic_fields  =  field_collection_type(name="prognostics")
    call prognostic_fields%add_field( theta )
    call prognostic_fields%add_field( rho )
    call prognostic_fields%add_field( u )
    call prognostic_fields%add_field( exner )

    ! The moisture mixing ratios are always created, but only add them
    ! to the prognostics collection if they are active
    if ( use_moisture )then
      do imr = 1,nummr
        tmp_ptr => mr(imr)
        call prognostic_fields%add_reference_to_field( tmp_ptr )
      end do
    end if

  end subroutine create_gungho_prognostics


end module create_gungho_prognostics_mod
