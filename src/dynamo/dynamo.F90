!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using 
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls Creates function spaces, then fields on those 
!> function spaces, before passing the fields to the algorithm layer

program dynamo

  use dynamo_algorithm_rk_timestep_mod, &
                               only : dynamo_algorithm_rk_timestep
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type, V0, V1, V2, V3
  use log_mod,                 only : log_event, LOG_LEVEL_INFO
  use set_up_mod,              only : set_up
  use gaussian_quadrature_mod, only : gaussian_quadrature_type

  implicit none

  type( function_space_type )      :: function_space
! coordinate fields
  type( field_type ) :: chi(3)
! prognostic fields    
  type( field_type ) :: u, rho, theta, exner, xi
                                     
  type( gaussian_quadrature_type ) :: gq
  integer                          :: coord

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  call set_up( )

  do coord = 1,3
    chi(coord) = field_type( vector_space = function_space%get_instance( V0 ),&
                             gq = gq%get_instance() )
  end do
               
  theta = field_type( vector_space = function_space%get_instance( V0 ),       &
                      gq = gq%get_instance() )
                    
  xi = field_type( vector_space = function_space%get_instance( V1 ),          &
                      gq = gq%get_instance() )
                    
  u = field_type( vector_space = function_space%get_instance( V2 ),           &
                      gq = gq%get_instance() )

  rho = field_type( vector_space = function_space%get_instance( V3 ),         &
                      gq = gq%get_instance() )

  exner = field_type( vector_space = function_space%get_instance( V3 ),       &
                      gq = gq%get_instance() )
                                           

  call dynamo_algorithm_rk_timestep( chi, u, rho, theta, exner, xi)                       

! do some i/o
  call rho%print_field( '   rho =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call exner%print_field( '   exner =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call theta%print_field( '   theta =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call u%print_field( '   u =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

end program dynamo
