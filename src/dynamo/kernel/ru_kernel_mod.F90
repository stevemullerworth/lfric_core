!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the momentum equation

!> @detail The kernel computes thr rhs of the momentum equation for the linear equations with
!>         no advection
module ru_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_inc, v3, v2, v0, fe, cells ! the enums
use constants_mod,           only : n_sq, gravity, cp, r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ru_kernel_type
  private
  type(arg_type) :: meta_args(7) = [  &
       arg_type(gh_inc  ,v2,fe,.true., .true.,.false.,.true.),         &
       arg_type(gh_read ,v3,fe,.true.,.false.,.false.,.false.),        &
       arg_type(gh_read ,v0,fe,.false.,.true.,.false., .false.),       &
       arg_type(gh_read ,v0,fe,.false.,.false.,.true.,.false.),        &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,v3,fe,.false.,.false.,.false.,.false.)        &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::ru_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface ru_kernel_type
   module procedure ru_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public ru_code
contains

type(ru_kernel_type) function ru_kernel_constructor() result(self)
  return
end function ru_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_v2 The number of degrees of freedom per cell for v2
!! @param[in] map_v2 Integer array holding the dofmap for the cell at the base of the column for v2
!! @param[in] v2_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] v2_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature points
!! @param[inout] r_u Real array the data 
!! @param[in] ndf_v0 The number of degrees of freedom per cell for v0
!! @param[in] map_v0 Integer array holding the dofmap for the cell at the base of the column for v0
!! @param[in] v0_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] v0_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Real array. the physical x coordinate in v0
!! @param[in] chi_2 Real array. the physical y coordinate in v0
!! @param[in] chi_3 Real array. the physical z coordinate in v0
!! @param[in] theta Real array. potential temperature
!! @param[in] ndf_v3 The number of degrees of freedom per cell for v3
!! @param[in] map_v3 Integer array holding the dofmap for the cell at the base of the column for v3
!! @param[in] v3_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] exner Real array. exner pressure
!! @param[inout] gq The gaussian quadrature rule 
!! @param[in] chi_3_v3 Real array. the physical z coordinate in v3
subroutine ru_code(nlayers,ndf_v2, map_v2, v2_basis, v2_diff_basis, gq,        &
                           boundary_value, r_u,                                &
                           ndf_v3, map_v3, v3_basis, exner,                    &
                           ndf_v0, map_v0, v0_basis, theta,                    &                           
                           v0_diff_basis, chi_1, chi_2, chi_3,                 &
                           chi_v3_3                                            &
                           )
                           
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile 
  use enforce_bc_mod,          only: enforce_bc_w2
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_v0, ndf_v2, ndf_v3
  integer, intent(in) :: map_v0(ndf_v0), map_v2(ndf_v2), map_v3(ndf_v3)
  integer, intent(in), dimension(ndf_v2,2) :: boundary_value
  real(kind=r_def), intent(in), dimension(1,ndf_v3,ngp_h,ngp_v) :: v3_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_v2,ngp_h,ngp_v) :: v2_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_v0,ngp_h,ngp_v) :: v0_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_v2,ngp_h,ngp_v) :: v2_diff_basis
  real(kind=r_def), intent(in), dimension(3,ndf_v0,ngp_h,ngp_v) :: v0_diff_basis   
  real(kind=r_def), intent(inout) :: r_u(*)
  real(kind=r_def), intent(in) ::  exner(*), theta(*)  
  real(kind=r_def), intent(in) :: chi_1(*), chi_2(*), chi_3(*), chi_v3_3(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def) :: exner_s(ndf_v3), rho_s(ndf_v3), theta_s(ndf_v0)
  real(kind=r_def), dimension(ndf_v0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_v3) :: chi_v3_3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)        :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v)    :: jac
  real(kind=r_def), dimension(ngp_h,ngp_v,ndf_v2) :: f
  real(kind=r_def) :: exner_e(ndf_v3), theta_e(ndf_v0)
  real(kind=r_def) :: exner_at_quad, theta_at_quad, theta_s_at_quad,              &
                   grad_term, bouy_term
  real(kind=r_def) ::k_vec(3), grad_theta_s_at_quad(3) 
  
  k_vec = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_v0
      loc = map_v0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
     do df = 1, ndf_v3
      chi_v3_3_e(df) = chi_v3_3( map_v3(df) + k )
    end do
    call coordinate_jacobian(ndf_v0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             v0_diff_basis, jac, dj)
    call reference_profile(ndf_v0, ndf_v3, exner_s, rho_s, theta_s, chi_3_e,   &
                           chi_v3_3_e)
    do df = 1, ndf_v3
      exner_e(df) = exner( map_v3(df) + k )
    end do    
    do df = 1, ndf_v0
      theta_e(df) = theta( map_v0(df) + k )
    end do    
  ! compute the RHS integrated over one cell
    do qp2 = 1, ngp_v
      do qp1 = 1, ngp_h
        exner_at_quad = 0.0_r_def
        do df = 1, ndf_v3
          exner_at_quad   = exner_at_quad + exner_e(df)*v3_basis(1,df,qp1,qp2)
        end do 
        
        theta_at_quad           = 0.0_r_def
        theta_s_at_quad         = 0.0_r_def
        grad_theta_s_at_quad(:) = 0.0_r_def
        do df = 1, ndf_v0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*v0_basis(1,df,qp1,qp2)
          theta_s_at_quad = theta_s_at_quad                                    &
                          + theta_s(df)*v0_basis(1,df,qp1,qp2)
          grad_theta_s_at_quad(:) = grad_theta_s_at_quad(:)                    &
                                  + theta_s(df)*v0_diff_basis(:,df,qp1,qp2) 
        end do
               
        do df = 1, ndf_v2
          bouy_term = dot_product(                                             &
                      matmul(jac(:,:,qp1,qp2),v2_basis(:,df,qp1,qp2)),         &
                      theta_at_quad/theta_s_at_quad*gravity*k_vec(:))
          grad_term = cp * (theta_s_at_quad * v2_diff_basis(1,df,qp1,qp2)      &
                           + dot_product(v2_basis(:,df,qp1,qp2),               &
                                         grad_theta_s_at_quad(:))              &
                           )*exner_at_quad
        
          f(qp1,qp2,df) = ( grad_term + bouy_term )
        end do
      end do
    end do
    do df = 1, ndf_v2
      r_u( map_v2(df) + k ) =  r_u( map_v2(df) + k ) + gq%integrate(f(:,:,df))
    end do 
  end do 
  
  call enforce_bc_w2(nlayers,ndf_v2,map_v2,boundary_value,r_u)
end subroutine ru_code

end module ru_kernel_mod
