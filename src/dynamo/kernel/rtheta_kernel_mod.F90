!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the thermodynamic equation

!> @detail The kernel computes the rhs of the thermodynamic equation for the linear equations with
!>         no advection
module rtheta_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_inc, v0, v2, v3, fe, cells ! the enums
use reference_profile_mod,   only : reference_profile
use constants_mod,           only : n_sq, gravity

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rtheta_kernel_type
  private
  type(arg_type) :: meta_args(5) = [  &
       arg_type(gh_inc  ,v0,fe,.true., .false.,.false.,.true.),        &
       arg_type(gh_read ,v2,fe,.true., .false.,.false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.true., .false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,v3,fe,.false., .false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::rtheta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rtheta_kernel_type
   module procedure rtheta_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rtheta_code
contains

type(rtheta_kernel_type) function rtheta_kernel_constructor() result(self)
  return
end function rtheta_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_v0 The number of degrees of freedom per cell for v3
!! @param[in] map_v0 Integer array holding the dofmap for the cell at the base of the column for v3
!! @param[in] v0_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] v0_diff_basis Real 5-dim array holding differential of the basis functions evaluated at gaussian quadrature points 
!! @param[inout] r_theta Real array the data 
!! @param[in] chi_1 Real array. the physical x coordinate in v0
!! @param[in] chi_2 Real array. the physical x coordinate in v0
!! @param[in] chi_3 Real array. the physical x coordinate in v0
!! @param[in] u Real array. the velocity
!! @param[inout] gq The gaussian quadrature rule 
subroutine rtheta_code(nlayers,ndf_v0, map_v0, v0_basis, gq, r_theta,          &
                               ndf_v2, map_v2, v2_basis, orientation, u,       &
                               v0_diff_basis, chi_1, chi_2, chi_3,             &
                               ndf_v3, map_v3, chi_v3_3                        &
                               )
                               
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile                               
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_v0, ndf_v2, ndf_v3
  integer, intent(in) :: map_v0(ndf_v0), map_v2(ndf_v2),  map_v3(ndf_v3)
  integer, intent(in), dimension(ndf_v2) :: orientation
  real(kind=r_def), intent(in), dimension(1,ndf_v0,ngp_h,ngp_v) :: v0_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_v0,ngp_h,ngp_v) :: v0_diff_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_v2,ngp_h,ngp_v) :: v2_basis 
  real(kind=r_def), intent(inout) :: r_theta(*)
  real(kind=r_def), intent(in) :: chi_1(*), chi_2(*), chi_3(*), u(*), chi_v3_3(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def) :: exner_s(ndf_v3), rho_s(ndf_v3), theta_s(ndf_v0)
  real(kind=r_def), dimension(ndf_v0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_v3) :: chi_v3_3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ngp_h,ngp_v,ndf_v0) :: f
  real(kind=r_def) :: u_at_quad(3), k_vec(3), vec_term, u_e(ndf_v2)
  real(kind=r_def) :: theta_s_at_quad
  real(kind=r_def) :: bouy_term
  
  k_vec(:) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  
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
    do df = 1, ndf_v2
      u_e(df) = u( map_v2(df) + k )*real(orientation(df))
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, ngp_v
      do qp1 = 1, ngp_h
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_v2
          u_at_quad(:)  = u_at_quad(:)  + u_e(df)*v2_basis(:,df,qp1,qp2)
        end do
        theta_s_at_quad = 0.0_r_def
        do df = 1, ndf_v0
          theta_s_at_quad = theta_s_at_quad                                    &
                          + theta_s(df)*v0_basis(1,df,qp1,qp2)
        end do

        vec_term = dot_product(k_vec,matmul(jac(:,:,qp1,qp2),u_at_quad))
        bouy_term =  - n_sq/gravity*theta_s_at_quad*vec_term
        
        do df = 1, ndf_v0
          f(qp1,qp2,df) = v0_basis(1,df,qp1,qp2)*bouy_term
        end do
      end do
    end do
    do df = 1, ndf_v0
      r_theta( map_v0(df) + k ) =  r_theta( map_v0(df) + k )                   &
                                + gq%integrate(f(:,:,df))
    end do 
  end do
  
end subroutine rtheta_code

end module rtheta_kernel_mod
