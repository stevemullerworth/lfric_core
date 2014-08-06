!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes exner from the equations of state

!> @detail The kernel computes exner from the linear equation of state given a 
!>         theta and rho profile
module calc_exner_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, v0, v3, fe, cells ! the enums

use matrix_invert_mod,       only : matrix_invert
use constants_mod,           only : kappa, r_def
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_exner_kernel_type
  private
  type(arg_type) :: meta_args(7) = [  &
       arg_type(gh_write,v3,fe,.true.,.false.,.false.,.true.),        &
       arg_type(gh_read ,v3,fe,.false.,.false.,.false.,.false.),      &       
       arg_type(gh_read ,v0,fe,.true.,.false.,.false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.true.,.false.,.false.),       &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read ,v0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read ,v3,fe,.false.,.false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::calc_exner_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface calc_exner_kernel_type
   module procedure calc_exner_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public calc_exner_code
contains

type(calc_exner_kernel_type) function calc_exner_kernel_constructor() result(self)
  return
end function calc_exner_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_v3 The number of degrees of freedom per cell for v3
!! @param[in] map_v3 Integer array holding the dofmap for the cell at the base of the column for v3
!! @param[in] v3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] exner Real array the data 
!! @param[in] chi_3_v3 Real array. the physical z coordinate in v3
!! @param[in] chi_1 Real array. the physical x coordinate in v0
!! @param[in] chi_2 Real array. the physical y coordinate in v0
!! @param[in] chi_3 Real array. the physical z coordinate in v0
!! @param[in] rho Real array.   the density
!! @param[in] theta Real array. the potential temperature
!! @param[inout] gq The gaussian quadrature rule 
subroutine calc_exner_code(nlayers,ndf_v3,map_v3,v3_basis,gq,exner, &
                                                             rho,   &                                                          
                                   ndf_v0,map_v0,v0_basis,   theta, &
                                            v0_diff_basis,   chi_1, &
                                                             chi_2, &
                                                             chi_3, &                                                                                                                   
                                                             chi_3_v3 &
                                                                       )

  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers, ndf_v0, ndf_v3
  integer, intent(in) :: map_v0(ndf_v0), map_v3(ndf_v3)
  real(kind=r_def), intent(in), dimension(1,ndf_v3,ngp_h,ngp_v) :: v3_basis  
  real(kind=r_def), intent(in), dimension(1,ndf_v0,ngp_h,ngp_v) :: v0_basis 
  real(kind=r_def), intent(in), dimension(3,ndf_v0,ngp_h,ngp_v) :: v0_diff_basis 
  real(kind=r_def), intent(inout) :: exner(*)
  real(kind=r_def), intent(in)    :: rho(*), theta(*)
  real(kind=r_def), intent(in)    :: chi_1(*), chi_2(*), chi_3(*), chi_3_v3(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_v3) :: exner_e, rho_e,  exner_s, rho_s, rhs_e   
  real(kind=r_def), dimension(ndf_v0) :: theta_e,  theta_s
  real(kind=r_def), dimension(ndf_v0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_v3) :: chi_3_v3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: f
  real(kind=r_def), dimension(ndf_v3,ndf_v3) :: mass_matrix_v3, inv_mass_matrix_v3
  real(kind=r_def) :: rho_at_quad, rho_s_at_quad,                                 &
                   theta_at_quad, theta_s_at_quad,                             &
                   exner_s_at_quad
  real(kind=r_def) :: rhs_eos, basis_func
  
  do k = 0, nlayers-1
  ! Extract element arrays of rho & theta
    do df1 = 1, ndf_v3
      rho_e(df1) = rho( map_v3(df1) + k )
      chi_3_v3_e(df1) = chi_3_v3( map_v3(df1) + k )
    end do
    do df1 = 1, ndf_v0
      theta_e(df1) = theta( map_v0(df1) + k )  
      chi_1_e(df1) = chi_1( map_v0(df1) + k )
      chi_2_e(df1) = chi_2( map_v0(df1) + k )
      chi_3_e(df1) = chi_3( map_v0(df1) + k )
    end do
    call coordinate_jacobian(ndf_v0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             v0_diff_basis, jac, dj)
    call reference_profile(ndf_v0, ndf_v3, exner_s, rho_s, theta_s, chi_3_e,   &
                           chi_3_v3_e)
  ! compute the RHS integrated over one cell
    do df1 = 1, ndf_v3
      do qp2 = 1, ngp_v
        do qp1 = 1, ngp_h
          rho_at_quad   = 0.0_r_def
          rho_s_at_quad = 0.0_r_def
          exner_s_at_quad = 0.0_r_def
          do df2 = 1, ndf_v3
            basis_func = v3_basis(1,df2,qp1,qp2)
            rho_at_quad      = rho_at_quad     + rho_e(df2)  * basis_func
            rho_s_at_quad    = rho_s_at_quad   + rho_s(df2)  * basis_func
            exner_s_at_quad  = exner_s_at_quad + exner_s(df2)* basis_func
          end do
          theta_at_quad   = 0.0_r_def
          theta_s_at_quad = 0.0_r_def
          do df2 = 1, ndf_v0
            basis_func = v0_basis(1,df2,qp1,qp2)
            theta_at_quad    = theta_at_quad   + theta_e(df2) * basis_func
            theta_s_at_quad  = theta_s_at_quad + theta_s(df2) * basis_func
          end do
          rhs_eos = kappa / (1.0_r_def - kappa) * exner_s_at_quad                 &
                  *( rho_at_quad/rho_s_at_quad + theta_at_quad/theta_s_at_quad )
          f(qp1,qp2) = v3_basis(1,df1,qp1,qp2) * rhs_eos * dj(qp1,qp2)
        end do
      end do
      rhs_e(df1) =  gq%integrate(f)
    end do
  ! compute the LHS integrated over one cell and solve  
    do df1 = 1, ndf_v3
       do df2 = 1, ndf_v3
          do qp2 = 1, ngp_v
             do qp1 = 1, ngp_h
                 f(qp1,qp2) = v3_basis(1,df1,qp1,qp2) * &
                              v3_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
             end do
          end do
          mass_matrix_v3(df1,df2) = gq%integrate(f)
       end do
    end do
    call matrix_invert(mass_matrix_v3,inv_mass_matrix_v3,ndf_v3)
    exner_e(:) = matmul(inv_mass_matrix_v3(:,:),rhs_e(:))    
    do df1 = 1,ndf_v3
      exner(map_v3(df1)+k) = exner_e(df1) 
    end do    
  end do
  
end subroutine calc_exner_code

end module calc_exner_kernel_mod
