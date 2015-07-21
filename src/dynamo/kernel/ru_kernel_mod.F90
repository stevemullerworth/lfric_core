!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the momentum equation for the nonlinear equations,
!>         written in the vector invariant form


!> @detail The kernel computes the rhs of the momentum equation for the nonlinear equations,
!>         written in the vector invariant form
!>         This consists of four terms:
!>         Pressure gradient: cp*theta*grad(pi)
!>         geopotential gradient: grad(Phi) ( \quiv g for some domains)
!>         gradient of kinetic energy: grad(1/2*u.u) = 
!>         vorticity advection: xi/rho \cross F (with vorticity xi and mass flux F)
!>         This results in:
!>         ru = -xi/rho x F - grad(Phi + 1/2*u.u) - cp*theta*grad(exner)
module ru_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    W0, W1, W2, W3, GH_BASIS, GH_DIFF_BASIS, &
                                    GH_ORIENTATION, CELLS
use constants_mod,           only : cp, r_def
use cross_product_mod,       only : cross_product

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ru_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD,   GH_READ, W1),                              &
       arg_type(GH_FIELD*3, GH_READ, W0)                               &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS, GH_ORIENTATION),         &
       func_type(W3, GH_BASIS),                                        &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W1, GH_BASIS, GH_ORIENTATION)                         &
       /)
  integer :: iterates_over = CELLS
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
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number unique of degrees of freedom  for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Real 4-dim array holding basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at  quadrature points
!! @param[in] boundary_value array of flags (= 0) for dofs that live on the
!!            vertical boundaries of the cell (=1 for other dofs)
!! @param[in] orientation_w2 Orientation array for w2 fields
!! @param[inout] r_u Real array the data 
!! @param[in] u The velocity array
!! @param[in] mass_flux The mass flux array F
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number unique of degrees of freedom  for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] rho Real array. The density
!! @param[in] ndf_w1 The number of degrees of freedom per cell for w1
!! @param[in] undf_w1 The number unique of degrees of freedom  for w1
!! @param[in] map_w1 Integer array holding the dofmap for the cell at the base of the column for w1
!! @param[in] w1_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] orientation_w1 Orientation array for w1 fields
!! @param[in] xi Real array. Vorticity
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] undf_w0 The number unique of degrees of freedom  for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] theta Real array. potential temperature
!! @param[in] phi The geopotential
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] nqp_h Integer, number of quadrature points in the horizontal
!! @param[in] nqp_v Integer, number of quadrature points in the vertical
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical

subroutine ru_code(nlayers,                                                    &
                   r_u, u, mass_flux, rho, theta, phi, xi,                     &
                   chi_1, chi_2, chi_3,                                        &
                   ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis,           &
                   orientation_w2, boundary_value_w2,                          &
                   ndf_w3, undf_w3, map_w3, w3_basis,                          &
                   ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis,           &
                   ndf_w1, undf_w1, map_w1, w1_basis, orientation_w1,          &
                   nqp_h, nqp_v, wqp_h, wqp_v                                  &
                   )

  use coordinate_jacobian_mod,  only: coordinate_jacobian, &
                                      coordinate_jacobian_inverse
  use enforce_bc_mod,           only: enforce_bc_w2
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  
  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w1, ndf_w2, ndf_w3
  integer, intent(in) :: undf_w0, undf_w1, undf_w2, undf_w3
  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_w1), intent(in) :: map_w1
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w3), intent(in) :: map_w3
  
  integer, dimension(ndf_w2,2), intent(in) :: boundary_value_w2
  integer, dimension(ndf_w2),   intent(in) :: orientation_w2
  integer, dimension(ndf_w1),   intent(in) :: orientation_w1


  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis  
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(3,ndf_w1,nqp_h,nqp_v), intent(in) :: w1_basis 
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_basis 
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_diff_basis   

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u, mass_flux
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_w1), intent(in)    :: xi
  real(kind=r_def), dimension(undf_w0), intent(in)    :: chi_1, chi_2, chi_3, theta, phi

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac, jac_inv
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_w2)          :: ru_e, u_e, f_e
  real(kind=r_def), dimension(ndf_w1)          :: xi_e
  real(kind=r_def), dimension(ndf_w0)          :: theta_e, phi_e

  real(kind=r_def) :: grad_theta_at_quad(3), &
                      u_at_quad(3), &
                      xi_at_quad(3), &
                      f_at_quad(3),     &
                      jac_v(3), &
                      vorticity(3), &
                      grad_phi_at_quad(3), &
                      v(3)
  real(kind=r_def) :: exner_at_quad, rho_at_quad, theta_at_quad, &
                      ke_at_quad,                                &
                      grad_term, geo_term, ke_term,              &
                      vorticity_term,                            &
                      dv
  
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      loc = map_w0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_w0, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    call coordinate_jacobian_inverse(nqp_h, nqp_v, jac, dj, jac_inv)   

    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do    
    do df = 1, ndf_w0
      theta_e(df) = theta( map_w0(df) + k )
      phi_e(df)   = phi(map_w0(df) + k)
    end do   
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
      f_e(df) = mass_flux( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do    
    do df = 1, ndf_w1
      xi_e(df) = xi( map_w1(df) + k )
    end do    
  ! compute the RHS integrated over one cell    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        grad_theta_at_quad(:) = 0.0_r_def
        grad_phi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*w0_basis(1,df,qp1,qp2)
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*w0_diff_basis(:,df,qp1,qp2) 
          grad_phi_at_quad(:)   = grad_phi_at_quad(:) &
                                + phi_e(df)*w0_diff_basis(:,df,qp1,qp2) 

        end do

        exner_at_quad = calc_exner_pointwise(rho_at_quad, theta_at_quad)

! k.e term
        u_at_quad(:) = 0.0_r_def
        f_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2) &
                         *real(orientation_w2(df),r_def)
          f_at_quad(:) = f_at_quad(:) &
                       + f_e(df)*w2_basis(:,df,qp1,qp2) &
                         *real(orientation_w2(df),r_def)
        end do
        ke_at_quad = 0.5_r_def*dot_product(matmul(jac(:,:,qp1,qp2),u_at_quad), &
                                           matmul(jac(:,:,qp1,qp2),u_at_quad))/(dj(qp1,qp2)**2)

! Vorticity advection term
        xi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w1
          xi_at_quad(:) = xi_at_quad(:) &
                        + xi_e(df)*w1_basis(:,df,qp1,qp2) &
                          *real(orientation_w1(df),r_def)
        end do
        vorticity = cross_product(matmul(transpose(jac_inv(:,:,qp1,qp2)),xi_at_quad),&
                                  matmul(jac    (:,:,qp1,qp2),f_at_quad))
        vorticity = vorticity/(dj(qp1,qp2)*rho_at_quad)

        do df = 1, ndf_w2
          v  = w2_basis(:,df,qp1,qp2)*real(orientation_w2(df),r_def)
          dv = w2_diff_basis(1,df,qp1,qp2)*real(orientation_w2(df),r_def)

          jac_v = matmul(jac(:,:,qp1,qp2),v)
! geopotential term
          geo_term = dot_product( grad_phi_at_quad(:), v)
!pressure gradient term
          grad_term = cp*exner_at_quad * (                           & 
                      theta_at_quad * dv                             &
                    + dot_product( grad_theta_at_quad(:),v)          &
                                         )

          ke_term =  dv*ke_at_quad

          vorticity_term = - dot_product(jac_v,vorticity)

          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*( grad_term &
                                                       - geo_term &
                                                       + ke_term &
                                                       + vorticity_term)

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do

  call enforce_bc_w2(nlayers,ndf_w2,undf_w2,map_w2,boundary_value_w2,r_u)
end subroutine ru_code

end module ru_kernel_mod
