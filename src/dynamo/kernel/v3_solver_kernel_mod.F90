!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes LHS of Galerkin projection and solves equation in W3 space

module v3_solver_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use gaussian_quadrature_mod, only : ngp_h, ngp_v, gaussian_quadrature_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, v0, v3, fe, cells ! the enums

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: v3_solver_kernel_type
  private
  type(arg_type) :: meta_args(5) = [  &
       arg_type(gh_write,v3,fe,.true.,.false.,.false.,.true.),        &
       arg_type(gh_read ,v3,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read, v0,fe,.false.,.true.,.false., .false.),      &
       arg_type(gh_read, v0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read, v0,fe,.false.,.false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::solver_v3_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface v3_solver_kernel_type
   module procedure v3_solver_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public solver_v3_code
contains

type(v3_solver_kernel_type) function v3_solver_kernel_constructor() result(self)
  return
end function v3_solver_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_v3 The number of degrees of freedom per cell
!! @param[in] map_v3 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] v3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] X Real array the data 
!! @param[in] rhs Real array. the data
!! @param[inout] gq The gaussian quadrature rule 
!! @param[in] ndf_v0 The number of degrees of freedom per cell
!! @param[in] map_v0 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] v0_diff_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[inout] chi_1 Real array, the x component of the v0 coordinate field
!! @param[inout] chi_2 Real array, the y component of the v0 coordinate field
!! @param[inout] chi_3 Real array, the z component of the v0 coordinate field
subroutine solver_v3_code(nlayers, ndf_v3, map_v3, v3_basis, x, rhs, gq, &
                          ndf_v0, map_v0, v0_diff_basis, chi_1, chi_2, chi_3 &
                         )
                         
   use matrix_invert_mod,       only : matrix_invert 
   use coordinate_jacobian_mod, only : coordinate_jacobian 

  ! needs to compute the integral of rho_df * P 
  ! P_analytic over a single column    
  
  !Arguments
  integer, intent(in) :: nlayers, ndf_v3, ndf_v0
  integer, intent(in) :: map_v3(ndf_v3), map_v0(ndf_v0)
  real(kind=r_def), intent(in), dimension(1,ndf_v3,ngp_h,ngp_v) :: v3_basis  
  real(kind=r_def), intent(inout) :: x(*)
  real(kind=r_def), intent(in) :: rhs(*)
  real(kind=r_def), intent(in)    :: chi_1(*), chi_2(*), chi_3(*)
  type(gaussian_quadrature_type), intent(in)                 :: gq
  real(kind=r_def), intent(in), dimension(3,ndf_v0,ngp_h,ngp_v) :: v0_diff_basis 

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2
  
  real(kind=r_def) :: x_e(ndf_v3), rhs_e(ndf_v3)
  real(kind=r_def), dimension(ngp_h,ngp_v) :: f
  real(kind=r_def), dimension(ndf_v3,ndf_v3) :: mass_matrix_v3, inv_mass_matrix_v3
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ndf_v0) :: chi_1_e, chi_2_e, chi_3_e
  
  ! compute the LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_v0
      chi_1_e(df1) = chi_1( map_v0(df1) + k)
      chi_2_e(df1) = chi_2( map_v0(df1) + k)
      chi_3_e(df1) = chi_3( map_v0(df1) + k)
    end do
    call coordinate_jacobian(ndf_v0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e, v0_diff_basis, jac, dj)
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
       rhs_e(df1) = rhs(map_v3(df1)+k)
    end do
    call matrix_invert(mass_matrix_v3,inv_mass_matrix_v3,ndf_v3)
    x_e = matmul(inv_mass_matrix_v3,rhs_e)
    do df1 = 1,ndf_v3
      x(map_v3(df1)+k) = x_e(df1) 
    end do
  end do
  
end subroutine solver_v3_code

end module v3_solver_kernel_mod
