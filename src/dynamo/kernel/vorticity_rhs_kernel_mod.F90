!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the absolute vorticity

!> @detail Compute the projection of the absolute vorticity into the w1 space. 
!>         The is: curl(u) + 2*Omega -> curl(c).u + 2*c.Omega
module vorticity_rhs_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    W0, W1, W2,                              &
                                    GH_BASIS, GH_DIFF_BASIS, GH_ORIENTATION, &
                                    CELLS
use constants_mod,           only : r_def
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vorticity_rhs_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W1),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, W0)                               &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W1, GH_DIFF_BASIS, GH_ORIENTATION),                   &
       func_type(W2, GH_BASIS,      GH_ORIENTATION),                   &
       func_type(W0, GH_DIFF_BASIS)                                    &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::vorticity_rhs_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface vorticity_rhs_kernel_type
   module procedure vorticity_rhs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public vorticity_rhs_code
contains

type(vorticity_rhs_kernel_type) function vorticity_rhs_kernel_constructor() result(self)
  return
end function vorticity_rhs_kernel_constructor

!> @brief 
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_xi The number of degrees of freedom per cell for w1
!! @param[in] undf_xi The unique number of degrees of freedom for w1
!! @param[in] map_xi Integer array holding the dofmap for the cell at the base of the column for w1
!! @param[in] orientation_xi Local orientation array for the vorticity field
!! @param[inout] rhs Real array the field to contain the right hand side to be computed
!! @param[in] ndf_u The number of degrees of freedom per cell for the velocity field
!! @param[in] undf_u The unique number of degrees of freedom for the velocity field
!! @param[in] map_u Integer array holding the dofmap for the cell at the base of the column for the velocity field
!! @param[in] basis_u Real 4-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] orientation_u Local orientation array for the velocity field
!! @param[in] u Real array the velocity field
!! @param[in] ndf_chi The number of degrees of freedom per cell for the function space containing chi
!! @param[in] undf_chi The unique number of degrees of freedom for the chi arrays
!! @param[in] map_chi Integer array holding the dofmap for the cell at the base of the column for the function space containing chi
!! @param[in] chi_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] nqp_h the number of horizontal quadrature points
!! @param[in] nqp_v the number of vertical quadrature points
!! @param[in] wqp_h the weights of the horizontal quadrature points
!! @param[in] wqp_v the weights of the vertical quadrature points
subroutine vorticity_rhs_code(nlayers,                                          &
                         rhs, u, chi_1, chi_2, chi_3,                           &
                         ndf_xi, undf_xi, map_xi, diff_basis_xi, orientation_xi,&
                         ndf_u, undf_u, map_u, basis_u, orientation_u,          &
                         ndf_chi, undf_chi, map_chi, diff_basis_chi,            &                           
                         nqp_h, nqp_v, wqp_h, wqp_v                             &
                         )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_chi, ndf_u, ndf_xi, undf_chi, undf_u, undf_xi  
  integer, intent(in) :: nqp_h, nqp_v
  integer, dimension(ndf_xi),  intent(in) :: map_xi
  integer, dimension(ndf_u),   intent(in) :: map_u
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_xi),  intent(in) :: orientation_xi
  integer, dimension(ndf_u),   intent(in) :: orientation_u
  real(kind=r_def), dimension(3,ndf_u,  nqp_h,nqp_v), intent(in) :: basis_u  
  real(kind=r_def), dimension(3,ndf_xi, nqp_h,nqp_v), intent(in) :: diff_basis_xi
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(undf_xi),               intent(inout) :: rhs
  real(kind=r_def), dimension(undf_u),                intent(in)    :: u 
  real(kind=r_def), dimension(undf_chi),              intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)        :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)    :: jac
  real(kind=r_def), dimension(ndf_u)  :: u_cell(ndf_u) 
  real(kind=r_def), dimension(ndf_xi) :: rhs_cell(ndf_xi)
  real(kind=r_def) :: u_at_quad(3), &
                      dc(3)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             diff_basis_chi, jac, dj)    
    do df = 1, ndf_u
      u_cell(df) = u( map_u(df) + k )*real(orientation_u(df),r_def)
    end do    
  ! compute the RHS integrated over one cell
    rhs_cell(:) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_u
          u_at_quad(:) = u_at_quad(:) + u_cell(df)*basis_u(:,df,qp1,qp2)
        end do        
        u_at_quad(:) = wqp_h(qp1)*wqp_v(qp2)*matmul(jac(:,:,qp1,qp2),u_at_quad(:))/dj(qp1,qp2)
        do df = 1, ndf_xi
          dc = matmul(jac(:,:,qp1,qp2), &
                      diff_basis_xi(:,df,qp1,qp2)*real(orientation_xi(df),r_def))
          rhs_cell(df) = rhs_cell(df) -  dot_product(dc,u_at_quad(:))
        end do
      end do
    end do
    do df = 1, ndf_xi
      rhs( map_xi(df) + k ) =  rhs( map_xi(df) + k ) + rhs_cell(df)
    end do 
  end do 
end subroutine vorticity_rhs_code

end module vorticity_rhs_kernel_mod
