!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module mass_matrices_mod

use constants_mod, only: r_def
use gaussian_quadrature_mod, only : gaussian_quadrature_type,              &
                                    ngp_h,ngp_v

implicit none

real(kind=r_def), allocatable, dimension(:,:,:) ::                                &
             w2_mass_matrix, w1_mass_matrix, w0_mass_matrix

contains 

!> @brief subroutine to allocate global arrays for the W0, W1, W2 mass matrices
!> @param[in] ndf_w0 the number of dofs for a element in W0 function space
!> @param[in] ndf_w1 the number of dofs for a element in W1 function space
!> @param[in] ndf_w2 the number of dofs for a element in W2 function space
!> @param[in] global_number_cells the total number of cells (ncells*nlayers)
subroutine mass_matrix_init(ndf_w0,ndf_w1,ndf_w2,global_number_cells)

  integer, intent(in) :: ndf_w0,ndf_w1,ndf_w2,global_number_cells

  allocate( w0_mass_matrix(ndf_w0,ndf_w0,global_number_cells) )
  allocate( w1_mass_matrix(ndf_w1,ndf_w1,global_number_cells) ) 
  allocate( w2_mass_matrix(ndf_w2,ndf_w2,global_number_cells) )

end subroutine mass_matrix_init

!> @brief subroutine to compute global mass matrices for the W0, W1, W2 spaces
!> @param[in] cell the horizontal cell index 
!> @param[in] nlayers the number of vertical layers
!> @param[in] ndf_w0 the number of dofs for a element in W0 function space
!> @param[in] basis_w0 the basis functions for W0
!> @param[in] ndf_w1 the number of dofs for a element in W1 function space
!> @param[in] basis_w1 the basis functions for W1
!> @param[in] ndf_w2 the number of dofs for a element in W2 function space
!> @param[in] basis_w2 the basis functions for W2
!> @param[in] gq the quadrature rule
!> @param[in] diff_basis_w0 the differential basis for W0
!> @param[in] map_w0 the dofmap for W0
!> @param[in] chi1 the first coordinate array in W0 space
!> @param[in] chi2 the second coordinate array in W0 space
!> @param[in] chi2 the third coordinate array in W0 space
subroutine compute_mass_matrix(cell,nlayers,  &
                               ndf_w0,basis_w0,&
                               ndf_w1,basis_w1,&
                               ndf_w2,basis_w2,gq, &                              
                               diff_basis_w0,map_w0,chi1,chi2,chi3)

  use coordinate_jacobian_mod, only: coordinate_jacobian,                      &
                                     coordinate_jacobian_inverse

  integer,                        intent(in)    :: cell, nlayers 
  integer,                        intent(in)    :: ndf_w0, ndf_w1, ndf_w2  
  type(gaussian_quadrature_type), intent(inout) :: gq
  integer,                        intent(in)    :: map_w0(ndf_w0)
  real(kind=r_def), dimension(1,ndf_w0,ngp_h,ngp_v), intent(in)    :: basis_w0
  real(kind=r_def), dimension(3,ndf_w0,ngp_h,ngp_v), intent(in)    :: diff_basis_w0
  real(kind=r_def), dimension(3,ndf_w1,ngp_h,ngp_v), intent(in)    :: basis_w1
  real(kind=r_def), dimension(3,ndf_w2,ngp_h,ngp_v), intent(in)    :: basis_w2  
  
  real(kind=r_def), intent(in)    :: chi1(*), chi2(*), chi3(*)
  
  integer :: df1, df2, qp1, qp2, k, ik, loc
  real(kind=r_def) :: f(ngp_h,ngp_v)
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac, jac_inv
  real(kind=r_def), dimension(ndf_w0)          :: chi_1_e, chi_2_e, chi_3_e                 
     
  do k = 1,nlayers
    ik = k + (cell-1)*nlayers      
    do df1 = 1, ndf_w0
      loc = map_w0(df1) + k - 1
      chi_1_e(df1) = chi1( loc )
      chi_2_e(df1) = chi2( loc )
      chi_3_e(df1) = chi3( loc )
    end do
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             diff_basis_w0, jac, dj)
    call coordinate_jacobian_inverse(ngp_h, ngp_v, jac, dj, jac_inv)

! W0 mass matrix
    do df2 = 1, ndf_w0
      do df1 = df2, ndf_w0  
        do qp2 = 1, ngp_v
          do qp1 = 1, ngp_h                                         
            f(qp1,qp2) = basis_w0(1,df1,qp1,qp2)*basis_w0(1,df2,qp1,qp2)       &
                       * dj(qp1,qp2)
          end do
        end do
        w0_mass_matrix(df1,df2,ik) = gq%integrate(f)
      end do
      do df1 = df2, 1, -1  
        w0_mass_matrix(df1,df2,ik) = w0_mass_matrix(df2,df1,ik)
      end do
    end do 

! W1 mass matrix
    do df2 = 1, ndf_w1
      do df1 = df2, ndf_w1
        do qp2 = 1, ngp_v
          do qp1 = 1, ngp_h                                         
            f(qp1,qp2) = dot_product(                                          &
                         matmul(jac_inv(:,:,qp1,qp2),basis_w1(:,df1,qp1,qp2)), &
                         matmul(jac_inv(:,:,qp1,qp2),basis_w1(:,df2,qp1,qp2))) &
                         *dj(qp1,qp2)           
          end do
        end do
        w1_mass_matrix(df1,df2,ik) = gq%integrate(f)
      end do
      do df1 = df2, 1, -1  
        w1_mass_matrix(df1,df2,ik) = w1_mass_matrix(df2,df1,ik)
      end do
    end do  

! W2 mass matrix
    do df2 = 1, ndf_w2
      do df1 = df2, ndf_w2  
        do qp2 = 1, ngp_v
          do qp1 = 1, ngp_h                                         
            f(qp1,qp2) = dot_product(                                          &
                         matmul(jac(:,:,qp1,qp2),basis_w2(:,df1,qp1,qp2)),     &
                         matmul(jac(:,:,qp1,qp2),basis_w2(:,df2,qp1,qp2)))     &
                         /dj(qp1,qp2)           
          end do
        end do
        w2_mass_matrix(df1,df2,ik) = gq%integrate(f)
      end do
      do df1 = df2, 1, -1  
        w2_mass_matrix(df1,df2,ik) = w2_mass_matrix(df2,df1,ik)
      end do
    end do 
              
  end do

end subroutine compute_mass_matrix

end module mass_matrices_mod
