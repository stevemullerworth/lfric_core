!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Subroutine to enforces the  boundary condition

module enforce_bc_mod

use constants_mod, only: r_def

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

!> @brief Subroutine to enforce zero flux bc in W2
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[in] boundary_value integer array holding flag =0 for dofs on top and bottom faces of an cell
!! @param[inout] u Real array, the actual data
subroutine enforce_bc_w2(nlayers,ndf,map,boundary_value,u)
  
  !Arguments
  integer, intent(in) :: nlayers, ndf
  integer, intent(in) :: map(ndf)
  integer, intent(in) :: boundary_value(ndf,2)
  real(kind=r_def), intent(inout) :: u(*)

  !Internal variables
  integer               :: df, k
  
  do df = 1,ndf
    k = 0
    u(map(df) + k) = u(map(df) + k)*real(boundary_value(df,1))
    k = nlayers - 1  
    u(map(df) + k) = u(map(df) + k)*real(boundary_value(df,2))
  end do
  
end subroutine enforce_bc_w2

end module enforce_bc_mod
