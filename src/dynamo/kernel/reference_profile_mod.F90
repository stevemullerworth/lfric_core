module reference_profile_mod
use constants_mod, only: r_def, n_sq, gravity, cp, rd, kappa, p_zero

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the reference profile for a single element
!! @param[in] ndf_v0     Integer. The size of the v0 arrays
!! @param[in] ndf_v3     Integer. The size of the v3 arrays
!! @param[in] exner_s    Real 1-dim array. Holds the exner reference profile
!! @param[in] rho_s      Real 1-dim array. Holds the rho reference profile
!! @param[in] theta_s    Real 1-dim array. Holds the theta reference profile
!! @param[in] z_v0       Real 1-dim array. Holds the z coordinate field for v0
!! @param[in] z_v3       Real 1-dim array. Holds the z coordinate field for v3
subroutine reference_profile(ndf_v0,ndf_v3,exner_s,rho_s,theta_s,z_v0,z_v3)

integer,       intent(in)     :: ndf_v0, ndf_v3
real(kind=r_def), intent(in)  :: z_v0(ndf_v0), z_v3(ndf_v3)
real(kind=r_def), intent(out) :: exner_s(ndf_v3), rho_s(ndf_v3), theta_s(ndf_v0)

real(kind=r_def), parameter :: theta_surf = 300.0_r_def
real(kind=r_def), parameter :: exner_surf = 1.0_r_def
real(kind=r_def), parameter :: rho_surf   = 1.0_r_def
real(kind=r_def)            :: theta_v3, nsq_over_g

integer :: df

nsq_over_g = n_sq/gravity

do df = 1, ndf_v0
  theta_s(df) = theta_surf * exp ( nsq_over_g * z_v0(df) )
end do
do df = 1, ndf_v3
  exner_s(df) = exner_surf - gravity**2/(cp*theta_surf*n_sq)   &
              * (1.0_r_def - exp ( - nsq_over_g * z_v3(df) ))
  
  theta_v3    = theta_surf * exp ( nsq_over_g * z_v3(df) )
  rho_s(df)   = p_zero/(rd*theta_v3) * exner_s(df) ** ((1.0_r_def - kappa)/kappa) 
end do

end subroutine reference_profile

end module reference_profile_mod
