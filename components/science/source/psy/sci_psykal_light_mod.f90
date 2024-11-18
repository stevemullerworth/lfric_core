!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module sci_psykal_light_mod

  use constants_mod,      only : i_def, r_def, r_double, r_solver
  use field_mod,          only : field_type, field_proxy_type
  use r_solver_field_mod, only : r_solver_field_type, r_solver_field_proxy_type

  implicit none

  public

contains

  !-------------------------------------------------------------------------------
  !> This PSyKAl-lite code is required because, currently, PSYclone does not support
  !> the output of scalar variables from kernels.(See PSyclone issue #1818)
  !> This subroutine recovers a scalar value from a field. This is required
  !> as scalars can't currently be written to checkpoint files. The workaround is to
  !> copy the scalar to a field, which may then be checkpointed. On a restart the
  !> scalar value needs to be recovered from the checkpointed field.
    subroutine invoke_getvalue(field, val)
      implicit none
      real(r_def), intent(out)     :: val
      type(field_type), intent(in) :: field
      type(field_proxy_type)       :: field_proxy
      field_proxy = field%get_proxy()
      val = field_proxy%data(1)
    end subroutine invoke_getvalue


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Perform innerproduct of a r_solver precision field in r_double precision
    subroutine invoke_rdouble_X_innerproduct_X(field_norm, field)

      use scalar_mod,         only: scalar_type
      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type
      use r_solver_field_mod, only: r_solver_field_type, r_solver_field_proxy_type

      implicit none

      real(kind=r_def), intent(out) :: field_norm
      type(r_solver_field_type), intent(in) :: field

      type(scalar_type)                           :: global_sum
      integer(kind=i_def)                         :: df
      real(kind=r_double), allocatable, dimension(:) :: l_field_norm
      integer(kind=i_def)                         :: th_idx
      integer(kind=i_def)                         :: loop0_start, loop0_stop
      integer(kind=i_def)                         :: nthreads
      type(r_solver_field_proxy_type)             :: field_proxy
      integer(kind=i_def)                         :: max_halo_depth_mesh
      type(mesh_type), pointer                    :: mesh => null()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise field and/or operator proxies
      !
      field_proxy = field%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => field_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = field_proxy%vspace%get_last_dof_owned()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      field_norm = 0.0_r_def
      ALLOCATE (l_field_norm(nthreads))
      l_field_norm = 0.0_r_double
      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        l_field_norm(th_idx) = l_field_norm(th_idx) + real(field_proxy%data(df),r_double)**2
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        field_norm = field_norm+real(l_field_norm(th_idx),r_def)
      END DO
      DEALLOCATE (l_field_norm)
      global_sum%value = field_norm
      field_norm = global_sum%get_sum()
      !
    end subroutine invoke_rdouble_X_innerproduct_X


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Psyclone does not currently have native support for builtins with mixed
    ! precision, this will be addressed in https://github.com/stfc/PSyclone/issues/1786
    ! Perform innerproduct of a r_solver precision field in r_def precision
    subroutine invoke_rdouble_X_innerproduct_Y(field_norm, field1, field2)

      use scalar_mod,         only: scalar_type
      use omp_lib,            only: omp_get_thread_num
      use omp_lib,            only: omp_get_max_threads
      use mesh_mod,           only: mesh_type
      use r_solver_field_mod, only: r_solver_field_type, r_solver_field_proxy_type

      implicit none

      real(kind=r_def), intent(out) :: field_norm
      type(r_solver_field_type), intent(in) :: field1, field2

      type(scalar_type)                           :: global_sum
      integer(kind=i_def)                         :: df
      real(kind=r_double), allocatable, dimension(:) :: l_field_norm
      integer(kind=i_def)                         :: th_idx
      integer(kind=i_def)                         :: loop0_start, loop0_stop
      integer(kind=i_def)                         :: nthreads
      type(r_solver_field_proxy_type)             :: field1_proxy, field2_proxy
      integer(kind=i_def)                         :: max_halo_depth_mesh
      type(mesh_type), pointer                    :: mesh => null()
      !
      ! Determine the number of OpenMP threads
      !
      nthreads = omp_get_max_threads()
      !
      ! Initialise field and/or operator proxies
      !
      field1_proxy = field1%get_proxy()
      field2_proxy = field2%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => field1_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = field1_proxy%vspace%get_last_dof_owned()
      !
      ! Call kernels and communication routines
      !
      !
      ! Zero summation variables
      !
      field_norm = 0.0_r_def
      ALLOCATE (l_field_norm(nthreads))
      l_field_norm = 0.0_r_double
      !
      !$omp parallel default(shared), private(df,th_idx)
      th_idx = omp_get_thread_num()+1
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        l_field_norm(th_idx) = l_field_norm(th_idx) + real(field1_proxy%data(df),r_double)*real(field2_proxy%data(df),r_double)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! sum the partial results sequentially
      !
      DO th_idx=1,nthreads
        field_norm = field_norm+real(l_field_norm(th_idx),r_def)
      END DO
      DEALLOCATE (l_field_norm)
      global_sum%value = field_norm
      field_norm = global_sum%get_sum()
      !
    end subroutine invoke_rdouble_X_innerproduct_Y


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine invoke_inc_rdefX_plus_rsolverY(X, Y)

      use mesh_mod, only: mesh_type

      implicit none

      type(field_type),          intent(inout) :: X
      type(r_solver_field_type), intent(in)    :: Y
      integer(kind=i_def) :: df
      integer(kind=i_def) :: loop0_start, loop0_stop
      type(field_proxy_type) :: X_proxy
      type(r_solver_field_proxy_type) :: Y_proxy
      integer(kind=i_def) :: max_halo_depth_mesh
      type(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      X_proxy = X%get_proxy()
      Y_proxy = Y%get_proxy()
      !
      ! Create a mesh object
      !
      mesh => X_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = X_proxy%vspace%get_last_dof_annexed()
      !
      ! Call kernels and communication routines
      !
      !$omp parallel default(shared), private(df)
      !$omp do schedule(static)
      DO df=loop0_start,loop0_stop
        X_proxy%data(df) = X_proxy%data(df) + real(Y_proxy%data(df),r_def)
      END DO
      !$omp end do
      !$omp end parallel
      !
      ! Set halos dirty/clean for fields modified in the above loop(s)
      !
      CALL X_proxy%set_dirty()
      !
      ! End of set dirty/clean section for above loop(s)
      !
      !
    end subroutine invoke_inc_rdefX_plus_rsolverY

  end module sci_psykal_light_mod
