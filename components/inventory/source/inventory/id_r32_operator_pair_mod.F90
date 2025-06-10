!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair operators with a unique identifier.
module id_r32_operator_pair_mod

  use constants_mod,         only: i_def
  use function_space_mod,    only: function_space_type
  use operator_real32_mod,   only: operator_real32_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Operator Pair
  ! ========================================================================== !

  !> @brief An object pairing an operator with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r32_operator_pair_type

    private

    type(operator_real32_type) :: operator_

  contains

    procedure, public :: initialise
    procedure, public :: copy_initialise
    procedure, public :: get_operator

    final :: destructor

  end type id_r32_operator_pair_type

contains

  !> @brief Initialises the id_r32_operator_pair object with a new operator
  !> @param[in] fs_target The function space of the target field of the operator
  !> @param[in] fs_source The function space of the source field of the operator
  !> @param[in] id        The integer ID to pair with the operator
  subroutine initialise(self, fs_target, fs_source, id)

    implicit none

    class(id_r32_operator_pair_type),   intent(inout) :: self
    type(function_space_type), pointer, intent(in)    :: fs_target
    type(function_space_type), pointer, intent(in)    :: fs_source
    integer(kind=i_def),                intent(in)    :: id

    call self%operator_%initialise(fs_target, fs_source)
    call self%set_id(id)

  end subroutine initialise

  !> @brief Initialises the id_r32_operator_pair object by copying in an operator
  !> @param[in] operator_in  The operator that will be stored in the paired object
  !> @param[in] id           The integer ID to pair with the operator
  subroutine copy_initialise(self, operator_in, id)

    implicit none

    class(id_r32_operator_pair_type), intent(inout) :: self
    type(operator_real32_type),       intent(in)    :: operator_in
    integer(kind=i_def),              intent(in)    :: id

    self%operator_ = operator_in%deep_copy()
    call self%set_id(id)

  end subroutine copy_initialise

  !> @brief Get the operator corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The operator
  function get_operator(self) result(operator_out)

    implicit none

    class(id_r32_operator_pair_type), target, intent(in) :: self
    type(operator_real32_type),               pointer    :: operator_out

    operator_out => self%operator_

  end function get_operator

  !> @breif Calls finaliser on field owned by the instance
  subroutine destructor(self)
    implicit none
    type(id_r32_operator_pair_type), intent(inout) :: self
    call self%operator_%operator_final()
  end subroutine destructor

end module id_r32_operator_pair_mod
