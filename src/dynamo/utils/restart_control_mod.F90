!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief   NML list IO
!!
!>  @details Primarily for checkpoint/restart to be replaced by cfg obj
!!           Reads a namelist with the data for controlling whether prognostic
!!           fields are read and or written, how many timesteps and diagnostic
!!           frequency.
!!           restart type data attributes are all private so it is used
!!           by calling procedures
!-------------------------------------------------------------------------------
module restart_control_mod
  use constants_mod, only : str_max_filename, str_long
  use log_mod,       only : log_event, log_scratch_space, &
                            LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  implicit none

  private

  type, public :: restart_type
     private
     logical :: checkpoint_read 
     integer :: timestep_start
     integer :: timestep_end 
     logical :: checkpoint_write
     character(len=str_max_filename) :: restart_stem_name 
     integer :: diagnostic_frequency
   contains

  !> @brief Returns the filename for start-1 + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] field_name character which holds the name of the field
     procedure :: startfname

  !> @brief Returns the filename for end + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] field_name character which holds the name of the field
     procedure :: endfname

  !> @brief Returns an integer the starting timestep 
  !> @param[in] self The calling restart_type
     procedure :: ts_start

  !> @brief Returns an integer the end timestep
  !> @param[in] self The calling restart_type
     procedure :: ts_end

  !> @brief Returns a logical, whether to read a checkpoint
  !> @param[in] self The calling restart_type
     procedure :: read_file

  !> @brief Returns a logical, whether to write a checkpoint
  !> @param[in] self The calling restart_type
     procedure :: write_file

  !> @brief Returns an integer the frequency of diagnostic measurement
  !> @param[in] self The calling restart_type
     procedure :: diag_freq

  !> @brief Returns the filename for ts + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] ts_fname character which holds the name of the field
  !> @param[in] ts integer the timestep
     procedure :: ts_fname
  end type restart_type

  interface restart_type
     module procedure restart_constructor
  end interface restart_type

contains

  !> @brief Consructor for the restart type
  !> @param[in] character fname the filename of the namelist file
  !> @detail Opens, reads and closes namelist file (with error reporting)
  !! sets the values of the data attributes from the namelist file
  !! does some sanity checking to help prevent user error for timestep start
  !! and end values.
  function restart_constructor(fname) result(self)
    implicit none
    character(len=str_max_filename), intent(in) :: fname ! the name of the nml file
    type(restart_type) :: self
    logical :: checkpoint_read, checkpoint_write
    integer :: timestep_start
    integer :: timestep_end    
    integer :: diagnostic_frequency
    integer :: ierr, funit
    character(len=str_long) :: ioerrmsg=''
    character(len=str_max_filename) :: restart_stem_name

    namelist /restart_nml/ checkpoint_read, timestep_start, & 
                           timestep_end, checkpoint_write,  &
                           restart_stem_name, diagnostic_frequency

    funit=555
    open(funit,file=trim(fname),iostat=ierr,status='old',iomsg=ioerrmsg)
    if(ierr/=0) then
       write( log_scratch_space,'(A)') "problems opening File..."
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call log_event(ioerrmsg,LOG_LEVEL_ERROR)
    end if
    
    read(555,nml=restart_nml,iostat=ierr,iomsg=ioerrmsg)
    if(ierr/=0) then
       write(*,*) checkpoint_read
       write(*,*) timestep_start
       write(*,*) timestep_end
       write(*,*) checkpoint_write
       write(*,*) trim(restart_stem_name)
       write(*,*) diagnostic_frequency
       write( log_scratch_space,'(A)') "problems reading File ..."
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call log_event(ioerrmsg,LOG_LEVEL_ERROR)
    end if

    close(555,iostat=ierr,iomsg=ioerrmsg)
    if(ierr/=0) then
       write( log_scratch_space,'(A,A)') "Closing File:",fname
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call log_event(ioerrmsg,LOG_LEVEL_ERROR)
    end if

    self%checkpoint_read = checkpoint_read
    self%timestep_start = timestep_start
    self%timestep_end = timestep_end
    self%checkpoint_write = checkpoint_write
    self%restart_stem_name = restart_stem_name
    self%diagnostic_frequency = diagnostic_frequency

    ! do some sanity checks 
    if(self%timestep_start >= self%timestep_end) then
       write(log_scratch_space,'(A,I6,",",I6)') & 
            "restart_control: timestep_start must be less than timestep_end:", &
            self%timestep_start,self%timestep_end
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    
    if(self%timestep_start<=0) then
       write(log_scratch_space,'(A,I6)') & 
            "restart_control: timestep_start cannot be negative or zero:",  &
            self%timestep_start
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if

    if(self%timestep_end<0) then
       write(log_scratch_space,'(A,I6)')    &
            "restart_control: timestep_end cannot be negative", self%timestep_end
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if

    if(self%timestep_start>=1000000) then
       write(log_scratch_space,'(A,I8)')  &
            "restart_type:Whoa, I can't count that far:",self%timestep_start
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if

    if(self%timestep_end>=1000000) then
       write(log_scratch_space,'(A,I8)')  &
            "restart_type:Whoa, I can't count that far:",self%timestep_end
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if
    if(self%checkpoint_read) then
       write(log_scratch_space,'(A,A)') "Re-starting Dynamo:", &
            trim(self%startfname(""))
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
    else 
       call log_event("Cold start, spinning up",LOG_LEVEL_INFO)
    end if
    
  end function restart_constructor

  !> @brief Returns the filename for start-1 + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] field_name character which holds the name of the field
  function startfname(self,field_name)
    class(restart_type), intent(in) :: self
    character(len=*),    intent(in) :: field_name
    character(len=str_max_filename) :: startfname
    integer :: ts_minus_1
    ts_minus_1 = self%timestep_start - 1

    write(startfname,'(A,A,A,A,I6.6)') trim(self%restart_stem_name),"_", &
         trim(field_name),".T",ts_minus_1
  end function startfname

  !> @brief Returns the filename for end + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] field_name character which holds the name of the field
  function endfname(self,field_name)
    class(restart_type), intent(in) :: self
    character(len=*),    intent(in) :: field_name
    character(len=str_max_filename) :: endfname

    write(endfname,'(A,A,A,A,I6.6)') trim(self%restart_stem_name),"_", &
         trim(field_name),".T",self%timestep_end

  end function endfname

  !> @brief Returns the filename for ts + string for filed name
  !> @param[in] self The calling restart_type
  !> @param[in] ts_fname character which holds the name of the field
  !> @param[in] ts integer the timestep
  function ts_fname(self,field_name,ts)
    class(restart_type), intent(in) :: self
    character(len=*),    intent(in) :: field_name
    integer,             intent(in) :: ts
    character(len=str_max_filename) :: ts_fname

    write(ts_fname,'(A,A,A,A,I6.6)') trim(self%restart_stem_name),"_", &
         trim(field_name),".T",ts

  end function ts_fname

  !> @brief Returns an integer the starting timestep 
  !> @param[in] self The calling restart_type
  function ts_start(self)
    class(restart_type), intent(in) :: self
    integer                         :: ts_start
    ts_start = self%timestep_start 
  end function ts_start

  !> @brief Returns an integer the end timestep
  !> @param[in] self The calling restart_type
  function ts_end(self)
    class(restart_type), intent(in) :: self
    integer                         :: ts_end
    ts_end = self%timestep_end
  end function ts_end

  !> @brief Returns an integer the frequency of diagnostic measurement
  !> @param[in] self The calling restart_type
  function diag_freq(self)
    class(restart_type), intent(in) :: self
    integer                         :: diag_freq
    diag_freq = self%diagnostic_frequency
  end function diag_freq

  !> @brief Returns a logical, whether to read a checkpoint
  !> @param[in] self The calling restart_type
  function read_file(self) 
    class(restart_type), intent(in) :: self
    logical read_file
    read_file = self%checkpoint_read
  end function read_file

  !> @brief Returns a logical, whether to write a checkpoint
  !> @param[in] self The calling restart_type
  function write_file(self) 
    class(restart_type), intent(in) :: self
    logical write_file
    write_file = self%checkpoint_write
  end function write_file
  
end module restart_control_mod
