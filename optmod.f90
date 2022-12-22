program optmod
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Implementation of optimization methods
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The way this code is structured:
! - Optimization methods are implemented as subroutines, defined
!   in the "contains" section. Their names are limited to 20
!   characters (``optimmethod`` variable).
! - Objective function definition, as well as its derivatives
!   (when necessary) must be defined in file ``funcmod.f90``.
! - The desired optimization method to be used, as well as its
!   corresponding hyper parameters (see each method's definition)
!   must be declared in the ``main.txt`` file. According to the
!   following formatting:
!
!   > Name of the optimization method (name of the subroutine)
!   > Name of the hyper parameter 1 [tab or space] Hyper
!     parameter 1 value
!   > Name of the hyper parameter 2 [tab or space] 
!     Hyper parameter 2 value
!   > ...
!   > Name of the next optimization method (and so on...)
!
! - Array hyper parameters must have a space or tab separating
!   each element, and the whole array must be contained in one
!   line. Example:
!
!   GD              [Gradient Descent method] 
!   x   1   0   2   [Initial guess hyper parameter]
!   LR  1D-3        [Learning Rate hyper parameter]
!   eps 1D-10       [Relative error hyper parameter]
!   freq 3          [How many this optim. should be repeated]
!   GD              [Another optimization, also with the GD method]
!   x   0   0   0   [Different initial guess]
!   LR  1D-3                
!   eps 1D-10               
!   freq 2                  [Different number of repetitions]
!
! - Subroutine ``handle_main_file(...)`` handles reading of this
!   file, calling the optimization methods and passing the
!   important hyperparameters to them. As of know, they are
!   written to a text file called ``hyperparams.txt``, which is
!   then read by the corresponding optimization method - this can
!   be optimized as described in the ``handle_main_file(...)``
!   subroutine documentation.
!
! - Subroutine ``handle_hyp_file(...)`` handles the reading of
!   the ``hyperparams.txt`` text file, attributing the correct
!   values to the correct lines in an array. It is generic and
!   called by each optimization method, which in turn interprets
!   its results.
!
! - Subroutine ``check_status(...)`` is only defined because it
!   is commonly used code, and serves to handle file reading
!   errors.
!
! - Subroutine ``check_hprparams(...)`` checks if all required
!   hyper parameters were informed by the user, and halts the
!   program if not. It also warns the user if optional hyper
!   parameters were not set, informing their default value. 
!
!   > NOTE: Default values for optional hyper parameters must be
!     defined in the optimization procedure before calling the
!     ``handle_hyp_file(...)`` subroutine, by initializing the
!     associated elements of the hprdata array.
!
! - It is best if optional hyper parameters are defined last (in
!   ``hprnames``), but this isn't a requirement. It is just an
!   organization guideline.
!
! - The arrays passed to ``handle_hyp_file(...)`` are defined as
!   ``allocatable`` in the optimization methods because they
!   depend on the dimension of the objective function's input ---
!   ``dimx`` --- this because one of the hyper parameters is
!   always the initial guess. It could be given special
!   treatment, but this approach generalizes the process and may
!   be useful in the future, in case multidimensional hyper
!   parameters are used in future methods.
!
!       TODO: Add a reporting subroutine
!
!
!
!
!
! ADDING NEW METHODS: Update ``handle_main_file(...)``
! subroutine's ``select case(optim_method)``. The method should
! also support the existing report parameters (check
! ``rprt_params``) in the main program. 
! Also: update nb_methods and method_names

use funcmod
implicit none

! Related to the methods themselves:
integer :: nb_methods=2
character*20, dimension(:), allocatable :: method_names
! Related to report parameters:
integer :: nb_rprt_params=3
character*256, dimension(:), allocatable :: rprt_params

! File reading related:
integer :: filestatus                   ! File reading status (indicates errors)
character*256 :: errormsg               ! Error message and read line
character*20 :: optimmethod             ! Line keyword (hyperparam. name)
integer :: file_number                  ! For generality, may be useful in the future.

allocate(method_names(nb_methods))
allocate(rprt_params(nb_rprt_params))

method_names = (/ "GD     ", &
                & "NR     "/)

rprt_params = (/ "freq  ", &            ! The number of times a method should be run
               & "format", &            ! Formatting string     
               & "print " /)            ! Whether or not to print results in command line
                                            ! This still has to
                                            ! be implemented, but
                                            ! should maybe
                                            ! integrate a few
                                            ! printing options

call handle_main_file(nb_methods, nb_rprt_params, method_names, rprt_params)

contains
    subroutine GD(format_str, output_str)
        ! ~~~~~~~~~~~~~~~~
        ! GRADIENT DESCENT
        ! ~~~~~~~~~~~~~~~~
        ! Gradient Descent method.
        ! This optimization method is based in the fact that the
        ! gradient of a function points in its maximum growth
        ! direction. In order to minimize it, we'll go in the
        ! opposite direction (-grad(f)). Of course, this is done
        ! by updating x (and grad(f)) iteratively.
        ! 
        ! x = x - LR * grad(f(x))
        ! 
        ! HYPER PARAMETERS
        ! ----------------
        ! > x : real*8, dimension(dimx)
        !       Initial guess. Must have as many dimensions as
        !       the input variable has (obviously).
        ! > LR : real*8
        !       Learning Rate. Related to the iteration step.
        ! > eps : real*8
        !       Minimum relative error for stopping. If a
        !       relative step (step divided by the current x) is
        !       smaller than eps, then the process stops and we
        !       consider the optimization as done.
        ! > Nmax : integer, optional. Default = 1,000,000
        !       Maximum number of iterations.
        !

        implicit none
        ! Hyper parameter subroutine:
        integer :: nb_hprparams = 4                                     ! Number of hyperparameters
        integer :: size_hprdata                                         ! Size of hprdata's second dimension
        character*256, dimension(:), allocatable :: hprnames              ! Hyperparameter names.
        real*8, dimension(:, :), allocatable :: hprdata                 ! Where data will be saved   
        integer, dimension(:), allocatable :: hprdim                    ! Dimension of each hprparam.
        logical, dimension(:), allocatable :: hprbool, requiredbool     ! Bool array: hprparam present or not

        ! Actual hyper parameters:
        real*8 :: LR            ! Learning Rate
        real*8 :: eps           ! Minimum relative error for stopping   
        integer :: Nmax         ! Max. number of iterations
        
        ! Actual method:
        integer :: it1, it2
        real*8, dimension(dimx) :: x        ! Input variable
        real*8, dimension(dimx) :: dy_dx    ! Gradient value at x
        real*8 :: y                         ! Function value at x
        real*8, dimension(dimx) :: re       ! Relative error
        logical :: stopbool
        
        ! Report variables:
        character*1000 :: value_str                 ! Actual x and y values for reportting.
        character*1000 :: output_str                ! Actual x and y values for reportting.
        character*1000 :: halting_str               ! Halting reason
        character*3 :: dimx_str                     ! dimx as a string
        character*20 :: format_str                  ! formatting string

        
        !
        size_hprdata = dimx
        allocate(hprnames(nb_hprparams))
        allocate(hprdata(nb_hprparams, size_hprdata))
        allocate(hprdim(nb_hprparams))
        allocate(hprbool(nb_hprparams))
        allocate(requiredbool(nb_hprparams))

        hprnames = (/ "x          ", &
                    & "LR         ", &
                    & "eps        ", &
                    & "Nmax       "/)                               ! Hpr param names
        hprdata = 0                                                 ! (not really necessary)
        hprdim = (/dimx, 1, 1, 1/)                                  ! Dimension of each of the hprparams
        hprbool = .false.                                           ! Initializing hprbool array as .false.
        requiredbool = (/.true., .true., .true., .false./) ! specifying which hprparams are optional
        
        ! DEFAULT VALUES FOR OPTIONAL HYPERPARAMETERS:
        hprdata(4, 1) = 1000000                                     ! Nmax

        call handle_hyp_file(nb_hprparams, size_hprdata, hprnames, hprdata, hprdim, hprbool, requiredbool)

        x = hprdata(1, :hprdim(1))
        LR = hprdata(2, 1)
        eps = hprdata(3, 1)
        Nmax = int(hprdata(4, 1))

        y = fn(x)
        
        call format_data(x, y, format_str, value_str)
        output_str = trim(value_str)
        halting_str=""
        do it1 = 1, Nmax
            dy_dx = dfn_dx(x)
            x = x - LR*dy_dx
            re = abs(LR*dy_dx/x)
            
            do it2 = 1, dimx
                stopbool = .true.
                if (re(it2) .ge. eps) then
                    stopbool = .false.
                    exit
                end if
            end do

            if (stopbool) then
                print *, "Stopping iteration: ", it1
                halting_str = "Precision condition reached;"
                exit
            end if

        end do
        if (halting_str == "") then
            halting_str = "Max. number of iterations reached;"
        endif
        y = fn(x)
        call format_data(x, y, format_str, value_str)
        output_str = trim(output_str)//"    "//trim(value_str)// &
                     & "    "//trim(halting_str)

    end subroutine

    subroutine handle_main_file(nb_methods, nb_rprt_params, method_names, rprt_params)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Subroutine for handling the optimization runs.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! This subroutine reads the ``main.txt`` file and
        ! interpets its contents.
        ! - This means reading the optimization methods, the
        !   report parameters (types of output, number of
        !   times a certain optimization should be run, etc.) and
        !   calling the corresponding optimization subroutine and
        !   passing it the important hyperparameter information.
        ! - As of know, in order to maintain the current
        !   implementation, the method's hyperparameters are
        !   written to the ``hyperparams.txt`` file. This is
        !   surely not the best way to do this (passing the
        !   ``character`` variable directly to the method ought
        !   to be more efficient). This is an improvement that
        !   can be done in the future (changing the
        !   ``handle_hyp_file`` subroutine) and how its called in
        !   the optimization methods.
        ! 
        ! ARGUMENTS
        ! ---------
        ! > nb_methods : integer
        !       The total number of optimization methods
        !       currently implemented
        ! > nb_rprt_params : integer
        !       Number of report parameters currently implemented.
        ! > method_names : character*20, dim(nb_methods)
        !       Array of optimization methods' names/labels.
        ! > rprt_params : character*256, dim(nb_rprt_params)
        !       Array of report parameters' names/labels.
 
        ! This subroutine's arguments:
        integer :: nb_methods, nb_rprt_params
        character*20, dimension(nb_methods) :: method_names
        character*256, dimension(nb_rprt_params) :: rprt_params
        ! Variables related to this subroutines
        integer, dimension(1000) :: line_counter=0                  ! "Saves" line numbers
        integer :: line_number=0, line_idx=0, nb_lines, nb_runs
        character*1000 :: hprparam_str              ! Values for the optim. methods' hpr params
        character*20 :: optim_method                ! Name of the optimization method for a run
        ! Report parameter related
        character*256, dimension(nb_rprt_params) :: rprt_data       ! Values of the rprt params
        integer :: repetitions, rep_iteration                       ! ``freq`` param.
        character*20 :: format_str                                  ! ``format`` param
        ! Report related:
        character*1000 :: hprparam_str2             ! Hyperparameter string
        character*100000 :: rprt_str                ! The report string
        character*1000 :: output_str                ! Opt. method output string
        character*1000 :: hdr_str, hdr_str2         ! Strings for the header
        character*5 :: repetition_str
        integer :: iter_dimx                        ! iterate in dimx
        character*3 :: dim_str                      ! dimension string

        ! File reading related:
        integer :: filestatus                   ! File reading status (indicates errors)
        character*256 :: errormsg, read_line    ! Error message and read line
        character*256 :: key_word                 ! Line keyword (hyperparam. name)
        integer, dimension(1) :: result1        ! Indexes corresponding to optim. methods
        integer, dimension(1) :: result2        ! Indexes corresponding to report params
        integer :: idx1, idx2                   ! The actual indexes
        integer :: main_file, hpr_file, rprt_file   ! For generality, may be useful in the future.


        !! Opening the file
        main_file = 00000
        open(main_file, file="main.txt", status="old", &    ! Opening the file
            & iostat = filestatus, iomsg = errormsg)
        call check_status(filestatus, errormsg)             ! Handling errors
        
        !! Finding where the optim. methods are defined.
        ! The next bit has been adapted from: the following StackOverflow response: 
        ! https://stackoverflow.com/questions/29125581/fortran-find-string-in-txt-file
        do  
            read(main_file, "(a)", iostat=filestatus) read_line     ! Read line
            if (filestatus /= 0) then                           ! When we've reached the end of the file
                exit 
            end if
            line_number = line_number + 1                       ! Increase line number (1st line is line 1).

            read(read_line, *) key_word                         ! Getting keyword (optim. method/param.)

            result1 = findloc(method_names, key_word)           ! Getting its index in the optimization method list
            idx1 = result1(1)                                   ! just for clarity

            if (idx1 == 0) then                                 ! If it isn't an optim. method label, then skip line
                cycle
            end if
            
            ! ELSE: Optimization method name found:
            line_idx = line_idx + 1                             ! Update index (starts at 1)
            line_counter(line_idx) = line_number                ! Save the line number which contains the opt. label
        end do
        nb_runs = line_idx
        nb_lines = line_number
        
        ! We add the end of the file to the list, will be useful
        ! later for halting the code:
        line_counter(nb_runs+1) = nb_lines+1
        
        !! Closing file (only to be reopened again at its beginning)
        close(main_file)                                        ! Closing file
        open(main_file, file="main.txt", status="old", &        ! Opening the file
            & iostat = filestatus, iomsg = errormsg)
        call check_status(filestatus, errormsg)
        
        ! Resetting line idx and number:
        line_idx = 1; line_number = 0;
        
        ! Writting the header for the report:
        rprt_str = "OPT. METHOD;    HYPER PARAMETERS;    SUB-RUN;"
        hdr_str = ""
        hdr_str2 = ""
        do iter_dimx = 1, dimx
            write(dim_str, "(I3)") iter_dimx
            hdr_str = trim(hdr_str)//"    X0("//trim(adjustl(dim_str))//");"
            hdr_str2 = trim(hdr_str2)//"    XF("//trim(adjustl(dim_str))//");"
        end do
        rprt_str = trim(rprt_str)//trim(hdr_str)//"    Y0;"// &
                   & trim(hdr_str2)//"    YF;    Stopping reason;"

        rprt_file = 10000
        open(rprt_file, file="report.txt", &    ! Opening the file
            & iostat = filestatus, iomsg = errormsg)
        call check_status(filestatus, errormsg)             ! Handling errors
        write(rprt_file, '(a)') trim(rprt_str)
        close(rprt_file)

        ! Now actually calling the optimization procedures:
        do  
            
            ! Reading current line
            read(main_file, "(a)", iostat=filestatus) read_line     ! Read line
            if (filestatus /= 0) then                               ! End if end of file
                exit
            end if
            line_number = line_number + 1                           ! Updating line nb.

            ! Tried using the "select case" statement for this,
            ! but apparently cases cannot be written in function
            ! of variales, so...
            if (line_number == line_counter(line_idx)) then         ! If it's the line of an optimization label:
                read(read_line, *) optim_method                     ! Save it as ``optim_method``
                hprparam_str=""                                     ! Reset method's hyperparameters
                hprparam_str2=""                                    ! Reset method's hyperparameters
                rprt_data = (/"1      ", &
                            & "(f11.4)", &
                            & ".false."/)                                  ! Reset report parameters

            else if (line_number == (line_counter(line_idx + 1) - 1)) then  ! If its the last line of a method:
                line_idx = line_idx + 1                     ! Update line_idx
                
                ! Check for report hprparams and adding them to
                ! the hprparam_str        
                read(read_line, *) key_word                 ! This part is the same code as the last (else) condition
                result2 = findloc(rprt_params, key_word)    ! Searching for report parameters
                idx2 = result2(1)

                if (idx2 /= 0) then                         ! If they're present, save.
                    read(read_line, *) key_word, rprt_data(idx2)
                else                                        ! If it's something else, they're optim. hpr params
                    hprparam_str = trim(hprparam_str)//NEW_LINE('A')//read_line
                    hprparam_str2 = trim(hprparam_str2)//", "//read_line
                end if
                
                ! Since this is this run's last line, then we can
                ! save the hyperparameters in a file.
                hpr_file = 00001
                open(hpr_file, file="hyperparams.txt", &                ! Opening the file
                    & iostat = filestatus, iomsg = errormsg)
                call check_status(filestatus, errormsg)                 ! Just in case
                write(hpr_file, '(a)') trim(hprparam_str(2:))
                close(hpr_file)
                
                ! Reading report parameters:
                read(rprt_data(1), *) repetitions
                read(rprt_data(2), "(a)") format_str

                do rep_iteration = 1, repetitions
                    select case(optim_method)
                        case("GD")
                            call GD(format_str, output_str)
                    end select

                    write(repetition_str, '(I5)') rep_iteration
                    rprt_str = trim(optim_method)//";    "//trim(hprparam_str2(3:)) &
                               & //";    "//trim(adjustl(repetition_str))//";    "//trim(output_str)
                    
                    ! Adding results to the report:
                    open(rprt_file, file="report.txt", position="append", &
                        & iostat = filestatus, iomsg = errormsg)
                    call check_status(filestatus, errormsg)
                    write(rprt_file, '(a)') trim(rprt_str)
                    close(rprt_file)
                enddo
            else
               ! Check for report hprparams and adding them to
                ! the hprparam_str        
                read(read_line, *) key_word                 ! This part is the same code as the last (else) condition
                result2 = findloc(rprt_params, key_word)    ! Searching for report parameters
                idx2 = result2(1)
                if (idx2 /= 0) then                         ! If they're present, save.
                    read(read_line, *) key_word, rprt_data(idx2)
                else                                        ! If it's something else, they're optim. hpr params
                    hprparam_str = trim(hprparam_str)//NEW_LINE('A')//read_line
                    hprparam_str2 = trim(hprparam_str2)//", "//read_line
                end if
            end if
        end do

        close(main_file)                                            ! Closing file   

    end subroutine

    subroutine handle_hyp_file(nb_hprparams, size_hprdata, hprnames, hprdata, hprdim, hprbool, requiredbool)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Subroutine for handling the method's hyper parameters.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! This subroutine reads the ``hyperparams.txt`` and
        ! stores the values informed by the user. 
        ! - It reads the file line by line. When the first word
        !   in a line matches with one of the hyper parameter's
        !   name (contained in ``hprnames``, defined in the
        !   method of choice), it saves its value in the
        !   ``hprdata`` array.
        ! 
        ! ARGUMENTS
        ! ---------
        ! > np_hprparams : integer
        !       The total number of hyper parameters possible for the
        !       optimization method.
        ! > size_hprdata : integer
        !       Size of hprdata second dimension (usually dimx,
        !       extended for generality).
        ! > hprnames : character*256, dim(nb_hprparams)
        !       Array of hyper parameter names/labels.
        ! > hprdata : real*8, dim(nb_hprparams, size_hprdata)
        !       Array of data (each hyperparameter value). Will
        !       have a lot of unused elements, since not all
        !       hyper parameters will have the same dimension.
        !       Each line is a different hyper parameter, in the
        !       same order as in hprnames.
        ! > hprdim : integer, dim(nb_hprparams)
        !       Dimension of each hyper parameter, again in order.
        ! > hprbool : logical, dim(nb_hprparams)
        !       Boolean array indicating which hyper parameters
        !       have been informed by the user (.true.) and which
        !       weren't (.false.).
        ! > requiredbool : logical, dim(nb_hprparams)
        !       Which hyper parameters are required (.true.) and
        !       which are optional (.false.). This argument is
        !       not used in this subroutine, but directly in
        !       ``check_hprparams``, which is called directly at
        !       the end of this subroutine, but was separated
        !       solely for clarity.

        ! This subroutine:
        integer :: nb_hprparams, size_hprdata
        character*256, dimension(nb_hprparams) :: hprnames
        real*8, dimension(nb_hprparams, size_hprdata) :: hprdata
        integer, dimension(nb_hprparams) :: hprdim
        logical, dimension(nb_hprparams) :: hprbool
        ! Hyper parameter check subroutine:
        logical, dimension(nb_hprparams) :: requiredbool    ! Bool array: hprparam optional or not

        ! File reading related:
        integer :: filestatus                   ! File reading status (indicates errors)
        character*256 :: errormsg, read_line    ! Error message and read line
        character*5 :: key_word                 ! Line keyword (hyperparam. name)
        integer, dimension(1) :: result         ! Indexes corresponding to the hyperparam. array
        integer :: idx                          ! The actual index
        integer :: file_number                  ! For generality, may be useful in the future.


        !! Opening the file
        file_number = 00001
        open(file_number, file="hyperparams.txt", status="old", &     ! Opening the file
            & iostat = filestatus, iomsg = errormsg)
        
        call check_status(filestatus, errormsg)                 ! When we've reached the end of the file
        
        !! Reading the hyperparameters:
        ! The next bit has been adapted from: the following StackOverflow response: 
        ! https://stackoverflow.com/questions/29125581/fortran-find-string-in-txt-file
        do
            read(file_number, "(a)", iostat=filestatus) read_line     ! Read line
            if (filestatus /= 0) then                           ! Error handling
                exit 
            end if
            read(read_line, *) key_word                         ! Getting keyword (hyperparam.)
            result = findloc(hprnames, key_word)                ! Getting its index
            idx = result(1)                                     ! just for clarity

            ! NOTE: Fortran runtime error: End of file will ensue
            ! if the hyperparameter informed has less dimensions
            ! than necessary.
            read (read_line, *) key_word, hprdata(idx, :hprdim(idx))
            hprbool(idx) = .true.
        end do
        
        !! Closing file
        close(file_number)                                            ! Closing file

        !! Checking required and optional hyperparameters:
        call check_hprparams(nb_hprparams, hprnames, hprdata, hprdim, hprbool, requiredbool)
    end subroutine

    subroutine check_status(filestatus, errormsg)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Subroutine for handling file reading errors.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! This subroutine takes directly the ``iostat`` and
        ! ``iomsg`` values from the file opening procedure.
        !
        ! ARGUMENTS
        ! ---------
        ! > filestatus : integer
        !       The ``iostat`` value.
        ! > errormsg : character*256
        !       The associated error message.
        !

        integer :: filestatus                   ! File reading status (indicates errors)
        character*256 :: errormsg               ! Error message and read line

        if (filestatus /= 0) then                               ! Handling error
            print *, "Error when reading the hyperparameter file:", &
                     & NEW_LINE('A')//" -> ", errormsg
            stop
        endif
    end subroutine

    subroutine check_hprparams(nb_hprparams, hprnames, hprdata, hprdim, hprbool, requiredbool)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Subroutine for checking which hyper parameters were informed by the user.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! This subroutine starts where the ``handle_hyp_file(...)`` 
        ! subroutine stops, and is written separately only for
        ! clarity. It checks if all required hyper parameters
        ! have been defined (and will halt execution if not). If
        ! optional hyper parameters haven't been specified, then
        ! it warns the user and shows them the default values
        ! used.
        !
        !   > NOTE: Default values for optional hyper parameters must be
        !     defined in the optimization procedure before calling the 
        !     ``handle_hyp_file(...)`` subroutine, by initializing the 
        !     associated elements of the hprdata array.
        !
        ! - ARGUMENTS are explained in the ``handle_hyp_file(...)``
        !   subroutine's documentation


        integer :: nb_hprparams, idx
        ! Hyperparameter arrays
        character*256, dimension(nb_hprparams) :: hprnames    ! Hyperparameter names.
        real*8, dimension(nb_hprparams, dimx) :: hprdata    ! Where data will be saved  
        integer, dimension(nb_hprparams) :: hprdim          ! Dimension of each hprparam.
        logical, dimension(nb_hprparams) :: hprbool         ! Bool array: hprparam present or not
        logical, dimension(nb_hprparams) :: requiredbool    ! Bool array: hprparam optional or not
        ! For checking the hyperparams:
        character*256 :: uninformed_required, uninformed_optional   ! Hyperparameters nor informed by the user
        character*256 :: default_value                              ! Default values of optional hprparams
        logical :: haltbool = .false.                               ! Whether the program should be halted or not

        !! Checking if all hyperparameters were informed:
        uninformed_required = ""
        uninformed_optional = ""
        do idx = 1, nb_hprparams
            if (.not. hprbool(idx)) then            ! If it wasn't informed
                if (requiredbool(idx)) then         ! And it was required
                    haltbool = .true.               ! Program will be halted
                    uninformed_required = trim(uninformed_required)//", "//trim(hprnames(idx))
                        ! Adding to the list of missing parameters
                                                    
                else                                ! If it was optional, then we'll warn the user
                    write(default_value, *) hprdata(idx, :hprdim(idx))
                    uninformed_optional = trim(uninformed_optional)//" -> "// &
                        & trim(hprnames(idx))// &               ! Name of the optional hprparam
                        & " = "//trim(default_value)// &        ! Default value of the optional hprparam
                        & NEW_LINE('A')  
                end if
            end if
        end do
        if (haltbool) then                          ! If a required hyperparam was missing: Halt code
            uninformed_required = uninformed_required
            print *, "The following REQUIRED hyperparameters were not informed by the user:"
            print *, " -> "//uninformed_required(3:)    ! 3: to take out the ", "
            stop
        end if
        if (uninformed_optional /= "") then         ! If an optional parameter was missing: Warn user
            print *, "The following OPTIONAL hyperparameters were not informed by the user &
                     &and their default values will be used."
            print *, " "//trim(uninformed_optional)
        end if
    end subroutine
    
    subroutine format_data(x, y, format_str, value_str)
        
        real*8, dimension(dimx), intent(in) :: x        ! Input variable
        real*8, intent(in) :: y                         ! Function value at x

        integer :: iter_dim             ! Iterator on x dimension
        character*100 :: var_str        ! Variable string
        character*20 :: format_str      ! Formatting string
        character*1000 :: value_str    ! Output string
        
        value_str = ""
        do iter_dim = 1, dimx
            write(var_str, format_str) x(iter_dim)
            value_str = trim(value_str)//"    "//trim(adjustl(var_str))//";"
        end do
        write(var_str, format_str) y
        value_str = trim(adjustl(value_str))//"    "//trim(adjustl(var_str))//";"

    end subroutine
end program