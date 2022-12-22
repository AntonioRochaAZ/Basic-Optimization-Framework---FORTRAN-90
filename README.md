# Basic Optimization Framework - FORTRAN 90
 Small project for relearning FORTRAN 90 and studying optimization concepts. This project is a generic framework over which simple optimization methods can be built upon and run.

> NOTE: this project is finished.

## Implementation of optimization methods

The way this code is structured:
- Optimization methods are implemented as subroutines, defined
  in the "contains" section. Their names are limited to 20
  characters (``optimmethod`` variable).
- Objective function definition, as well as its derivatives
  (when necessary) must be defined in file ``funcmod.f90``.
  Looking at the source file should give you a good idea of how
  it should be defined.
- The desired optimization method to be used, as well as its
  corresponding hyper parameters (see each method's definition)
  must be declared in the ``main.txt`` file. According to the
  following formatting:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  > Name of the optimization method (name of the subroutine)
  > Name of the hyper parameter 1 [tab or space] Hyper parameter 1 value
  > Name of the hyper parameter 2 [tab or space] Hyper parameter 2 value
  > ...
  > Name of the next optimization method (and so on...)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Array hyper parameters must have a space or tab separating each element, and the whole array must be contained in one line. Example:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  GD              [Gradient Descent method] 
  x   1   0   2   [Initial guess hyper parameter]
  LR  1D-3        [Learning Rate hyper parameter]
  eps 1D-10       [Relative error hyper parameter]
  freq 3          [How many this optim. should be repeated]
  GD              [Another optimization, also with the GD method]
  x   0   0   0   [Different initial guess]
  LR  1D-3                
  eps 1D-10               
  freq 2                  [Different number of repetitions]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Running the optimization
With gfortran:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> gfortran -c funcmod.f90 optmod.f90
> gfortran funcmod.o optmod.o -o optim.out
> ./optim.out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- The first command compiles ``funcmod.f90``and ``optmod.f90``,
generating files ``funcmod.o``, ``funcmod.mod`` and ``optmod.o`` in the same
folder.
- The second command turns them both into a single executable here named as ``optim.out``.
- The third command runs the executable. The intermediary file
  ``hyperparams.txt`` is created and is not automatically deleted
  after runtime (although it could). The ``report.txt`` file is
  generated with the optimization results. The current report is
  simple, but should be relatively easily expandable in the
  context of this framework.

  # "Documentation"

- Subroutine ``handle_main_file(...)`` handles reading of this
  file, calling the optimization methods and passing the
  important hyperparameters to them. As of know, they are
  written to a text file called ``hyperparams.txt``, which is
  then read by the corresponding optimization method - this can
  be optimized as described in the ``handle_main_file(...)``
  subroutine documentation.

- Subroutine ``handle_hyp_file(...)`` handles the reading of
  the ``hyperparams.txt`` text file, attributing the correct
  values to the correct lines in an array. It is generic and
  called by each optimization method, which in turn interprets
  its results.

- Subroutine ``check_status(...)`` is only defined because it
  is commonly used code, and serves to handle file reading
  errors.

- Subroutine ``check_hprparams(...)`` checks if all required
  hyper parameters were informed by the user, and halts the
  program if not. It also warns the user if optional hyper
  parameters were not set, informing their default value. 

  > NOTE: Default values for optional hyper parameters must be
    defined in the optimization procedure before calling the
    ``handle_hyp_file(...)`` subroutine, by initializing the
    associated elements of the hprdata array.

- It is best if optional hyper parameters are defined last (in
  ``hprnames``), but this isn't a requirement. It is just an
  organization guideline.

- The arrays passed to ``handle_hyp_file(...)`` are defined as
  ``allocatable`` in the optimization methods because they
  depend on the dimension of the objective function's input ---
  ``dimx`` --- this because one of the hyper parameters is
  always the initial guess. It could be given special
  treatment, but this approach generalizes the process and may
  be useful in the future, in case multidimensional hyper
  parameters are used in future methods.
    
    TODO: Add a reporting subroutine

> FOR ADDING NEW METHODS: Update ``handle_main_file(...)``
subroutine's ``select case(optim_method)``. The method should
also support the existing report parameters (check
``rprt_params``) in the main program. 
Also: update nb_methods and method_names