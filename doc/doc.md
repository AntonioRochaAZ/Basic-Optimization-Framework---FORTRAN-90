---
project: Basic Optimization Framework - FORTRAN 90
author: AntonioRochaAZ
src_dir: ../src
preprocess: False
exclude_dir: ./objs
project_github: https://github.com/AntonioRochaAZ/Basic-Optimization-Framework---FORTRAN-90
github: https://github.com/AntonioRochaAZ/
website: https://antoniorochaaz.github.io/
favicon: favicon.ico
graph: true
---

Small project for learning modern FORTRAN 90 and studying optimization concepts. This project is a generic framework over which simple optimization methods can be built upon and run.

**NOTE: this project is finished.**

# Code structure:

The way this code is structured:
- Optimization methods are implemented as subroutines, defined
  in the "contains" section of module [[optmod]]. Their names are limited to 20
  characters (``method_names`` variable).
- Objective function definition, as well as its derivatives
  (when necessary) must be defined in file ``funcmod.f90``, in the [[funcmod]]
  module.
- The desired optimization method to be used, as well as its
  corresponding hyper parameters (see each method's definition)
  must be declared in the ``main.txt`` file. According to the
  following formatting:

```
   > Name of the optimization method (name of the subroutine)
   > Name of the hyper parameter 1 [tab or space] Hyper
     parameter 1 value
   > Name of the hyper parameter 2 [tab or space] 
     Hyper parameter 2 value
   > ...
   > Name of the next optimization method (and so on...)
```

- Array hyper parameters must have a space or tab separating
  each element, and the whole array must be contained in one
  line. Example:

```
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
```

# How the code works

The mentioned procedures are all defined in the [[optmod]] module.

- Subroutine ``handle_main_file(...)`` handles reading of this
  file, calling the optimization methods and passing the
  important hyperparameters to them. As of know, they are
  written to a text file called ``output/hyperparams.txt`` (deleted
  by the end of execution), which is
  then read by the corresponding optimization method - this could
  be eventually optimized as described in the ``handle_main_file(...)``
  subroutine documentation.

- Subroutine ``handle_hyp_file(...)`` handles the reading of
  the ``output/hyperparams.txt`` text file, attributing the correct
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

**NOTE**: Default values for optional hyper parameters must be
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

# How to add to it

ADDING NEW METHODS: Other than updating [[optmod]]
variables relating to the number of methods and their names
(``nb_methods`` and ``method_names``), one
must then update ``handle_main_file(...)``
subroutine's ``select case(optim_method)`` with the name
of the new method. The method should
also support the existing report parameters (check
``rprt_params``) in the main program. 

TODO: Add a reporting subroutine