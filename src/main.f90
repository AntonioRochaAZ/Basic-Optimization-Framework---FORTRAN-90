!* Implementation of program (just calls [[handle_main_file]]).

program main
use optmod
implicit none

call handle_main_file(nb_methods, nb_rprt_params, method_names, rprt_params)
    
end program main