!* Implementation of the objetcive function and derivatives.

module funcmod
!* Implementation of the objetcive function and derivatives.
implicit none
integer, parameter, public :: dimx = 2

contains
    function fn(x)  !! Objective function
        real*8, dimension(dimx), intent(in) :: x
        real*8 :: fn
        fn = sin(2*x(1)-1) + 2*sin(x(2)-4)
    end function

    function dfn_dx(x)  !! Objective function's derivative
        real*8, dimension(dimx), intent(in) :: x
        real*8, dimension(dimx) :: dfn_dx
        dfn_dx(1) = 2*cos(2*x(1)-1)
        dfn_dx(2) = 2*cos(x(2)-4)
    endfunction
end module