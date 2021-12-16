module roots
  use iso_fortran_env, only: real32
  ! This module contains methods for finding roots of linear and non-
  ! linear equations.
  ! Methods:
  !          root_searchD - uses exhaustive searching together with
  !                         "hybridD" routine to locate all roots. 
  !          hybridD - uses either the Newton-Raphson or Bisection method
  !                   to find a root.
  !
  !          root_search - uses exhaustive searching together with
  !                         "hybridSA" routine to locate all roots. 
  !          hybridSA - uses either the accelerated Secant or Bisection method
  !                     to find a root.
  implicit none
  
  private
  
  public :: root_searchD, hybridD          ! routines if the derivative is available
  public :: root_search, hybridSA          ! routines if the derivative is unavailable
  
    
contains

  
  subroutine root_searchD(lbound, rbound, tol, f, df, max_iters, step)
    ! Uses exhausive searching for bracketing all the roots within a user-defined range,
    ! and then calls the hybridD routine to find the root in each bracket.  The step-size
    ! for the search is an optional user input.  If none is supplied, 0.1 is used. The
    ! maximum number of iterations is also optional. If none is supplied, 30 iterations
    ! are allowed.
    implicit none
    
    ! arguments and interfaces
    real(real32), intent(inout) :: lbound, rbound
    real(real32), intent(in) :: tol
    
    interface
       pure function f(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function f
       
       pure function df(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function df
    end interface
    
    ! optional arguments
    integer, intent(in), optional :: max_iters
    real(real32), intent(in), optional :: step
    
    ! locals
    integer :: i, iters, max_iterations, root_number, total_steps
    real(real32) :: lbrac, rbrac, root, step_size
    logical :: root_present = .false.
    
    ! begin
    max_iterations = 30; step_size = 0.1d0
    root_number = 0
    
    ! Check for optional arguments
    if(present(max_iters)) max_iterations = max_iters
    if(present(step)) step_size = step
    
    ! Perform exhaustive search
    ! ensure steps are not too small
    if(step_size <= tol) goto 100
    
    total_steps = int((rbound-lbound)/step_size)
    
    print *, 'Total steps:', total_steps
    
    do i = 1, total_steps-1
       ! Select search interval for root i
       lbrac = lbound+step_size*i
       rbrac = lbound+step_size*(i+1)
       ! Verify that the root is bracketed
       if (f(lbrac)*f(rbrac) < 0.d0) then
          root_present = .true.
          root_number = root_number + 1
          ! Find the root within interval
          call hybridD(lbrac, rbrac, root, tol, iters, max_iters, f, df)
          ! Print results
          if (iters .ge. max_iters) then
             print '(a)', "Maximum iterations reached."
             print '(a,i2,a,f11.8,a,i3,a)', 'Last approximation to the root', root_number, &
                  ' found at x =', root, ' after', max_iters,' iterations.'
          else
             print '(a,i3,a,f11.8,a,i3,a)', 'Root', root_number, ' found at x =', root, &
                  ' in', iters,' iterations.'
          end if
       end if 
    end do
    
    if (root_present .eqv. .false.) goto 200
    
    ! normal exit
    return
    
    ! error trap
100 print '(a,f11.8,a)', 'The specified step-size dx =', step_size, ' is too small.'
    return
200 print '(a)', "No roots are bracketed."
    return
    
  end subroutine root_searchD
  
  
  subroutine hybridD(lbound, rbound, root, tolerance, iterations, max_iterations, &
       f, df)
    ! A hybrid Bisection/Newton-Raphson method for finding the root of a
    ! function F, with derivative Fprime. The root is known to be between "lbound"
    ! and "rbound", and the result is returned in "root". If the next NR guess
    ! is within the known bounds, the step is accepted; otherwise, a bisection
    ! step is taken.
    implicit none
    
    ! arguments and interfaces
    real(real32), intent(inout) :: lbound, rbound, root
    real(real32), intent(in) :: tolerance
    integer, intent(out) :: iterations
    integer, intent(in) :: max_iterations
    
    interface
       pure function f(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function f
       
       pure function df(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function df
    end interface
    
    ! locals
    real(real32) :: delta, dfxmid, error, fxmid, fx0, fx1, x0, x1, xmid
    character(len=6) :: type
    
    ! Verify that the root is bracketed
    if (f(lbound) * f(rbound) .gt. 0.d0) goto 100

    ! Initialization
    x0 = lbound; x1 = rbound
    fx0 = f(x0)
    fx1 = f(x1)
    if (abs(fx0) .le. abs(fx1)) then
       xmid = x0
       fxmid = fx0
    else
       xmid = x1
       fxmid = fx1
    end if
    dfxmid = df(xmid)
    
    ! Begin root finding
    do iterations = 1, max_iterations
       ! Determine if NR or Bisection step
       if ((dfxmid * (xmid-x0) - fxmid) * &
            (dfxmid * (xmid-x1) -fxmid) .le. 0) then
          ! OK to take a NR step
          delta = -fxmid / dfxmid
          xmid = xmid + delta
          type = '(NR)'
       else
          ! take a bisection step instead
          delta = (x1 - x0) / 2.d0
          xmid = (x0 + x1) / 2.d0
          type = '(BS)'
       end if
       
       ! Compare the relative error ti the Tolerance.
       if (abs(delta/xmid) .le. tolerance) then
          ! Error is BELOW tolerance, the root has been found,
          ! so exit loop.
          exit
       else
          ! The relative error is too big, loop again.
          fxmid = f(xmid)
          dfxmid = df(xmid)
          ! Adjust brackets.
          if (fx0 * fx1 .le. 0) then
             ! The root is in left subinterval...
             x1 = xmid
             fx1 = fxmid
          else
             ! The root is in right subinterval
             x0 = xmid
             fx0 = fxmid
          end if
       end if
       print '(a,i2,a8,a10,f12.8)', 'Iteration:', iterations, type, 'Root: ', xmid
    end do

    ! normal exit
    root = xmid
    return
    
    ! error trap
100 print '(a)', "Root not bracketed."
    return
  end subroutine hybridD


  subroutine root_search(lbound, rbound, tol, f, max_iterations, step)
    ! Uses exhausive searching for bracketing all the roots within a user-defined range,
    ! and then calls the "hybridD" routine to find the root in each bracket, where no
    ! derivative is used.  The step-size for the search is an optional user input.  If
    ! none is supplied, 0.1 is used. The maximum number of iterations is also optional.
    ! If none is supplied, 30 is used.
    implicit none

    ! arguments and interfaces
    real(real32), intent(in) :: lbound, rbound, tol

    interface
       pure function f(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function f
    end interface
    
    ! optional arguments
    integer, intent(in), optional :: max_iterations
    real(real32), intent(in), optional :: step

    ! locals
    integer :: i, iters, max_iters, root_number, total_steps
    real(real32) :: lbrac, rbrac, root, step_size
    logical :: root_present = .false.

    ! begin
    max_iters = 30; step_size = 0.1d0
    root_number = 0
    
    ! Check for optional arguments
    if(present(max_iterations)) max_iters = max_iterations
    if(present(step)) step_size = step

    ! Perform exhaustive search
    ! ensure steps are not too small
    if(step_size <= tol) goto 100

    total_steps = int((rbound-lbound)/step_size)

    print *, 'Total steps:', total_steps

    do i = 1, total_steps-1
       ! Select search interval for root i
       lbrac = lbound+step_size*i
       rbrac = lbound+step_size*(i+1)
       ! Verify that the root is bracketed
       if (f(lbrac)*f(rbrac) < 0.d0) then
          root_present = .true.
          root_number = root_number + 1
          ! Find the root within interval
          call hybridSA(lbrac, rbrac, root, tol, iters, max_iters, f)
          ! Print results
          if (iters .ge. max_iters) then
             print '(a)', "Maximum iterations reached."
             print '(a,i2,a,f11.8,a,i3,a)', 'Last approximation to the root', root_number, &
                  ' found at x =', root, ' after', max_iters,' iterations.'
          else
             print '(a,i3,a,f11.8,a,i3,a)', 'Root', root_number, ' found at x =', root, &
                  ' in', iters,' iterations.'
          end if
       end if
    end do
    
    if (root_present .eqv. .false.) goto 200

    ! normal exit
    return

    ! error trap
100 print '(a,f11.8,a)', 'The specified step-size dx =', step_size, ' is too small.'
    return
200 print '(a)', "No roots are bracketed."
    return        
  end subroutine root_search


  
  subroutine hybridSA(lbound, rbound, root, tolerance, iterations, max_iterations, f)
    ! A hybrid Bisection/Secant method for finding the root of a
    ! function F when the derivative of F is not available. The Secant method
    ! uses acceleration to speed converdenge. The root is known
    ! to be between "lbound" and "rbound", and the result is returned in "root".
    ! If the next SA guess is within the known bounds, the step is accepted;
    ! otherwise, a bisection step is taken.
    implicit none
  
    ! arguments and interfaces
    real(real32), intent(inout) :: lbound, rbound, root
    real(real32), intent(in) :: tolerance
    integer, intent(out) :: iterations
    integer, intent(in) :: max_iterations
    
    interface
       pure function f(x) result(answer)
         import
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function f
    end interface
    
    ! locals
    real(real32) :: a, b, c, d, dx, error
    real(real32) :: x0, x1, x2, x3, f0, f1, f2, f3
    character(len=8) :: type
    
    ! Verify that the root is bracketed
    if (f(lbound) * f(rbound) .gt. 0.d0) goto 100
    
    ! Initialization
    x0 = lbound; x1 = rbound
    x2 = (x0 + x1) / 2.d0
    f0 = f(x0)
    f1 = f(x1)
    f2 = f(x2)
    
    a = ((f0-f1) * (x2-x1) - (f2-f1) * (x0-x1)) &
         / ((x0-x2) * (x0-x1) * (x2-x1))
    
    b = ((f2-f1) * (x0-x1)**2 - (f0-f1) * (x2-x1)**2) &
         / ((x0-x2) * (x0-x1) * (x2-x1))
    
    c = f1
    
    d = sqrt(b*b - 4.d0*a*c)
    
    ! Begin root finding
    do iterations = 1, max_iterations
       ! Determine if Secant or Bisection step
       if (b .ge. 0) then
          if (((x0-x1)*(b+d)+2*c)*2*c .le. 0) then
             ! OK to take a Secant step
             dx = -2.d0 * c / (b + d)
             x3 = x1 + dx
             type = '(Secant)'
          else
             ! take a bisection step instead
             dx = (x1 - x0) / 2.d0
             x3 = (x0 + x1) / 2.d0
             type = '(Bisect)'
          end if
       else
          if (((x0-x1)*(b-d)+2*c)*2*c .le. 0) then
             ! OK to take a Secant step
             dx = -2.d0 * c / (b - d)
             x3 = x1 + dx
             type = '(SA)'
          else
             ! take a bisection step instead
             dx = (x1 - x0) / 2.d0
             x3 = (x0 + x1) / 2.d0
             type = '(BS)'
          end if
       end if
       ! Compare the relative error ti the Tolerance.
       if (abs(dx/x3) .le. tolerance) then
          ! Error is BELOW tolerance, the root
          ! has been found, so exit loop.
          exit
       else
          ! The relative error is too big, loop again.
          print '(a,i3,a10,a10,f12.8)', 'Iteration:', iterations, type, 'Root: ', x3
          f3 = f(x3)
          
          ! Adjust brackets.
          if (f0 * f2 .le. 0) then
             ! The root is in left subinterval...
             x1 = x3
             f1 = f3
          else
             ! The root is in right subinterval
             x0 = x3
             f0 = f3
          end if
          ! Recalculate terms
          x2 = (x0 + x1) / 2.d0  
          f2 = f(x2)
          
          a = ((f0-f1) * (x2-x1) - (f2-f1) * (x0-x1)) &
               / ((x0-x2) * (x0-x1) * (x2-x1))
          
          b = ((f2-f1) * (x0-x1)**2 - (f0-f1) * (x2-x1)**2) &
               / ((x0-x2) * (x0-x1) * (x2-x1))
          
          c = f1
          
          d = sqrt(b*b - 4.d0*a*c)
       end if
    end do

    ! Final extimate of root
    root = x3
    
    ! normal exit
    return
    ! error trap
100 print '(a)', "Root not bracketed."
    return
  end subroutine hybridSA

end module roots
