MODULE roots
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
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: root_searchD, hybridD          ! routines if the derivative is available
  PUBLIC :: root_search, hybridSA          ! routines if the derivative is unavailable

  
CONTAINS

  
  SUBROUTINE root_searchD(lbound, rbound, tol, f, df, max_iters, step)
    ! Uses exhausive searching for bracketing all the roots within a user-defined range,
    ! and then calls the hybridD routine to find the root in each bracket.  The step-size
    ! for the search is an optional user input.  If none is supplied, 0.1 is used. The
    ! maximum number of iterations is also optional. If none is supplied, 30 iterations
    ! are allowed.
    IMPLICIT NONE

    ! arguments and interfaces
    REAL(8), INTENT(inout) :: lbound, rbound
    REAL(8), INTENT(in) :: tol
    INTERFACE
       PURE FUNCTION f(x) RESULT(answer)
         REAL(8), INTENT(in) :: x
         REAL(8) :: answer
       END FUNCTION f
       
       PURE FUNCTION df(x) RESULT(answer)
         REAL(8), INTENT(in) :: x
         REAL(8) :: answer
       END FUNCTION df
    END INTERFACE
    ! optional arguments
    INTEGER, INTENT(in), OPTIONAL :: max_iters
    REAL(8), INTENT(in), OPTIONAL :: step
    ! locals
    INTEGER :: i, iters, max_iterations, root_number, total_steps
    REAL(8) :: lbrac, rbrac, root, step_size
    LOGICAL :: root_present = .FALSE.

    ! begin
    max_iterations = 30; step_size = 0.1d0
    root_number = 0
    
    ! Check for optional arguments
    IF(PRESENT(max_iters)) max_iterations = max_iters
    IF(PRESENT(step)) step_size = step

    ! Perform exhaustive search
    ! ensure steps are not too small
    IF(step_size <= tol) GOTO 100
    total_steps = INT((rbound-lbound)/step_size)
    PRINT *, 'Total steps:', total_steps
    DO i = 1, total_steps-1
       ! Select search interval for root i
       lbrac = lbound+step_size*i
       rbrac = lbound+step_size*(i+1)
       ! Verify that the root is bracketed
       IF (f(lbrac)*f(rbrac) < 0.d0) THEN
          root_present = .TRUE.
          root_number = root_number + 1
          ! Find the root within interval
          CALL hybridD(lbrac, rbrac, root, tol, iters, max_iters, f, df)
          ! Print results
          IF (iters .GE. max_iters) THEN
             PRINT '(a)', "Maximum iterations reached."
             PRINT '(a,i2,a,f11.8,a,i3,a)', 'Last approximation to the root', root_number, &
                  ' found at x =', root, ' after', max_iters,' iterations.'
          ELSE
             PRINT '(a,i3,a,f11.8,a,i3,a)', 'Root', root_number, ' found at x =', root, &
                  ' in', iters,' iterations.'
          END IF
       END IF
    END DO
       
    IF (root_present .EQV. .FALSE.) GOTO 200

    ! normal exit
    RETURN

    ! error trap
100 PRINT '(a,f11.8,a)', 'The specified step-size dx =', step_size, ' is too small.'
    RETURN
200 PRINT '(a)', "No roots are bracketed."
    RETURN
    
  END SUBROUTINE root_searchD

  
  SUBROUTINE hybridD(lbound, rbound, root, tolerance, iterations, max_iterations, &
       f, df)
    ! A hybrid Bisection/Newton-Raphson method for finding the root of a
    ! function F, with derivative Fprime. The root is known to be between "lbound"
    ! and "rbound", and the result is returned in "root". If the next NR guess
    ! is within the known bounds, the step is accepted; otherwise, a bisection
    ! step is taken.
    !
    ! The root is initially bracketed between X0 and X1:
    !
    !    x :        X0          Xmid           X1
    !  f(x):        fX0         fXmid          fX1
    ! f'(x):          ---       DerfXmid            ---
    !
    ! To Do: All errors should be handled by Error Handler
    ! arguments and interfaces
    REAL(8), INTENT(inout) :: lbound, rbound, root
    real(8), intent(in) :: tolerance
    integer, intent(out) :: iterations
    integer, intent(in) :: max_iterations
    interface
       pure function f(x) result(answer)
         real(8), intent(in) :: x
         real(8) :: answer
       end function f
       
       pure function df(x) result(answer)
         real(8), intent(in) :: x
         real(8) :: answer
       end function df
    end interface
    
    ! locals
    REAL(8) :: delta, dfxmid, error, fxmid, fx0, fx1, x0, x1, xmid
    character(len=6) :: type
    
    ! Verify that the root is bracketed
    IF (f(lbound) * f(rbound) .GT. 0.d0) GOTO 100

    ! Initialization
    x0 = lbound; x1 = rbound
    fx0 = f(x0)
    fx1 = f(x1)
    IF (ABS(fx0) .LE. ABS(fx1)) THEN
       xmid = x0
       fxmid = fx0
    ELSE
       xmid = x1
       fxmid = fx1
    END IF
    dfxmid = df(xmid)
    
    ! Begin root finding
    DO iterations = 1, max_iterations
       
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
    END DO

    ! normal exit
    root = xmid
    RETURN
    
    ! error trap
100 PRINT '(a)', "Root not bracketed."
    RETURN
  END SUBROUTINE hybridD


  SUBROUTINE root_search(lbound, rbound, tol, f, max_iterations, step)
    ! Uses exhausive searching for bracketing all the roots within a user-defined range,
    ! and then calls the "hybridD" routine to find the root in each bracket, where no
    ! derivative is used.  The step-size for the search is an optional user input.  If
    ! none is supplied, 0.1 is used. The maximum number of iterations is also optional.
    ! If none is supplied, 30 is used.
    IMPLICIT NONE
    ! arguments and interfaces
    REAL(8), INTENT(in) :: lbound, rbound, tol
    INTERFACE
       PURE FUNCTION f(x) RESULT(answer)
         REAL(8), INTENT(in) :: x
         REAL(8) :: answer
       END FUNCTION f
    END INTERFACE
    ! optional arguments
    INTEGER, INTENT(in), OPTIONAL :: max_iterations
    REAL(8), INTENT(in), OPTIONAL :: step
    ! locals
    INTEGER :: i, iters, max_iters, root_number, total_steps
    REAL(8) :: lbrac, rbrac, root, step_size
    LOGICAL :: root_present = .FALSE.

    ! begin
    max_iters = 30; step_size = 0.1d0
    root_number = 0
    
    ! Check for optional arguments
    IF(PRESENT(max_iterations)) max_iters = max_iterations
    IF(PRESENT(step)) step_size = step

    ! Perform exhaustive search
    ! ensure steps are not too small
    IF(step_size <= tol) GOTO 100
    total_steps = INT((rbound-lbound)/step_size)
    PRINT *, 'Total steps:', total_steps

    DO i = 1, total_steps-1
       ! Select search interval for root i
       lbrac = lbound+step_size*i
       rbrac = lbound+step_size*(i+1)
       ! Verify that the root is bracketed
       IF (f(lbrac)*f(rbrac) < 0.d0) THEN
          root_present = .TRUE.
          root_number = root_number + 1
          ! Find the root within interval
          CALL hybridSA(lbrac, rbrac, root, tol, iters, max_iters, f)
          ! Print results
          IF (iters .GE. max_iters) THEN
             PRINT '(a)', "Maximum iterations reached."
             PRINT '(a,i2,a,f11.8,a,i3,a)', 'Last approximation to the root', root_number, &
                  ' found at x =', root, ' after', max_iters,' iterations.'
          ELSE
             PRINT '(a,i3,a,f11.8,a,i3,a)', 'Root', root_number, ' found at x =', root, &
                  ' in', iters,' iterations.'
          END IF
       END IF
    END DO
       
    IF (root_present .EQV. .FALSE.) GOTO 200

    ! normal exit
    RETURN

    ! error trap
100 PRINT '(a,f11.8,a)', 'The specified step-size dx =', step_size, ' is too small.'
    RETURN
200 PRINT '(a)', "No roots are bracketed."
    RETURN
        
  END SUBROUTINE root_search


  
SUBROUTINE hybridSA(lbound, rbound, root, tolerance, iterations, max_iterations, f)
  ! A hybrid Bisection/Secant method for finding the root of a
  ! function F when the derivative of F is not available. The Secant method
  ! uses acceleration to speed converdenge. The root is known
  ! to be between "lbound" and "rbound", and the result is returned in "root".
  ! If the next SA guess is within the known bounds, the step is accepted;
  ! otherwise, a bisection step is taken.
  !
  ! arguments and interfaces
  REAL(8), INTENT(inout) :: lbound, rbound, root
  real(8), intent(in) :: tolerance
  integer, intent(out) :: iterations
  integer, intent(in) :: max_iterations
  interface
     pure function f(x) result(answer)
       real(8), intent(in) :: x
       real(8) :: answer
     end function f
  end interface

  ! locals
  real(8) :: a, b, c, d, dx, error
  real(8) :: x0, x1, x2, x3, f0, f1, f2, f3
  character(len=8) :: type
  
  ! Verify that the root is bracketed
  IF (f(lbound) * f(rbound) .GT. 0.d0) GOTO 100

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
    RETURN

    ! error trap
100 PRINT '(a)', "Root not bracketed."
    return
        
END SUBROUTINE hybridSA

  
end module roots
