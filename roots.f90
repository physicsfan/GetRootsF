MODULE roots
  USE iso_fortran_env, only: real32
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
    REAL(real32), INTENT(inout) :: lbound, rbound
    REAL(real32), INTENT(in) :: tol
    
    INTERFACE
       PURE FUNCTION f(x) RESULT(answer)
         IMPORT
         REAL(real32), INTENT(in) :: x
         REAL(real32) :: answer
       END FUNCTION f
       
       PURE FUNCTION df(x) RESULT(answer)
         IMPORT
         REAL(real32), INTENT(in) :: x
         REAL(real32) :: answer
       END FUNCTION df
    END INTERFACE
    
    ! optional arguments
    INTEGER, INTENT(in), OPTIONAL :: max_iters
    REAL(real32), INTENT(in), OPTIONAL :: step
    
    ! locals
    INTEGER :: i, iters, max_iterations, root_number, total_steps
    REAL(real32) :: lbrac, rbrac, root, step_size
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
    IMPLICIT NONE
    
    ! arguments and interfaces
    REAL(real32), INTENT(inout) :: lbound, rbound, root
    REAL(real32), INTENT(in) :: tolerance
    INTEGER, INTENT(out) :: iterations
    integer, intent(in) :: max_iterations

    INTERFACE
       PURE FUNCTION f(x) RESULT(answer)
         IMPORT
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function f
       
       PURE FUNCTION df(x) RESULT(answer)
         IMPORT
         real(real32), intent(in) :: x
         real(real32) :: answer
       end function df
    end interface
    
    ! locals
    REAL(real32) :: delta, dfxmid, error, fxmid, fx0, fx1, x0, x1, xmid
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
       IF (ABS(delta/xmid) .LE. tolerance) THEN
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

       END IF

       PRINT '(a,i2,a8,a10,f12.8)', 'Iteration:', iterations, TYPE, 'Root: ', xmid

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
    REAL(real32), INTENT(in) :: lbound, rbound, tol

    INTERFACE
       PURE FUNCTION f(x) RESULT(answer)
         import
         REAL(real32), INTENT(in) :: x
         REAL(real32) :: answer
       END FUNCTION f
    END INTERFACE

    ! optional arguments
    INTEGER, INTENT(in), OPTIONAL :: max_iterations
    REAL(real32), INTENT(in), OPTIONAL :: step

    ! locals
    INTEGER :: i, iters, max_iters, root_number, total_steps
    REAL(real32) :: lbrac, rbrac, root, step_size
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

  IMPLICIT NONE
  
  ! arguments and interfaces
  REAL(real32), INTENT(inout) :: lbound, rbound, root
  REAL(real32), INTENT(in) :: tolerance
  INTEGER, INTENT(out) :: iterations
  INTEGER, INTENT(in) :: max_iterations
  
  INTERFACE
     PURE FUNCTION f(x) RESULT(answer)
       import
       REAL(real32), INTENT(in) :: x
       REAL(real32) :: answer
     END FUNCTION f
  END INTERFACE
  
  ! locals
  REAL(real32) :: a, b, c, d, dx, error
  REAL(real32) :: x0, x1, x2, x3, f0, f1, f2, f3
  CHARACTER(len=8) :: TYPE
  
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
  
  d = SQRT(b*b - 4.d0*a*c)
  
  ! Begin root finding
  DO iterations = 1, max_iterations
     
     ! Determine if Secant or Bisection step
     IF (b .GE. 0) THEN
        IF (((x0-x1)*(b+d)+2*c)*2*c .LE. 0) THEN
           ! OK to take a Secant step
           dx = -2.d0 * c / (b + d)
           x3 = x1 + dx
           TYPE = '(Secant)'
        ELSE
           ! take a bisection step instead
           dx = (x1 - x0) / 2.d0
           x3 = (x0 + x1) / 2.d0
           TYPE = '(Bisect)'
        END IF
     ELSE
        IF (((x0-x1)*(b-d)+2*c)*2*c .LE. 0) THEN
           ! OK to take a Secant step
           dx = -2.d0 * c / (b - d)
           x3 = x1 + dx
           TYPE = '(SA)'
        ELSE
           ! take a bisection step instead
           dx = (x1 - x0) / 2.d0
           x3 = (x0 + x1) / 2.d0
           TYPE = '(BS)'
        END IF
     END IF
     
     ! Compare the relative error ti the Tolerance.
     IF (ABS(dx/x3) .LE. tolerance) THEN
        ! Error is BELOW tolerance, the root
        ! has been found, so exit loop.
        EXIT
     ELSE
        ! The relative error is too big, loop again.
        PRINT '(a,i3,a10,a10,f12.8)', 'Iteration:', iterations, TYPE, 'Root: ', x3
        f3 = f(x3)
        
        ! Adjust brackets.
        IF (f0 * f2 .LE. 0) THEN
           ! The root is in left subinterval...
           x1 = x3
           f1 = f3
        ELSE
           ! The root is in right subinterval
           x0 = x3
           f0 = f3
        END IF
        
        ! Recalculate terms
        x2 = (x0 + x1) / 2.d0  
        f2 = f(x2)
        
        a = ((f0-f1) * (x2-x1) - (f2-f1) * (x0-x1)) &
             / ((x0-x2) * (x0-x1) * (x2-x1))
        
        b = ((f2-f1) * (x0-x1)**2 - (f0-f1) * (x2-x1)**2) &
             / ((x0-x2) * (x0-x1) * (x2-x1))
        
        c = f1
        
        d = SQRT(b*b - 4.d0*a*c)
     END IF
  END DO
  ! Final extimate of root
  root = x3
  
  ! normal exit
  RETURN
  
  ! error trap
100 PRINT '(a)', "Root not bracketed."
  RETURN
  
END SUBROUTINE hybridSA

  
END MODULE roots
