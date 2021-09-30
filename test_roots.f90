PROGRAM troots
  USE iso_fortran_env, ONLY: real32
  USE test_functions
  USE roots
 ! This program tests the exhustive search root finder.
  IMPLICIT NONE

  REAL(real32) :: xleft, xright, step_size
  REAL(real32), PARAMETER :: tol = 5.d-08
  INTEGER :: iters, max_iters = 30
  
  ! Prompt user for x0 value
  WRITE(*,'(a)', advance = 'no') "Enter a starting xL, xR: "
  READ(*,*) xleft, xright
  WRITE(*,'(a)', advance = 'no') "Enter a step-size: "
  READ(*,*) step_size
  ! Run root_search algorithm with derivatives
  CALL root_searchD(xleft, xright, tol, P8, dP8, max_iters, step_size)
  ! Run root_search algorithm without exhaustive search
  WRITE(*,'(a)', advance = 'no') "Enter a starting xL, xR: "
  READ(*,*) xleft, xright
  CALL root_search(xleft, xright, tol, P8, max_iters)
  
END PROGRAM troots
