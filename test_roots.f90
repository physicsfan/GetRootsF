PROGRAM troots
  use test_functions
  use roots
 ! This program tests the exhustive search root finder.
  implicit none

  real(8) :: xleft, xright, step_size
  real(8), parameter :: tol = 5.d-08
  integer :: iters, max_iters = 30

  ! Prompt user for x0 value
  write(*,'(a)', advance = 'no') "Enter a starting xL, xR: "
  read(*,*) xleft, xright
  write(*,'(a)', advance = 'no') "Enter a step-size: "
  read(*,*) step_size
  ! Run root_search algorithm with derivatives
  call root_searchD(xleft, xright, tol, P8, dP8, max_iters, step_size)
  ! Run root_search algorithm without exhaustive search
  WRITE(*,'(a)', advance = 'no') "Enter a starting xL, xR: "
  READ(*,*) xleft, xright
  call root_search(xleft, xright, tol, P8, max_iters)

END PROGRAM troots
