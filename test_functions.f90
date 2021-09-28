module test_functions
  ! This module contains non-linear test functions for testing the
  ! root finding routines


contains
  
  pure function f1(x) result(answer)
    ! This is a nonlinear function whose root
    ! cannot be found analytically.
    real(8), intent(in) :: x
    real(8) :: answer
    answer = x - exp(1/x)
  end function f1

!===============================================================

  pure function df1(x) result(answer)
    ! This function is the derivative of "f1".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 1.d0 + exp(1/x) / x**2
  end function df1
  
!===============================================================  
  
  pure function f2(x) result(answer)
    ! This is an axample of a nonlinear function whose root
    ! cannot be found analytically.
    real(8), intent(in) :: x
    real(8) :: answer
    answer = cos(x) - x
  end function f2

!===============================================================

  pure function df2(x) result(answer)
    ! This function is the derivative of "f2".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = -sin(x) - 1.d0
  end function df2

!===============================================================  
  
  pure function f3(x) result(answer)
    real(8), intent(in) :: x
    real(8) :: answer
    answer = x**3 - 169.d0
  end function f3

!===============================================================

  pure function df3(x) result(answer)
    ! This function is the derivative of "f3".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 3.d0 * x**2
  end function df3

!===============================================================  
  
  pure function f4(x) result(answer)
    real(8), intent(in) :: x
    real(8) :: answer
    answer = x**2 + 1.d0
  end function f4

!===============================================================

  pure function df4(x) result(answer)
    ! This function is the derivative of "f3".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 2.d0 * x
  end function df4
  
!===============================================================

  pure function P8(x) result(answer)
    ! Legendre polynomial
    real(8), intent(in) :: x
    real(8) :: answer
    answer = (6435.d0*x**8 - 12012.d0*x**6 + 6930.d0*x**4 &
         - 1260.d0*x**2 + 35.d0) / 128.d0
  end function P8

!===============================================================

  pure function dP8(x) result(answer)
    ! Derivative of Legendre polynomial
    real(8), intent(in) :: x
    real(8) :: answer
    answer = (51480.d0*x**7 - 72072.d0*x**5 + 27720.d0*x**3 &
         - 2520.d0*x) / 128.d0
  end function dP8

!===============================================================  
  
  pure function f5(x) result(answer)
    real(8), intent(in) :: x
    real(8) :: answer
    answer = x**2 - 2.d0*x - 2.d0 
  end function f5

!===============================================================

  pure function df5(x) result(answer)
    ! This function is the derivative of "f5".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 2.d0 * x - 2.d0
  end function df5
  
!===============================================================

  pure function f6(x) result(answer)
    ! This is an axample of a nonlinear function whose root
    ! cannot be found analytically.
    real(8), intent(in) :: x
    real(8) :: answer
    answer = x - cos(x)
  end function f6

!===============================================================

  pure function df6(x) result(answer)
    ! This function is the derivative of "f2".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 1.d0 - sin(x)
  end function df6

!===============================================================

  pure function f7(x) result(answer)
    ! This is an axample of a nonlinear function whose root
    ! cannot be found analytically.
    real(8), intent(in) :: x
    real(8) :: answer
    answer = sin(x) - x/2.d0
  end function f7

!===============================================================

  pure function df7(x) result(answer)
    ! This function is the derivative of "f2".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = cos(x) - 1.d0/2.d0
  end function df7

!===============================================================

  pure function f8(x) result(answer)
    ! This is an axample of a nonlinear function whose root
    ! cannot be found analytically.
    real(8), intent(in) :: x
    real(8) :: answer
    answer = cos(x) - x * sin(x)
  end function f8

!===============================================================

  pure function df8(x) result(answer)
    ! This function is the derivative of "f2".
    real(8), intent(in) :: x
    real(8) :: answer
    answer = 2.d0*sin(x) + x * cos(x)
  end function df8
  
end module test_functions
