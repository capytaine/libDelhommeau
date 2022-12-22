! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/libDelhommeau>

MODULE CONSTANTS

  USE FLOATING_POINT_PRECISION, ONLY: PRE

  IMPLICIT NONE

  REAL(KIND=PRE), PARAMETER :: ZERO = 0
  REAL(KIND=PRE), PARAMETER :: ONE = 1

  REAL(KIND=PRE), PARAMETER :: PI = 3.141592653588979
  COMPLEX(KIND=PRE), PARAMETER :: II = (0, 1)         ! Imaginary unit

END MODULE CONSTANTS
