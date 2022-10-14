! Copyright (C) 2017-2022 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
module green_rankine

  use constants
  use face_type_module

  implicit none

  public
contains

  pure subroutine integral_of_rankine_source(x, face, n, g, dgdn)

    real(kind=pre), dimension(3), intent(in) :: x
    type(face_type),              intent(in) :: face
    real(kind=pre), dimension(3), intent(in) :: n

    real(kind=pre), intent(out) :: g
    real(kind=pre), intent(out) :: dgdn

    real(kind=pre) :: R

    R = norm2(x(1:3) - face%center(1:3)) ! Distance from center of face to M

    if (R < 7*face%radius) then
      call exact_integral_of_rankine_source(x, face, n, g, dgdn)
    else
      call one_point_integration_of_rankine_source(x, face, n, g, dgdn)
    endif
  end subroutine

  ! ====================================

  pure subroutine one_point_integration_of_rankine_source &
      (x, face, n, g, dgdn)
    ! Estimate ∫∫_Γ 1/||x-ξ|| dξ ~ area(Γ)/||x - ξ_0||
    ! and its gradient with respect to x.

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: x
    type(face_type),              intent(in) :: face
    real(kind=pre), dimension(3), intent(in) :: n

    ! Outputs
    real(kind=pre), intent(out) :: g
    real(kind=pre), intent(out) :: dgdn

    ! Local variables
    real(kind=pre) :: R

    R    = norm2(x(1:3) - face%center(1:3)) ! Distance from center of face to M
    g    = face%area/R
    dgdn = g*dot_product(face%center(1:3) - x, n)/R**2

  end subroutine one_point_integration_of_rankine_source

  ! ====================================

  pure subroutine exact_integral_of_rankine_source &
      (x, face, n, g, dgdn)
    ! Compute the integral ∫∫_Γ 1/||x-ξ|| dξ over a face Γ
    ! and its derivative with respect to x.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: x
    type(face_type),              intent(in) :: face
    real(kind=pre), dimension(3), intent(in) :: n

    ! Outputs
    real(kind=pre), intent(out) :: g
    real(kind=pre), intent(out) :: dgdn

    ! The index of the following node when going around a face.
    integer, dimension(4), parameter :: next_node = [2, 3, 4, 1]

    ! Local variables
    INTEGER                         :: L
    REAL(KIND=PRE)                  :: GZ, DK, GY
    REAL(KIND=PRE), DIMENSION(4)    :: RR
    REAL(KIND=PRE), DIMENSION(3, 4) :: DRX
    REAL(KIND=PRE)                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    REAL(KIND=PRE), DIMENSION(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    GZ = DOT_PRODUCT(x(1:3) - face%center(1:3), face%normal(1:3)) ! Called Z in [Del]

    DO CONCURRENT (L = 1:4)
      RR(L) = NORM2(x(1:3) - face%vertices(L, 1:3))       ! Distance from vertices of Face to M.
      DRX(:, L) = (x(1:3) - face%vertices(L, 1:3))/RR(L)  ! Normed vector from vertices of Face to M.
    END DO

    DO L = 1, 4
      DK = NORM2(face%vertices(NEXT_NODE(L), :) - face%vertices(L, :))    ! Distance between two consecutive points, called d_k in [Del]
      IF (DK >= REAL(1e-3, PRE)*face%radius) THEN
        PJ(:) = (face%vertices(NEXT_NODE(L), :) - face%vertices(L, :))/DK ! Normed vector from one corner to the next
        ! The following GYX(1:3) are called (a,b,c) in [Del]
        GYX(1) = face%normal(2)*PJ(3) - face%normal(3)*PJ(2)
        GYX(2) = face%normal(3)*PJ(1) - face%normal(1)*PJ(3)
        GYX(3) = face%normal(1)*PJ(2) - face%normal(2)*PJ(1)
        GY = DOT_PRODUCT(x - face%vertices(L, :), GYX)                                    ! Called Y_k in  [Del]

        ANT = 2*GY*DK                                                                  ! Called N^t_k in [Del]
        DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2*ABS(GZ)*(RR(NEXT_NODE(L))+RR(L)) ! Called D^t_k in [Del]
        ANL = RR(NEXT_NODE(L)) + RR(L) + DK                                            ! Called N^l_k in [Del]
        DNL = RR(NEXT_NODE(L)) + RR(L) - DK                                            ! Called D^l_k in [Del]
        ALDEN = LOG(ANL/DNL)

        IF (ABS(GZ) >= REAL(1e-4, PRE)*face%radius) THEN
          AT = ATAN(ANT/DNT)
        ELSE
          AT = 0.
        ENDIF

        ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

        ANTX(:) = 2*DK*GYX(:)                                         ! Called N^t_k_{x,y,z} in [Del]
        DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + ABS(GZ))*ANLX(:) &
          + 2*SIGN(ONE, GZ)*(RR(NEXT_NODE(L)) + RR(L))*face%normal(:) ! Called D^t_k_{x,y,z} in [Del]

        IF (ABS(GY) < 1e-5) THEN
          ! Edge case where the singularity is on the boundary of the face (GY = 0, ALDEN = infty).
          ! This case seems to only occur when computating the free surface elevation,
          ! so no fix has been implemented for VS0, which is not needed then.
          G = G - 2*AT*ABS(GZ)
        ELSE
          ! General case
          G = G + GY*ALDEN - 2*AT*ABS(GZ)
        END IF

        dGdn = dGdn + dot_product(n, ALDEN*GYX(:)     &
          - 2*SIGN(ONE, GZ)*AT*face%normal(:)   &
          + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:) &
          - 2*ABS(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT))
      END IF
    END DO

  end subroutine exact_integral_of_rankine_source

END MODULE GREEN_RANKINE
