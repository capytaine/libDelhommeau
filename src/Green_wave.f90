! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/libDelhommeau>
MODULE GREEN_WAVE

  USE ieee_arithmetic
  USE CONSTANTS
  USE FACE_TYPE_MODULE
  USE DELHOMMEAU_INTEGRALS

  IMPLICIT NONE

CONTAINS

  subroutine integral_of_wave_part &
      (x, face, n, k, depth, &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      nexp, ambda, ar, &
      g, dgdn_sym, dgdn_antisym)

    real(kind=pre), dimension(3), intent(in)  :: x, n
    type(face_type),              intent(in)  :: face
    real(kind=pre),               intent(in)  :: k, depth
    complex(kind=pre),            intent(out) :: g, dgdn_sym, dgdn_antisym

    real(kind=pre), dimension(:),          intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),          intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :, :), intent(in) :: tabulated_integrals

    integer,                          intent(in) :: nexp
    real(kind=pre), dimension(nexp),  intent(in) :: ambda, ar

    integer :: q
    complex(kind=pre) :: g_q, dgdn_sym_q, dgdn_antisym_q

    g = zero
    dgdn_sym = zero
    dgdn_antisym = zero

    do q = 1, size(face%quad_points, 1)

      if (is_infinity(depth)) then
        call wave_part_infinite_depth &
          (x, face%quad_points(q, :), n, k, &
          tabulated_r_range, tabulated_z_range, tabulated_integrals, &
          g_q, dgdn_sym_q, dgdn_antisym_q)
      else
        call wave_part_finite_depth &
          (x, face%quad_points(q, :), n, k, depth, &
          tabulated_r_range, tabulated_z_range, tabulated_integrals, &
          NEXP, AMBDA, AR, &
          g_q, dgdn_sym_q, dgdn_antisym_q)
      endif

      g = g + g_q * face%quad_weights(q)
      dgdn_sym = dgdn_sym + dgdn_sym_q * face%quad_weights(q)
      dgdn_antisym = dgdn_antisym + dgdn_antisym_q * face%quad_weights(q)
    end do

  contains

    pure logical function is_infinity(x)
      real(kind=pre), intent(in) :: x
      is_infinity = (.not. ieee_is_finite(x))
    end function
  end subroutine

  ! =====================================================================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH &
      (x, xi, n, wavenumber,          &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      g, dgdn_sym, dgdn_antisym)

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN)  :: wavenumber
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: x, xi, n

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),          INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),          INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :, :), INTENT(IN) :: tabulated_integrals

    ! Outputs
    complex(kind=pre), intent(out) :: g, dgdn_sym, dgdn_antisym

    ! Local variables
    complex(kind=pre)               :: calG
    complex(kind=pre), dimension(3) :: dcalGdx
#ifndef XIE_CORRECTION
    real(kind=pre) :: two_over_r1_cube
    real(kind=pre), dimension(3) :: reflected_xi
#endif

    ! The integrals
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      x, xi, wavenumber,                                         &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      calG, dcalGdx)
    g  = 2*wavenumber*calG
    dgdn_antisym = 2*wavenumber**2*(dcalGdx(1)*n(1) + dcalGdx(2)*n(2))
    dgdn_sym = 2*wavenumber**2*dcalGdx(3)*n(3)

#ifndef XIE_CORRECTION
    ! In the original Delhommeau method, a singularity is missing in the derivative
    reflected_xi = [xi(1), xi(2), -xi(3)]
    two_over_r1_cube = 2.0/(norm2(x-reflected_xi)**3)
    dgdn_antisym = dgdn_antisym - ((x(1)-reflected_xi(1))*n(1) + (x(2)-reflected_xi(2))*n(2))*two_over_r1_cube
    dgdn_sym = dgdn_sym - (x(3) - reflected_xi(3))*n(3)*two_over_r1_cube
#endif

  end subroutine wave_part_infinite_depth

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (X0I, X0J, n, wavenumber, depth, &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      NEXP, AMBDA, AR,              &
      g, dgdn_sym, dgdn_antisym)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J, n
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),          INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),          INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :, :), INTENT(IN) :: tabulated_integrals

    ! Prony decomposition for finite depth
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Outputs
    complex(kind=pre), intent(out) :: g, dgdn_sym, dgdn_antisym

    ! Local variables
    INTEGER                              :: KE
    REAL(KIND=PRE)                       :: AMH, AKH, A
    REAL(KIND=PRE)                       :: AQT, R
    REAL(KIND=PRE),    DIMENSION(3)      :: XI, XJ
    REAL(KIND=PRE),    DIMENSION(4)      :: FTS, PSR
    REAL(KIND=PRE),    DIMENSION(3, 4)   :: VTS
    COMPLEX(KIND=PRE)                    :: SP
    COMPLEX(KIND=PRE), DIMENSION(3)      :: VSP_SYM, VSP_ANTISYM
    COMPLEX(KIND=PRE), DIMENSION(4)      :: FS
    COMPLEX(KIND=PRE), DIMENSION(3, 4)   :: VS

    !========================================
    ! Part 1: Solve 4 infinite depth problems
    !========================================

    XI(:) = X0I(:)
    XJ(:) = X0J(:)

    ! Distance in xOy plane
    R = NORM2(XI(1:2) - XJ(1:2))

    ! 1.a First infinite depth problem
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(1), VS(:, 1))

#ifndef XIE_CORRECTION
    ! In the original Delhommeau method, the integrals are Re[ ∫(J(ζ) - 1/ζ)dθ ]/π + i Re[ ∫(e^ζ)dθ ]
    ! whereas we need Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ]
    ! So the term PSR is the difference because  Re[ ∫ 1/ζ dθ ] = - π/sqrt(r² + z²)
    !
    ! Note however, that the derivative part of Delhommeau integrals is the derivative of
    ! Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ] so no fix is needed for the derivative.
    PSR(1) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

#ifndef XIE_CORRECTION
    PSR(2) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(3), VS(:, 3))

#ifndef XIE_CORRECTION
    PSR(3) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

#ifndef XIE_CORRECTION
    PSR(4) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! Add up the results of the four problems
#ifdef XIE_CORRECTION
    SP               = SUM(FS(1:4))
#else
    SP               = SUM(FS(1:4)) - SUM(PSR(1:4))
#endif
    VSP_SYM(1:3)     = VS(1:3, 1) + VS(1:3, 4)
    VSP_ANTISYM(1:3) = VS(1:3, 2) + VS(1:3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(2*depth*(AMH**2-AKH**2+AKH))

    SP          = A*SP
    VSP_ANTISYM = A*wavenumber*VSP_ANTISYM
    VSP_SYM     = A*wavenumber*VSP_SYM

    !=====================================================
    ! Part 2: Integrate (NEXP+1)×4 terms of the form 1/MM'
    !=====================================================

    DO KE = 1, NEXP
      XI(:) = X0I(:)

      ! 2.a Shift observation point and compute integral
      XI(3) =  X0I(3) + depth*AMBDA(KE) - 2*depth
      FTS(1) = one/norm2(XI - X0J)
      VTS(:, 1) = (XI - X0J)*FTS(1)**3

      ! 2.b Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) - depth*AMBDA(KE)
      FTS(2) = one/norm2(XI - X0J)
      VTS(:, 2) = (XI - X0J)*FTS(2)**3
      VTS(3, 2) = -VTS(3, 2) ! Reflection of the output vector

      ! 2.c Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) + depth*AMBDA(KE) - 4*depth
      FTS(3) = one/norm2(XI - X0J)
      VTS(:, 3) = (XI - X0J)*FTS(3)**3
      VTS(3, 3) = -VTS(3, 3) ! Reflection of the output vector

      ! 2.d Shift observation point and compute integral
      XI(3) =  X0I(3) - depth*AMBDA(KE) + 2*depth
      FTS(4) = one/norm2(XI - X0J)
      VTS(:, 4) = (XI - X0J)*FTS(4)**3

      AQT = AR(KE)/2

      ! Add all the contributions
      SP               = SP               + AQT*SUM(FTS(1:4))
      VSP_ANTISYM(1:3) = VSP_ANTISYM(1:3) + AQT*(VTS(1:3, 1) + VTS(1:3, 4))
      VSP_SYM(1:3)     = VSP_SYM(1:3)     + AQT*(VTS(1:3, 2) + VTS(1:3, 3))

    END DO

    g = sp
    dgdn_sym = dot_product(VSP_SYM, n)
    dgdn_antisym = dot_product(VSP_ANTISYM, n)
  END SUBROUTINE

  ! ======================

  SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS                        &
      (X0I, X0J, wavenumber,                                     &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS, VS)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: FS  ! the integral
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: VS  ! its gradient

    ! Local variables
    REAL(KIND=PRE) :: r, z, r1, drdx, drdy
    REAL(KIND=PRE), dimension(2, 2) :: integrals

    r = wavenumber * NORM2(X0I(1:2) - X0J(1:2))
    z = wavenumber * (X0I(3) + X0J(3))
    r1 = hypot(r, z)

    IF (ABS(r) > 16*EPSILON(r)) THEN
      drdx = wavenumber * (X0I(1) - X0J(1))/r
      drdy = wavenumber * (X0I(2) - X0J(2))/r
    ELSE
      ! Limit when r->0 is not well defined...
      drdx = ZERO
      drdy = ZERO
    END IF

    IF (z > -1e-8) THEN
      PRINT*, "Error: Impossible to compute the wave part of the Green function due to panels on the free surface (z=0) or above."
      ERROR STOP
    ENDIF

    !=======================================================
    ! Evaluate the elementary integrals depending on z and r
    !=======================================================
    IF ((MINVAL(tabulated_z_range) < z) .AND. (r < MAXVAL(tabulated_r_range))) THEN
      ! Within the range of tabulated data
      integrals = pick_in_default_tabulation(r, z, tabulated_r_range, tabulated_z_range, tabulated_integrals)

    ELSE
      ! Asymptotic expression for distant panels
      integrals = asymptotic_approximations(MAX(r, 1e-10), z)
    ENDIF

    !================================================
    ! Add the elementary integrals to build FS and VS
    !================================================

    FS    = CMPLX(integrals(1, 2)/PI, integrals(2, 2), KIND=PRE)
    VS(1) = -drdx * CMPLX(integrals(1, 1)/PI, integrals(2, 1), KIND=PRE)
    VS(2) = -drdy * CMPLX(integrals(1, 1)/PI, integrals(2, 1), KIND=PRE)
#ifdef XIE_CORRECTION
    VS(3) = CMPLX(integrals(1, 2)/PI + ONE/r1, integrals(2, 2), KIND=PRE)
#else
    VS(3) = CMPLX(integrals(1, 2)/PI, integrals(2, 2), KIND=PRE)
#endif

    RETURN
  END SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS

  ! =========================

END MODULE GREEN_WAVE
