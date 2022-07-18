program benchmarck_omp

use ieee_arithmetic

use matrices, only: build_matrices
use delhommeau_integrals, only: default_r_spacing, default_z_spacing, construct_tabulation
use constants, only: pi, pre  ! Floating point precision

implicit none

integer(kind=8) :: starting_time, final_time, clock_rate

integer, parameter :: nb_faces = 100
integer, parameter :: nb_vertices = 4*nb_faces
integer, parameter :: nb_quadrature_points = 1

real(kind=pre) :: wavenumber, depth

! Geometry of the mesh
real(kind=pre), dimension(nb_vertices, 3) :: vertices
integer, dimension(nb_faces, 4) :: faces
real(kind=pre), dimension(nb_faces, 3) :: face_center
real(kind=pre), dimension(nb_faces, 3) :: face_normal
real(kind=pre), dimension(nb_faces) :: face_area
real(kind=pre), dimension(nb_faces) :: face_radius
real(kind=pre), dimension(nb_faces, nb_quadrature_points, 3) :: quadrature_points
real(kind=pre), dimension(nb_faces, nb_quadrature_points) :: quadrature_weights

! Tabulation of the integrals used in the Green function
integer, parameter :: tabulation_nr = 328
integer, parameter :: tabulation_nz = 46
real(kind=pre), dimension(tabulation_nr)                       :: tabulated_r
real(kind=pre), dimension(tabulation_nz)                       :: tabulated_z
real(kind=pre), dimension(tabulation_nr, tabulation_nz, 2, 2)  :: tabulated_integrals

! Prony decomposition for the finite depth Green function
integer, parameter    :: nexp = 31
real(kind=pre), dimension(nexp) :: ambda, ar

! The interaction matrices to be computed
complex(kind=pre), dimension(nb_faces, nb_faces) :: S, K

tabulated_r(:) = default_r_spacing(tabulation_nr)
tabulated_z(:) = default_z_spacing(tabulation_nz)
tabulated_integrals(:, :, :, :) = construct_tabulation(tabulated_r, tabulated_z, 251)

wavenumber = 1.0
depth = ieee_value(depth, ieee_positive_inf)

call random_panels(nb_faces, vertices, faces, face_center, face_normal, face_area, face_radius)

! ! For debugging
! open (unit=4, file='vertices.dat', form='formatted')
! do i = 1, 4*nb_faces
!    write (4, *) vertices(i, :)
! end do

quadrature_points = reshape(face_center, shape(quadrature_points))
quadrature_weights = reshape(face_area, shape(quadrature_weights))

call system_clock(count_rate=clock_rate)
call system_clock(starting_time)

call build_matrices(                                           &
  nb_faces, face_center, face_normal,                          &
  nb_vertices, nb_faces, vertices, faces,                      &
  face_center, face_normal, face_area, face_radius,            &
  nb_quadrature_points, quadrature_points, quadrature_weights, &
  wavenumber, depth,                                           &
  [1d0, -1d0, 1d0],                                            &
  tabulated_r, tabulated_z, tabulated_integrals,               &
  nexp, ambda, ar,                                             &
  .true.,                                                      &
  S, K)

call system_clock(final_time)

print*, "Elapsed time:", real(final_time - starting_time)/clock_rate


contains

  subroutine random_panels(nb_faces, vertices, faces, face_center, face_normal, face_area, face_radius)
    integer, intent(in) :: nb_faces

    integer, dimension(nb_faces, 4), intent(out) :: faces
    real(kind=pre), dimension(4*nb_faces, 3), intent(out) :: vertices
    real(kind=pre), dimension(nb_faces, 3), intent(out) :: face_center
    real(kind=pre), dimension(nb_faces, 3), intent(out) :: face_normal
    real(kind=pre), dimension(nb_faces), intent(out) :: face_area
    real(kind=pre), dimension(nb_faces), intent(out) :: face_radius

    real(kind=pre), dimension(3, 2) :: vertex_shifts
    integer :: i

    faces = transpose(reshape([(i, i=1,(nb_faces*4),1)], [4, nb_faces]))

    call random_number(face_area)
    face_area = 1.0 + face_area

    face_radius = face_area * sqrt(2.0)/2.0

    call random_number(face_center)
    face_center(:, 1) = 20*(face_center(:, 1) - 0.5)
    face_center(:, 2) = 20*(face_center(:, 2) - 0.5)
    face_center(:, 3) = -10*face_center(:, 3) - face_radius(:)

    call random_number(face_normal)
    face_normal(:, :) = face_normal(:, :) - 0.5
    do i = 1, nb_faces
      face_normal(i, :) = face_normal(i, :)/norm2(face_normal(i, :))
    enddo

    do i = 1, nb_faces
      vertex_shifts = face_radius(i) * two_orthogonal_vector(face_normal(i, :))
      vertices(4*(i-1)+1, :) = face_center(i, :) + vertex_shifts(:, 1)
      vertices(4*(i-1)+2, :) = face_center(i, :) + vertex_shifts(:, 2)
      vertices(4*(i-1)+3, :) = face_center(i, :) - vertex_shifts(:, 1)
      vertices(4*(i-1)+4, :) = face_center(i, :) - vertex_shifts(:, 2)
    enddo
  end subroutine

  pure function two_orthogonal_vector(n) result(vecs)
    ! Given a normal vector `n`, returns two oter vectors
    ! such that the three of them is an orthonormal basis
    real(kind=pre), dimension(3), intent(in) :: n
    real(kind=pre), dimension(3, 2) :: vecs

    vecs(:, 1) = [n(2), -n(1), 0d0]
    if (norm2(vecs(:, 1)) < abs(1e-5)) then
      vecs(:, 1) = [n(3), 0d0, -n(1)]
    endif

    ! Cross product
    vecs(1, 2) = n(2) * vecs(3, 1) - n(3) * vecs(2, 1)
    vecs(2, 2) = n(3) * vecs(1, 1) - n(1) * vecs(3, 1)
    vecs(3, 2) = n(1) * vecs(2, 1) - n(2) * vecs(1, 1)

    vecs(:, 1) = vecs(:, 1)/norm2(vecs(:, 1))
    vecs(:, 2) = vecs(:, 2)/norm2(vecs(:, 2))
  end function

end program benchmarck_omp