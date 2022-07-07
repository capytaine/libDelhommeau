program test

  use ieee_arithmetic

  use matrices, only: build_matrices
  use delhommeau_integrals, only: default_r_spacing, default_z_spacing, construct_tabulation
  use constants, only: pre  ! Floating point precision

  implicit none

  integer, parameter :: nb_faces = 2
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

  vertices = reshape([  &
    0.0,  1.0,  1.0,  0.0,  1.0,  2.0,  2.0,  1.0,  &
    0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  1.0,  1.0,  &
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0  &
    ], shape(vertices))

  faces = reshape([  &
    1, 5, &
    2, 6, &
    3, 7, &
    4, 8  &
    ], shape(faces))

  face_center = reshape([  &
    0.5,  1.5,  &
    0.5,  0.5,  &
    -1.0, -1.0  &
    ], shape(face_center))

  face_normal = reshape([  &
    0.0, 0.0,  &
    0.0, 0.0,  &
    1.0, 1.0   &
    ], shape(face_normal))

  face_area = [1.0, 1.0]
  face_radius = [0.71, 0.71]

  quadrature_points = reshape(face_center, shape(quadrature_points))
  quadrature_weights = reshape(face_area, shape(quadrature_weights))

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

  print*, S
  print*, K

end program test
