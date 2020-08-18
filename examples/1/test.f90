program test

  use matrices, only: build_matrices
  use initialize_green_wave, only: initialize_tabulated_integrals
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
  integer, parameter :: tabulation_x = 328
  integer, parameter :: tabulation_z = 46
  real(kind=pre), dimension(tabulation_x)                      :: x
  real(kind=pre), dimension(tabulation_z)                      :: z
  real(kind=pre), dimension(tabulation_x, tabulation_z, 2, 2)  :: tabulation

  ! Prony decomposition for the finite depth Green function
  integer, parameter    :: nexp = 31
  real(kind=pre), dimension(nexp) :: ambda, ar

  ! The interaction matrices to be computed
  complex(kind=pre), dimension(nb_faces, nb_faces) :: S, K

  wavenumber = 1.0
  depth = 0.0  ! means infinite depth

  vertices = reshape([  &
                      0.0, 0.0, -1.0,  &
                      1.0, 0.0, -1.0,  &
                      1.0, 1.0, -1.0,  &
                      0.0, 1.0, -1.0,  &
                      1.0, 0.0, -1.0,  &
                      2.0, 0.0, -1.0,  &
                      2.0, 1.0, -1.0,  &
                      1.0, 1.0, -1.0   &
                      ], shape(vertices))
  faces = reshape([1, 2, 3, 4, 5, 6, 7, 8], shape(faces))
                    
  face_center = reshape([0.5, 0.5, -1.0,  &
                         1.5, 0.5, -1.0], &
                        shape(face_center))
  face_normal = reshape([0.0, 0.0, 1.0,  &
                         0.0, 0.0, 1.0], &
                        shape(face_normal))

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
    [0d0, -0d0, 1d0],                                            &
    x, z, tabulation,                                            &
    nexp, ambda, ar,                                             &
    .false.,                                                      &
    S, K)

  print*, S
  print*, K

end program test
