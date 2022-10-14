module face_type_module

  use constants

  implicit none

  public

  type face_type
    real(kind=pre), dimension(4, 3) :: vertices
    real(kind=pre), dimension(3)    :: center
    real(kind=pre), dimension(3)    :: normal
    real(kind=pre)                  :: area
    real(kind=pre)                  :: radius
    real(kind=pre), dimension(:), allocatable    :: quad_weights
    real(kind=pre), dimension(:, :), allocatable :: quad_points

  end type face_type
end module face_type_module
