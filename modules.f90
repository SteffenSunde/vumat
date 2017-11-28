module utilities
  implicit none
  public dp, sp, qp, timesTwo
  integer, parameter :: dp=kind(1.d0)
  integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  contains


    function timesTwo(A) result(Adouble)
        ! Returns the matrix A scaled by 2.0
        real(dp), intent(in) :: A(3,3)
        real(dp) :: Adouble(3,3)
        Adouble = A * 2.0
    end function timesTwo

    pure function l2norm(vec) result(norm)
      real(dp), intent(in) :: vec(:)
      real(dp) norm
      norm = sqrt(dot_product(vec, vec))
    end function l2norm

end module utilities

module isohard
  use utilities, only: dp
  implicit none
  private
  public isotropic_stiffness, elastic_increment, deviatoric, l2norm, &
        tensordot, calc_eldmg
  contains
    pure function isotropic_stiffness(E, v) result(C)
      ! Calculates a vector represeting the three stiffnesses
      ! for isotropic material
      real(dp), intent(in) :: E, v
      real(dp) :: C(3), C11, C12, C44

      C11 = E*(1.0-v)/((1.0+v)*(1.0-2.0*v))
      C12 = E*v/((1.0+v)*(1.0-2.0*v))
      C44 = E/(2.0*(1.0+v))
      C = (/C11, C12, C44/)

      return
    end function isotropic_stiffness


    pure function elastic_increment(C, de) result(S)
      real(dp), intent(in) :: C(3), de(6)
      real(dp) :: S(6)

      S(1) = C(1)*de(1) + C(2)*de(2) + C(2)*de(3)
      S(2) = C(2)*de(1) + C(1)*de(2) + C(2)*de(3)
      S(3) = C(2)*de(1) + C(2)*de(2) + C(1)*de(3)
      S(4:6) = 2 * C(3) * de(4:6)

      return
    end function elastic_increment

    pure function deviatoric(tensor) result(dev_tensor)
      ! Calculates the deviatoric part of the second order tensor.
      ! Assumes the tensor to be on Voigt vector notation
      ! Returns a vector using Voigt notation
      real(dp), intent(in) :: tensor(6)
      real(dp) :: dev_tensor(6), hyd

      hyd = sum(tensor(1:3))/3.d0
      dev_tensor(1:3) = tensor(1:3) - hyd
      dev_tensor(4:6) = tensor(4:6)

      return
    end function deviatoric

    pure function l2norm(tensor) result(scalar)
      real(dp), intent(in), dimension(6) :: tensor
      real(dp) :: scalar
      integer :: i

      scalar = 0.d0
      do i=1,6
        scalar = scalar + tensor(i)*tensor(i)
      end do

      scalar = dsqrt(scalar)

    end function l2norm

    pure function tensordot(tensor1, tensor2) result(scalar)
      real(dp), intent(in), dimension(6) :: tensor1, tensor2
      real(dp) :: scalar
      integer :: i

      scalar = 0.d0
      do i=1,3
        scalar = scalar + tensor1(i) * tensor2(i)
        scalar = scalar + 2 * tensor1(i+3) * tensor2(i+3)
      end do

    end function tensordot

    pure function calc_eldmg(stress) result(scalar)
      real(dp), intent(in), dimension(6) :: stress
      real(dp) :: scalar

      integer :: i

      scalar = 0.1
    end function calc_eldmg
end module isohard
