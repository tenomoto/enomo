module fft_module

	use type_module, only: i4b, i8b, dp
	include "fftw3.f"
!	implicit none
	private

	integer(kind=i8b), private :: plan_forward, plan_backward

	public :: fft_init, fft_clean, fft_analysis, fft_synthesis, fft_derivative

contains

	subroutine fft_init(n)
		implicit none

		integer(kind=i4b), intent(in) ::  n

    real(kind=dp), dimension(n) :: g
    complex(kind=dp), dimension(n) :: w

		call dfftw_plan_dft_r2c_1d(plan_forward, n, g, w, FFTW_ESTIMATE)
		call dfftw_plan_dft_c2r_1d(plan_backward, n, w, g, FFTW_ESTIMATE)

	end subroutine fft_init

	subroutine fft_clean()
		implicit none

		call dfftw_destroy_plan(plan_forward)
		call dfftw_destroy_plan(plan_backward)

	end subroutine fft_clean

	subroutine fft_analysis(g, w)
		implicit none

		real(kind=dp), dimension(:,:), intent(in) :: g 
		complex(kind=dp), dimension(:,:), intent(inout) :: w 
		integer(kind=i4b) :: j

		do j=1, size(g,2)
			call dfftw_execute_dft_r2c(plan_forward, g(:,j), w(:,j))
		end do
		w = w/size(g,1)

	end subroutine fft_analysis

	subroutine fft_synthesis(w, g)
		implicit none

		complex(kind=dp), dimension(:,:), intent(in) :: w 
		real(kind=dp), dimension(:,:), intent(inout) :: g 
		integer(kind=i4b) :: j

		do j=1, size(g,2)
			call dfftw_execute_dft_c2r(plan_backward, w(:,j), g(:,j))
		end do

	end subroutine fft_synthesis

  subroutine fft_derivative(g, gx, trunc)
    implicit none

		real(kind=dp), dimension(:,:), intent(in) :: g 
		real(kind=dp), dimension(:,:), intent(inout) :: gx
    integer(kind=i4b), optional, intent(in) :: trunc

    complex(kind=dp), dimension(size(g,1)) :: w
		integer(kind=i4b) :: n, m, j, nt

    n = size(g,1)
    if (present(trunc).and.(trunc<n/2)) then
      nt = trunc
    else
      nt = trunc
    end if
		do j=1, size(g,2)
			call dfftw_execute_dft_r2c(plan_forward, g(:,j), w)
			do m=0, nt
        w(m+1) = -cmplx(0.0_dp,m)*w(m+1)/n
      end do
      w(nt+2:n) = 0.0_dp 
			call dfftw_execute_dft_c2r(plan_backward, w, gx(:,j))
		end do

  end subroutine fft_derivative

end module fft_module
