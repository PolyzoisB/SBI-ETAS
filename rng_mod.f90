module rng_mod
  use iso_fortran_env, only: real64
  
  implicit none
  private

  public :: rng_init, rng_uniform, rng_poisson

  interface
    subroutine zbqlini(seed)
      integer, intent(in) :: seed
    end subroutine zbqlini

    function zbqlu01(dummy) result(value)
      import :: real64
      real(real64), intent(in) :: dummy
      real(real64) :: value
    end function zbqlu01

    function zbqlpoi(mu) result(value)
      import :: real64
      real(real64), intent(in) :: mu
      integer :: value
    end function zbqlpoi
  end interface

contains

  subroutine rng_init(seed)
    integer, intent(in) :: seed

    ! Seed zero invokes an old Unix-specific clock routine in randgen.f.
    if (seed == 0) then
      error stop 'The RNG seed must be nonzero.'
    end if

    call zbqlini(abs(seed))
  end subroutine rng_init

  function rng_uniform() result(value)
    real(real64) :: value

    value = zbqlu01(0.0_real64)
  end function rng_uniform

  function rng_poisson(mu) result(value)
    real(real64), intent(in) :: mu
    integer :: value

    if (mu < 0.0_real64) then
      error stop 'Poisson mean cannot be negative.'
    end if

    value = zbqlpoi(mu)
  end function rng_poisson

end module rng_mod