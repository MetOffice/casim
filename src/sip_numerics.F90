!> Numerical utilities used by the Phillips et al. secondary ice
!> production (SIP) schemes in ice_multiplication:
!>
!>  * fixed-order (10-point) Gauss-Legendre quadrature in one dimension
!>    and over a two-dimensional region with variable inner limits;
!>  * the regularised lower incomplete gamma function P(a,x) and its
!>    inverse in x, used to set the upper truncation of the gamma
!>    particle size distributions in the collision integrals.
!>
!> References
!>  NIST Digital Library of Mathematical Functions (DLMF),
!>    https://dlmf.nist.gov/ : sections 3.5 (Gauss-Legendre quadrature),
!>    3.10 (continued fractions), 8.2, 8.7 and 8.9 (incomplete gamma).
!>  Abramowitz, M. and Stegun, I. A. (1964), Handbook of Mathematical
!>    Functions: 26.2.23 (normal quantile), 26.4.17 (Wilson-Hilferty).
!>
!> Contributed by the University of Manchester under the Horizon Europe
!> project CERTAINTY (grant agreement 101137680).
module sip_numerics

  use variable_precision, only: wp

  implicit none

  private

  ! Positive abscissae and weights of the 10-point Gauss-Legendre rule on
  ! [-1,1] (the rule is symmetric about zero).
  integer, parameter :: n_half = 5
  real(wp), parameter :: gl_node(n_half) = (/                                  &
       0.14887433898163122_wp, 0.43339539412924720_wp, 0.67940956829902440_wp, &
       0.86506336668898450_wp, 0.97390652851717170_wp /)
  real(wp), parameter :: gl_weight(n_half) = (/                                &
       0.29552422471475280_wp, 0.26926671930999650_wp, 0.21908636251598200_wp, &
       0.14945134915058040_wp, 0.06667134430868814_wp /)

  integer, parameter :: max_terms = 1000   ! cap on series / continued fraction
  real(wp), parameter :: rel_tol = 1.0e-14_wp

  abstract interface
    !> Integrand of one variable, evaluated at a vector of abscissae
    function integrand_1d(x) result(f)
      import :: wp
      real(wp), intent(in) :: x(:)
      real(wp) :: f(size(x))
    end function integrand_1d

    !> Integrand of two variables: scalar outer coordinate x and a vector
    !> of inner abscissae y
    function integrand_2d(x, y) result(f)
      import :: wp
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y(:)
      real(wp) :: f(size(y))
    end function integrand_2d

    !> Inner integration limit as a function of the outer coordinate
    function limit_function(x) result(y)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: y
    end function limit_function
  end interface

  public :: integrand_1d, integrand_2d, limit_function
  public :: gl_quad_1d, gl_quad_2d, gamma_p, inverse_gamma_p

contains

  !> Integral of f(x) from a to b with the 10-point Gauss-Legendre rule.
  function gl_quad_1d(f, a, b) result(total)

    implicit none

    procedure(integrand_1d) :: f
    real(wp), intent(in) :: a, b
    real(wp) :: total

    real(wp) :: centre, half_width
    real(wp) :: abscissa(2*n_half), weight(2*n_half)

    centre = 0.5_wp*(a + b)
    half_width = 0.5_wp*(b - a)

    abscissa(1:n_half) = centre - half_width*gl_node
    abscissa(n_half+1:2*n_half) = centre + half_width*gl_node
    weight(1:n_half) = gl_weight
    weight(n_half+1:2*n_half) = gl_weight

    total = half_width*sum(weight*f(abscissa))

  end function gl_quad_1d

  !> Integral over x from a to b, and over y from y_low(x) to y_high(x),
  !> of f(x,y), using the 10-point Gauss-Legendre rule in each direction.
  !> The outer coordinate is passed to the integrand explicitly, so the
  !> routine holds no module state and is safe to call from OpenMP threads.
  function gl_quad_2d(f, y_low, y_high, a, b) result(total)

    implicit none

    procedure(integrand_2d) :: f
    procedure(limit_function) :: y_low, y_high
    real(wp), intent(in) :: a, b
    real(wp) :: total

    real(wp) :: centre, half_width, x_outer, inner
    real(wp) :: y_centre, y_half_width
    real(wp) :: y_abscissa(2*n_half), weight(2*n_half)
    integer :: i

    weight(1:n_half) = gl_weight
    weight(n_half+1:2*n_half) = gl_weight

    centre = 0.5_wp*(a + b)
    half_width = 0.5_wp*(b - a)

    total = 0.0_wp
    do i = 1, 2*n_half
      if (i <= n_half) then
        x_outer = centre - half_width*gl_node(i)
      else
        x_outer = centre + half_width*gl_node(i-n_half)
      end if

      y_centre = 0.5_wp*(y_high(x_outer) + y_low(x_outer))
      y_half_width = 0.5_wp*(y_high(x_outer) - y_low(x_outer))
      y_abscissa(1:n_half) = y_centre - y_half_width*gl_node
      y_abscissa(n_half+1:2*n_half) = y_centre + y_half_width*gl_node

      inner = y_half_width*sum(weight*f(x_outer, y_abscissa))
      total = total + weight(i)*inner
    end do
    total = half_width*total

  end function gl_quad_2d

  !> Regularised lower incomplete gamma function P(a,x), a > 0, x >= 0.
  !> For x < a+1 the power series DLMF 8.7.1 is summed directly; otherwise
  !> P = 1 - Q with Q(a,x) from the continued fraction DLMF 8.9.2,
  !> evaluated with the modified Lentz algorithm (DLMF 3.10).
  function gamma_p(a, x) result(p)

    implicit none

    real(wp), intent(in) :: a, x
    real(wp) :: p

    real(wp), parameter :: tiny_value = 1.0e-300_wp
    real(wp) :: log_prefactor, term, series, denom
    real(wp) :: b_n, a_n, c, d, delta, fraction
    integer :: n

    if (x <= 0.0_wp .or. a <= 0.0_wp) then
      p = 0.0_wp
      return
    end if

    log_prefactor = a*log(x) - x - log_gamma(a)

    if (x < a + 1.0_wp) then
      ! P(a,x) = x^a e^-x / Gamma(a+1) * sum_n x^n / ((a+1)...(a+n))
      term = 1.0_wp/a
      series = term
      denom = a
      do n = 1, max_terms
        denom = denom + 1.0_wp
        term = term*x/denom
        series = series + term
        if (abs(term) < rel_tol*abs(series)) exit
      end do
      p = min(1.0_wp, series*exp(log_prefactor))
    else
      ! Q(a,x) = x^a e^-x / Gamma(a) * 1/(x+1-a- 1(1-a)/(x+3-a- 2(2-a)/(...)))
      b_n = x + 1.0_wp - a
      c = 1.0_wp/tiny_value
      d = 1.0_wp/b_n
      fraction = d
      do n = 1, max_terms
        a_n = -real(n, wp)*(real(n, wp) - a)
        b_n = b_n + 2.0_wp
        d = b_n + a_n*d
        if (abs(d) < tiny_value) d = tiny_value
        c = b_n + a_n/c
        if (abs(c) < tiny_value) c = tiny_value
        d = 1.0_wp/d
        delta = c*d
        fraction = fraction*delta
        if (abs(delta - 1.0_wp) < rel_tol) exit
      end do
      p = max(0.0_wp, 1.0_wp - fraction*exp(log_prefactor))
    end if

  end function gamma_p

  !> Inverse of P(a,x) in x: returns x such that P(a,x) = prob,
  !> for a > 0 and 0 < prob < 1.  The starting value is the Wilson-Hilferty
  !> approximation (A&S 26.4.17) with the normal quantile from A&S 26.2.23;
  !> the root is then bracketed and refined by Newton's method, falling
  !> back to bisection whenever a Newton step would leave the bracket.
  function inverse_gamma_p(prob, a) result(x)

    implicit none

    real(wp), intent(in) :: prob, a
    real(wp) :: x

    integer, parameter :: max_iter = 200
    real(wp) :: q, t, z, x_low, x_high, residual, density, x_new
    integer :: iter

    if (prob <= 0.0_wp .or. a <= 0.0_wp) then
      x = 0.0_wp
      return
    end if

    ! Standard normal quantile z with Phi(z) = prob
    q = min(prob, 1.0_wp - prob)
    q = max(q, 1.0e-300_wp)
    t = sqrt(-2.0_wp*log(q))
    z = t - (2.515517_wp + t*(0.802853_wp + t*0.010328_wp)) /                  &
            (1.0_wp + t*(1.432788_wp + t*(0.189269_wp + t*0.001308_wp)))
    if (prob < 0.5_wp) z = -z

    ! Wilson-Hilferty starting value
    x = a*(1.0_wp - 1.0_wp/(9.0_wp*a) + z/(3.0_wp*sqrt(a)))**3
    if (x <= 0.0_wp) x = 0.5_wp*a

    ! Bracket the root
    x_low = 0.0_wp
    x_high = max(x, a, 1.0_wp)
    do iter = 1, max_iter
      if (gamma_p(a, x_high) >= prob) exit
      x_low = x_high
      x_high = 2.0_wp*x_high
    end do
    if (prob >= 1.0_wp) then
      x = x_high
      return
    end if
    x = min(max(x, x_low), x_high)

    ! Safeguarded Newton iteration
    do iter = 1, max_iter
      residual = gamma_p(a, x) - prob
      if (residual > 0.0_wp) then
        x_high = x
      else
        x_low = x
      end if

      density = 0.0_wp
      if (x > 0.0_wp) density = exp((a - 1.0_wp)*log(x) - x - log_gamma(a))

      if (density > 0.0_wp) then
        x_new = x - residual/density
      else
        x_new = 0.5_wp*(x_low + x_high)
      end if
      if (x_new <= x_low .or. x_new >= x_high) x_new = 0.5_wp*(x_low + x_high)

      if (abs(x_new - x) <= 1.0e-12_wp*max(x_new, 1.0e-30_wp)) then
        x = x_new
        exit
      end if
      x = x_new
    end do

  end function inverse_gamma_p

end module sip_numerics
