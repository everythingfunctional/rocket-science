module refurbished_mod1
  implicit none

  integer, parameter :: precision = 15, range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)
  real(dp), parameter :: gravity = 9.81d0
  real(dp), parameter :: pi = 3.1415926539
  real(dp), parameter :: RU = 8314d0
  real(dp), parameter :: cd = 1.1, rhob = 1.225, tamb = 300d0
  ! assuminng a 1.1 drag coefficient and
  real(dp), parameter :: mwair = 28.96, surfrocket = pi/4

  real(dp) :: cp, cv, g, rgas, mw, dia, cf, id, od, length, rref, rhos, psipa, pref
  real(dp) :: dt, tmax, Tflame
  real(dp) :: area, r, n, surf, mdotgen, mdotout, edotgen, edotout, energy
  real(dp) :: edotos, texit, dsigng, pamb, p, t
  real(dp) :: mcham, echam
  integer nsteps, i
  real(dp), allocatable :: output(:,:)
  real(dp) den ! air density
end module
