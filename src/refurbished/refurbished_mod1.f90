module refurbished_mod1
  implicit none

  integer, parameter :: precision = 15, range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)
  real(dp), parameter :: gravity = 9.81d0
  real(dp), parameter :: pi = 3.1415926539
  real(dp), parameter :: RU = 8314d0
  real(dp), parameter :: cd = 1.1, rhob = 1.225, tamb = 300d0
  real(dp), parameter :: mwair = 28.96, surfrocket = pi/4
  ! assuminng a 1.1 drag coefficient and

  real(dp) :: cp, cv, g, rgas, mw, vol, dia, cf, id, od, length, rref, rhos, psipa, pref
  real(dp) :: db, dt, tmax, Tflame
  real(dp) :: thrust, area, r, n, surf, mdotgen, mdotout, edotgen, edotout, energy
  real(dp) :: mdotos, edotos, texit, dsigng, pamb, p, t
  real(dp) :: mcham, echam, time, propmass, drag, netthrust
  integer nsteps, i
  real(dp) :: accel, vel, altitude, rocketmass
  real(dp), allocatable :: output(:,:)
  real(dp) den ! air density
end module
