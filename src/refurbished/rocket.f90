module refurbished
  implicit none
  private

  public :: rocket, dp

  integer, parameter :: precision = 15
  integer, parameter :: range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)

  real(dp), parameter :: AMBIENT_TEMPERATURE = 300.0_dp
  real(dp), parameter :: ATMOSPHERIC_PRESSURE = 101325.0_dp
  real(dp), parameter :: GRAVITY = 9.81_dp
  real(dp), parameter :: PI = 3.1415926539_dp
  real(dp), parameter :: UNIVERSAL_GAS_CONSTANT = 8314.0_dp
contains
  subroutine initialize_fuel_and_rocket_mass( &
      fuel_inner_diameter, &
      fuel_length, &
      fuel_outer_diameter, &
      fuel_density,  &

      fuel_mass, &
      rocket_mass)
    real(dp), intent(in) :: fuel_inner_diameter
    real(dp), intent(in) :: fuel_length
    real(dp), intent(in) :: fuel_outer_diameter
    real(dp), intent(in) :: fuel_density
    real(dp), intent(out) :: fuel_mass
    real(dp), intent(out) :: rocket_mass

    associate( &
        id => fuel_inner_diameter, &
        od => fuel_outer_diameter, &
        l => fuel_length)
      associate(fuel_volume => PI * (od**2 - id**2) / 4.0_dp * l)
        fuel_mass = fuel_volume * fuel_density
      end associate
    end associate
    rocket_mass = 0.15_dp * fuel_mass
  end subroutine

  subroutine burnrate(dt, n, p, rref,  db,  r)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: n
    real(dp), intent(in) :: p
    real(dp), intent(in) :: rref
    real(dp), intent(inout) :: db
    real(dp), intent(out) :: r

    real(dp), parameter :: PASCALS_PER_PSI = 6894.76_dp
    real(dp), parameter :: REFERENCE_PRESSURE = 3000.0_dp * PASCALS_PER_PSI

    r = rref * (p/REFERENCE_PRESSURE)**n ! calculate burn rate
    db = db + r*dt ! calculate incremental burn distance
  end subroutine

  subroutine calcsurf(db, dt, id, length, od,  vol,  r, surf)
    ! cylinder burning from id outward and from both ends along the length
    real(dp), intent(in) :: db
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: id
    real(dp), intent(in) :: length
    real(dp), intent(in) :: od
    real(dp), intent(inout) :: vol
    real(dp), intent(out) :: r
    real(dp), intent(out) :: surf

    surf = PI * (id + 2.0_dp*db) * (length - 2.0_dp*db) + PI * (od**2 - (id + 2.0_dp*db)**2) * 0.5_dp

    if (id + 2.0_dp*db .gt. od .or. db.gt.length/2.0_dp) then
      surf = 0.0_dp ! we hit the wall and burned out
      r = 0.0_dp ! turn off burn rate so burn distance stops increasing
    end if

    vol = vol + r*surf*dt ! increment the interior volume of the chamber a little
  end subroutine

  subroutine calmdotgen(cp, r, rhos, surf, Tflame,  edotgen, mdotgen)
    real(dp), intent(in) :: cp
    real(dp), intent(in) :: r
    real(dp), intent(in) :: rhos
    real(dp), intent(in) :: surf
    real(dp), intent(in) :: Tflame
    real(dp), intent(out) :: edotgen
    real(dp), intent(out) :: mdotgen

    mdotgen = rhos * r * surf
    edotgen = mdotgen * cp * Tflame
  end subroutine

  subroutine massflow(area, cp, g, p, rgas, t,  edotos, mdotos)
    real(dp), intent(in) :: area
    real(dp), intent(in) :: cp
    real(dp), intent(in) :: g
    real(dp), intent(in) :: p
    real(dp), intent(in) :: rgas
    real(dp), intent(in) :: t
    real(dp), intent(out) :: edotos
    real(dp), intent(out) :: mdotos

    real(dp) :: cstar
    real(dp) :: dsigng
    real(dp) :: engyx
    real(dp) :: facx
    real(dp) :: hx
    real(dp) :: mdtx
    real(dp) :: pcrit
    real(dp) :: pratio
    real(dp) :: px
    real(dp) :: term1
    real(dp) :: term2
    real(dp) :: tx

    mdotos = 0.0_dp
    edotos = 0.0_dp ! initially set them to zero prior to running this loop

    if (p.GT.ATMOSPHERIC_PRESSURE) then
      dsigng = 1.0_dp
      tx = t
      px = p
      hx = cp * t
      pratio = p / ATMOSPHERIC_PRESSURE
    else
      dsigng = -1.0_dp
      tx = AMBIENT_TEMPERATURE
      px = ATMOSPHERIC_PRESSURE
      hx = cp * AMBIENT_TEMPERATURE
      pratio = ATMOSPHERIC_PRESSURE / p
    end if

    pcrit = (2.0_dp / (g + 1.0_dp))**(g / (g - 1.0_dp))
    if ((1.0_dp / pratio) .LT. pcrit) then
      ! choked flow
      cstar = sqrt((1.0_dp / g) * ((g + 1.0_dp) / 2.0_dp)**((g + 1.0_dp) / (g - 1.0_dp)) * rgas * tx)
      mdtx= px * area / cstar
    else
      ! unchoked flow
      facx = pratio**((g - 1.0_dp) / g)
      term1 = sqrt(g * rgas * tx / facx)
      term2 = sqrt((facx - 1.0_dp) / (g - 1.0_dp))
      mdtx = sqrt(2.0_dp) * px / pratio / rgas / tx * facx * term1 * term2 * area
    end if
    engyx = mdtx * hx  ! reformulate based on enthalpy of the chamber
    mdotos = mdtx * dsigng ! exiting mass flow (could be negative "dsigng")
    edotos = engyx * dsigng ! exiting enthalpy
  end subroutine

  subroutine addmass(dt, edotgen, edotos, mdotgen, mdotos,  echam, mcham)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: edotgen
    real(dp), intent(in) :: edotos
    real(dp), intent(in) :: mdotgen
    real(dp), intent(in) :: mdotos
    real(dp), intent(inout) :: echam
    real(dp), intent(inout) :: mcham

    mcham = mcham + (mdotgen - mdotos) * dt
    echam = echam + (edotgen - edotos) * dt
  end subroutine

  subroutine calct(cv, echam, mcham,  t)
    real(dp), intent(in) :: cv
    real(dp), intent(in) :: echam
    real(dp), intent(in) :: mcham
    real(dp), intent(out) :: t

    t = echam / mcham / cv
  end subroutine

  subroutine calcp(mcham, rgas, t, vol,  p)
    real(dp), intent(in) :: mcham
    real(dp), intent(in) :: rgas
    real(dp), intent(in) :: t
    real(dp), intent(in) :: vol
    real(dp), intent(out) :: p

    p = mcham * rgas * t / vol
  end subroutine

  subroutine calcthrust( &
      altitude, area, cf, p, vel,  drag, netthrust, thrust)
    real(dp), intent(in) :: altitude
    real(dp), intent(in) :: area
    real(dp), intent(in) :: cf
    real(dp), intent(in) :: p
    real(dp), intent(in) :: vel
    real(dp), intent(out) :: drag
    real(dp), intent(out) :: netthrust
    real(dp), intent(out) :: thrust

    real(dp), parameter :: DRAG_COEFFICIENT = 1.1_dp
    real(dp), parameter :: MOLECULAR_WEIGHT_OF_AIR = 28.96_dp
    real(dp), parameter :: REFERENCE_AIR_DENSITY = 1.225_dp
    real(dp), parameter :: ROCKET_SURFACE_AREA = PI / 4.0_dp
    real(dp) :: den

    thrust = (p - ATMOSPHERIC_PRESSURE) * area * cf ! correction to thrust (actual vs vacuum thrust)
    den = REFERENCE_AIR_DENSITY * exp(-GRAVITY * MOLECULAR_WEIGHT_OF_AIR * altitude / UNIVERSAL_GAS_CONSTANT / AMBIENT_TEMPERATURE)
    drag = -DRAG_COEFFICIENT * 0.5_dp * den * vel * abs(vel) * ROCKET_SURFACE_AREA
    netthrust=thrust+drag
  end subroutine

  subroutine height( &
      dt, &
      mcham, &
      mdotgen, &
      netthrust, &
      rocketmass, &

      altitude, &
      propmass, &
      vel, &

      accel)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: mcham
    real(dp), intent(in) :: mdotgen
    real(dp), intent(in) :: netthrust
    real(dp), intent(in) :: rocketmass
    real(dp), intent(inout) :: altitude
    real(dp), intent(inout) :: propmass
    real(dp), intent(inout) :: vel
    real(dp), intent(out) :: accel

    propmass = propmass - mdotgen*dt ! incremental change in propellant mass
    accel = netthrust / (propmass + rocketmass + mcham) - GRAVITY
    vel = vel + accel*dt
    altitude = altitude + vel*dt
  end subroutine

  function rocket( &
      dt, &
      tmax, &
      cp, &
      mw, &
      temperature_, &
      pressure_, &
      Tflame, &
      rref, &
      n, &
      id, &
      od, &
      length, &
      rhos, &
      dia, &
      cf)
    !! this is a basic program of a single stage
    !! rocket motor flowing out of a nozzle, assuming
    !! a thrust coefficient and ignoring the complexities of
    !! what happens to thrust at low pressures, i.e. shock in the nozzle
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: tmax
    real(dp), intent(in) :: cp
    real(dp), intent(in) :: mw
    real(dp), intent(in) :: temperature_
    real(dp), intent(in) :: pressure_
    real(dp), intent(in) :: Tflame
    real(dp), intent(in) :: rref
    real(dp), intent(in) :: n
    real(dp), intent(in) :: id
    real(dp), intent(in) :: od
    real(dp), intent(in) :: length
    real(dp), intent(in) :: rhos
    real(dp), intent(in) :: dia
    real(dp), intent(in) :: cf
    real(dp), allocatable :: rocket(:,:)

    real(dp), parameter :: ONE = 1.0_dp
    real(dp), parameter :: ZERO = 0.0_dp
    real(dp) :: accel = ZERO
    real(dp) :: altitude = ZERO
    real(dp) :: area
    real(dp) :: cv
    real(dp) :: db = ZERO
    real(dp) :: drag = ZERO
    real(dp) :: echam
    real(dp) :: edotgen
    real(dp) :: edotos
    real(dp) :: g
    integer :: i
    real(dp) :: mcham
    real(dp) :: mdotgen
    real(dp) :: mdotos = ZERO
    real(dp) :: netthrust = ZERO
    integer :: nsteps
    real(dp) :: p
    real(dp) :: propmass = ZERO
    real(dp) :: rgas
    real(dp) :: r
    real(dp) :: rocketmass = ZERO
    real(dp) :: surf
    real(dp) :: t
    real(dp) :: thrust = ZERO
    real(dp) :: time = ZERO
    real(dp) :: vol = ONE
    real(dp) :: vel = ZERO

    t = temperature_
    p = pressure_

    ! propellent grain is a cylinder burning radially outward and axially inward from one end.
    ! the other end is considered inhibited.
    ! outer diameter is inhibited because this is a cast propellent: it was poured
    ! into the tube/chamber and only the inner diameter burns when ignited.

    ! propellant burn rate information

    nsteps = nint(tmax / dt) ! number of time steps

    ! preallocate an output file for simulation infomration
    allocate(rocket(0:nsteps, 11))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

    !! now begin calculating and initializing
    ! gas variables
    rgas = UNIVERSAL_GAS_CONSTANT / mw
    cv = cp - rgas
    g = cp / cv

    area = PI / 4.0_dp * dia**2 ! nozzle area

    ! calculate initial mass and energy in the chamber
    mcham = p * vol / rgas / t ! use ideal gas law to determine mass in chamber
    echam = mcham * cv * t ! initial internal energy in chamber

    rocket(0,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]

    call initialize_fuel_and_rocket_mass( &
        id, length, od, rhos,  propmass, rocketmass)
    do i=1,nsteps
      call burnrate(dt, n, p, rref,  db,  r)
      call calcsurf(db, dt, id, length, od,  vol,  r, surf)
      call calmdotgen(cp, r, rhos, surf, Tflame,  edotgen, mdotgen) ! [mdot,engy,dsign]= massflow(p1,ATMOSPHERIC_PRESSURE,t1,AMBIENT_TEMPERATURE,cp,cp,rgas,rgas,g,g,area)
      call massflow(area, cp, g, p, rgas, t,  edotos, mdotos)
      call addmass(dt, edotgen, edotos, mdotgen, mdotos,  echam, mcham)
      call calct(cv, echam, mcham,  t)
      call calcp(mcham, rgas, t, vol,  p)
      call calcthrust( &
          altitude, area, cf, p, vel,  drag, netthrust, thrust)
      call height( &
          dt, &
          mcham, &
          mdotgen, &
          netthrust, &
          rocketmass, &

          altitude, &
          propmass, &
          vel, &

          accel)
      time = time + dt
      rocket(i,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]
    end do
  end function
end module
