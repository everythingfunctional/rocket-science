module refurbished
  implicit none
  private

  public :: rocket, dp

  integer, parameter :: precision = 15
  integer, parameter :: range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)

  real(dp), parameter :: cd = 1.1_dp
  real(dp), parameter :: gravity = 9.81_dp
  real(dp), parameter :: mwair = 28.96_dp
  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: pi = 3.1415926539_dp
  real(dp), parameter :: rhob = 1.225_dp
  real(dp), parameter :: RU = 8314.0_dp
  real(dp), parameter :: surfrocket = pi / 4.0_dp
  real(dp), parameter :: tamb = 300.0_dp
  real(dp), parameter :: zero = 0.0_dp
contains
  subroutine propwt(id, length, od, rhos,  propmass, rocketmass)
    ! calculate weight of propellent
    real(dp), intent(in) :: id
    real(dp), intent(in) :: length
    real(dp), intent(in) :: od
    real(dp), intent(in) :: rhos
    real(dp), intent(out) :: propmass
    real(dp), intent(out) :: rocketmass

    propmass = pi / 4.0_dp * (od**2 - id**2) * length * rhos
    rocketmass = 0.15_dp * propmass ! assume 85% propellant loading and 15% extra wt of rocket
  end subroutine

  subroutine burnrate(dt, n, p, pref, rref,  db,  r)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: n
    real(dp), intent(in) :: p
    real(dp), intent(in) :: pref
    real(dp), intent(in) :: rref
    real(dp), intent(inout) :: db
    real(dp), intent(out) :: r

    r = rref * (p/pref)**n ! calculate burn rate
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

    surf = pi * (id + 2.0_dp*db) * (length - 2.0_dp*db) + pi * (od**2 - (id + 2.0_dp*db)**2) * 0.5_dp

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

  subroutine massflow(area, cp, g, p, pamb, rgas, t,  edotos, mdotos)
    real(dp), intent(in) :: area
    real(dp), intent(in) :: cp
    real(dp), intent(in) :: g
    real(dp), intent(in) :: p
    real(dp), intent(in) :: pamb
    real(dp), intent(in) :: rgas
    real(dp), intent(in) :: t
    real(dp), intent(out) :: edotos
    real(dp), intent(out) :: mdotos

    real(dp) :: ax
    real(dp) :: cpx
    real(dp) :: cstar
    real(dp) :: dsigng
    real(dp) :: engyx
    real(dp) :: facx
    real(dp) :: gx
    real(dp) :: hx
    real(dp) :: mdtx
    real(dp) :: p1
    real(dp) :: p2
    real(dp) :: pcrit
    real(dp) :: pratio
    real(dp) :: px
    real(dp) :: rx
    real(dp) :: term1
    real(dp) :: term2
    real(dp) :: tx

    mdotos = 0.0_dp
    edotos = 0.0_dp ! initially set them to zero prior to running this loop

    p1 = p
    p2 = pamb
    ax = area
    if (p1.GT.p2) then
      dsigng = 1.0_dp
      tx = t
      gx = g
      rx = rgas
      px = p
      cpx = cp
      hx = cp * t
      pratio = p1 / p2
    else
      dsigng = -1.0_dp
      tx = tamb
      gx = g
      rx = rgas
      px = pamb
      cpx = cp
      hx = cp * tamb
      pratio = p2 / p1
    end if

    pcrit = (2.0_dp / (gx + 1.0_dp))**(gx / (gx - 1.0_dp))
    if ((1.0_dp / pratio) .LT. pcrit) then
      ! choked flow
      cstar = sqrt((1.0_dp / gx) * ((gx + 1.0_dp) / 2.0_dp)**((gx + 1.0_dp) / (gx - 1.0_dp)) * rx * tx)
      mdtx= px * ax / cstar
    else
      ! unchoked flow
      facx = pratio**((gx - 1.0_dp) / gx)
      term1 = sqrt(gx * rx * tx / facx)
      term2 = sqrt((facx - 1.0_dp) / (gx - 1.0_dp))
      mdtx = sqrt(2.0_dp) * px / pratio / rx / tx * facx * term1 * term2 * ax
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
      altitude, area, cf, p, pamb, vel,  drag, netthrust, thrust)
    real(dp), intent(in) :: altitude
    real(dp), intent(in) :: area
    real(dp), intent(in) :: cf
    real(dp), intent(in) :: p
    real(dp), intent(in) :: pamb
    real(dp), intent(in) :: vel
    real(dp), intent(out) :: drag
    real(dp), intent(out) :: netthrust
    real(dp), intent(out) :: thrust

    real(dp) :: den

    thrust = (p - pamb) * area * cf ! correction to thrust (actual vs vacuum thrust)
    den = rhob * exp(-gravity * mwair * altitude / RU / tamb)
    drag = -cd * 0.5_dp * den * vel * abs(vel) * surfrocket
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
    accel = netthrust / (propmass + rocketmass + mcham) - gravity
    vel = vel + accel*dt
    altitude = altitude + vel*dt
  end subroutine

  function rocket( &
      dt_, &
      t_max_, &
      c_p_, &
      MW_, &
      temperature_, &
      pressure_, &
      T_flame_, &
      r_ref_, &
      n_, &
      id_, &
      od_, &
      length_, &
      rho_solid_, &
      dia_, &
      C_f_)
    !! this is a basic program of a single stage
    !! rocket motor flowing out of a nozzle, assuming
    !! a thrust coefficient and ignoring the complexities of
    !! what happens to thrust at low pressures, i.e. shock in the nozzle
    real(dp), intent(in) :: dt_
    real(dp), intent(in) :: t_max_
    real(dp), intent(in) :: c_p_
    real(dp), intent(in) :: MW_
    real(dp), intent(in) :: temperature_
    real(dp), intent(in) :: pressure_
    real(dp), intent(in) :: T_flame_
    real(dp), intent(in) :: r_ref_
    real(dp), intent(in) :: n_
    real(dp), intent(in) :: id_
    real(dp), intent(in) :: od_
    real(dp), intent(in) :: length_
    real(dp), intent(in) :: rho_solid_
    real(dp), intent(in) :: dia_
    real(dp), intent(in) :: C_f_
    real(dp), allocatable :: rocket(:,:)

    real(dp) :: accel = zero
    real(dp) :: altitude = zero
    real(dp) :: area
    real(dp) :: cf
    real(dp) :: cp
    real(dp) :: cv
    real(dp) :: db = zero
    real(dp) :: dia
    real(dp) :: drag = zero
    real(dp) :: dt
    real(dp) :: echam
    real(dp) :: edotgen
    real(dp) :: edotos
    real(dp) :: g
    integer :: i
    real(dp) :: id
    real(dp) :: length
    real(dp) :: mcham
    real(dp) :: mdotgen
    real(dp) :: mdotos = zero
    real(dp) :: mw
    real(dp) :: n
    real(dp) :: netthrust = zero
    integer :: nsteps
    real(dp) :: od
    real(dp), allocatable :: output(:,:)
    real(dp) :: p
    real(dp) :: pamb
    real(dp) :: pref
    real(dp) :: propmass = zero
    real(dp) :: psipa
    real(dp) :: rgas
    real(dp) :: r
    real(dp) :: rhos
    real(dp) :: rocketmass = zero
    real(dp) :: rref
    real(dp) :: surf
    real(dp) :: t
    real(dp) :: Tflame
    real(dp) :: thrust = zero
    real(dp) :: time = zero
    real(dp) :: tmax
    real(dp) :: vol = one
    real(dp) :: vel = zero

    dt = dt_
    tmax = t_max_
    cp = c_p_
    mw = MW_
    t = temperature_
    p = pressure_
    Tflame = T_flame_
    rref = r_ref_
    n = n_
    id = id_
    od = od_
    length = length_
    rhos = rho_solid_
    dia = dia_
    cf = C_f_

    ! propellent grain is a cylinder burning radially outward and axially inward from one end.
    ! the other end is considered inhibited.
    ! outer diameter is inhibited because this is a cast propellent: it was poured
    ! into the tube/chamber and only the inner diameter burns when ignited.

    ! propellant burn rate information
    psipa = 6894.76_dp ! pascals per psi (constant)
    pref = 3000.0_dp * psipa ! reference pressure (constant)

    nsteps = nint(tmax / dt) ! number of time steps

    ! preallocate an output file for simulation infomration
    allocate(output(0:nsteps, 11))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

    !! now begin calculating and initializing
    ! gas variables
    rgas = ru / mw
    cv = cp - rgas
    g = cp / cv

    area = pi / 4.0_dp * dia**2.0_dp ! nozzle area

    pamb = 101325.0_dp ! atmospheric pressure

    ! calculate initial mass and energy in the chamber
    mcham = p * vol / rgas / t ! use ideal gas law to determine mass in chamber
    echam = mcham * cv * t ! initial internal energy in chamber

    output(0,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]

    call propwt(id, length, od, rhos,  propmass, rocketmass)
    do i=1,nsteps
      call burnrate(dt, n, p, pref, rref,  db,  r)
      call calcsurf(db, dt, id, length, od,  vol,  r, surf)
      call calmdotgen(cp, r, rhos, surf, Tflame,  edotgen, mdotgen) ! [mdot,engy,dsign]= massflow(p1,pamb,t1,tamb,cp,cp,rgas,rgas,g,g,area)
      call massflow(area, cp, g, p, pamb, rgas, t,  edotos, mdotos)
      call addmass(dt, edotgen, edotos, mdotgen, mdotos,  echam, mcham)
      call calct(cv, echam, mcham,  t)
      call calcp(mcham, rgas, t, vol,  p)
      call calcthrust( &
          altitude, area, cf, p, pamb, vel,  drag, netthrust, thrust)
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
      output(i,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]
    end do
    rocket = output
  end function
end module
