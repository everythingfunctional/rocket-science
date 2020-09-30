module refurbished_mod1
  implicit none

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
  ! assuminng a 1.1 drag coefficient and

  real(dp) :: accel = zero
  real(dp) :: altitude = zero
  real(dp) :: area
  real(dp) :: cf
  real(dp) :: cp
  real(dp) :: cv
  real(dp) :: db = zero
  real(dp) :: den ! air density
  real(dp) :: dia
  real(dp) :: drag = zero
  real(dp) :: dsigng
  real(dp) :: dt
  real(dp) :: echam
  real(dp) :: edotgen
  real(dp) :: edotos
  real(dp) :: edotout
  real(dp) :: energy
  real(dp) :: g
  integer :: i
  real(dp) :: id
  real(dp) :: length
  real(dp) :: mcham
  real(dp) :: mdotgen
  real(dp) :: mdotos = zero
  real(dp) :: mdotout
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
  real(dp) :: texit
  real(dp) :: Tflame
  real(dp) :: thrust = zero
  real(dp) :: time = zero
  real(dp) :: tmax
  real(dp) :: vol = one
  real(dp) :: vel = zero
end module

module refurbished
  implicit none
contains
  subroutine propwt ! calculate weight of propellent
    use refurbished_mod1, only: &
        dp, id, length, od, pi, propmass, rocketmass, rhos

    propmass = pi / 4.0_dp * (od**2 - id**2) * length * rhos
    rocketmass = 0.15_dp * propmass ! assume 85% propellant loading and 15% extra wt of rocket
  end subroutine

  subroutine burnrate
    use refurbished_mod1, only: db, dt, n, p, pref, r, rref

    r = rref * (p/pref)**n ! calculate burn rate
    db = db + r*dt ! calculate incremental burn distance
  end subroutine

  subroutine calcsurf
    ! cylinder burning from id outward and from both ends along the length
    use refurbished_mod1, only: db, dp, dt, id, length, od, pi, r, surf, vol

    surf = pi * (id + 2.0_dp*db) * (length - 2.0_dp*db) + pi * (od**2 - (id + 2.0_dp*db)**2) * 0.5_dp

    if (id + 2.0_dp*db .gt. od .or. db.gt.length/2.0_dp) then
      surf = 0.0_dp ! we hit the wall and burned out
      r = 0.0_dp ! turn off burn rate so burn distance stops increasing
    end if

    vol = vol + r*surf*dt ! increment the interior volume of the chamber a little
  end subroutine

  subroutine calmdotgen
    use refurbished_mod1, only: cp, edotgen, mdotgen, r, rhos, surf, Tflame

    mdotgen = rhos * r * surf
    edotgen = mdotgen * cp * Tflame
  end subroutine

  subroutine massflow
    use refurbished_mod1, only: &
        area, cp, dp, dsigng, edotos, g, mdotos, p, pamb, rgas, t, tamb

    real(dp) :: ax
    real(dp) :: cpx
    real(dp) :: cstar
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

  subroutine addmass
    use refurbished_mod1, only: &
        dt, echam, edotgen, edotos, mcham, mdotgen, mdotos

    mcham = mcham + (mdotgen - mdotos) * dt
    echam = echam + (edotgen - edotos) * dt
  end subroutine

  subroutine calct
    use refurbished_mod1, only: cv, echam, mcham, t

    t = echam / mcham / cv
  end subroutine

  subroutine calcp
    use refurbished_mod1, only: mcham, p, rgas, t, vol

    p = mcham * rgas * t / vol
  end subroutine

  subroutine calcthrust
    use refurbished_mod1, only: &
        altitude, &
        area, &
        cf, &
        cd, &
        den, &
        dp, &
        drag, &
        gravity, &
        mwair, &
        netthrust, &
        p, &
        pamb, &
        rhob, &
        RU, &
        surfrocket, &
        tamb, &
        thrust, &
        vel

    thrust = (p - pamb) * area * cf ! correction to thrust (actual vs vacuum thrust)
    den = rhob * exp(-gravity * mwair * altitude / RU / tamb)
    drag = -cd * 0.5_dp * den * vel * abs(vel) * surfrocket
    netthrust=thrust+drag
  end subroutine

  subroutine height
    use refurbished_mod1, only: &
        accel, &
        altitude, &
        dt, &
        gravity, &
        mcham, &
        mdotgen, &
        netthrust, &
        propmass, &
        rocketmass, &
        vel

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

    use refurbished_mod1, only: &
        accel, &
        altitude, &
        area, &
        cf, &
        cp, &
        cv, &
        dia, &
        dp, &
        drag, &
        dt, &
        echam, &
        g, &
        i, &
        id, &
        length, &
        mcham, &
        mdotos, &
        mw, &
        n, &
        netthrust, &
        nsteps, &
        od, &
        output, &
        p, &
        pamb, &
        pi, &
        pref, &
        psipa, &
        rgas, &
        rhos, &
        rref, &
        ru, &
        t, &
        Tflame, &
        thrust, &
        time, &
        tmax, &
        vel, &
        vol

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

    call propwt
    do i=1,nsteps
      call burnrate
      call calcsurf
      call calmdotgen  ! [mdot,engy,dsign]= massflow(p1,pamb,t1,tamb,cp,cp,rgas,rgas,g,g,area)
      call massflow
      call addmass
      call calct
      call calcp
      call calcthrust
      call height
      time = time + dt
      output(i,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]
    end do
    rocket = output
  end function
end module
