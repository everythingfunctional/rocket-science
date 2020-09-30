module refurbished_mod1
  implicit none

  integer, parameter :: precision = 15
  integer, parameter :: range = 307
  integer, parameter :: dp = selected_real_kind(precision, range)
  real(dp), parameter :: gravity = 9.81d0
  real(dp), parameter :: pi = 3.1415926539
  real(dp), parameter :: RU = 8314d0
  real(dp), parameter :: zero = 0._dp
  real(dp), parameter :: one = 1._dp
  real(dp), parameter :: cd = 1.1
  real(dp), parameter :: rhob = 1.225
  real(dp), parameter :: tamb = 300d0
  real(dp), parameter :: mwair = 28.96
  real(dp), parameter :: surfrocket = pi/4
  ! assuminng a 1.1 drag coefficient and

  real(dp) :: cp
  real(dp) :: cv
  real(dp) :: g
  real(dp) :: rgas
  real(dp) :: mw
  real(dp) :: vol = one
  real(dp) :: dia
  real(dp) :: cf
  real(dp) :: id
  real(dp) :: od
  real(dp) :: length
  real(dp) :: rref
  real(dp) :: rhos
  real(dp) :: psipa
  real(dp) :: pref
  real(dp) :: db = zero
  real(dp) :: dt
  real(dp) :: tmax
  real(dp) :: Tflame
  real(dp) :: thrust = zero
  real(dp) :: area
  real(dp) :: r
  real(dp) :: n
  real(dp) :: surf
  real(dp) :: mdotgen
  real(dp) :: mdotout
  real(dp) :: edotgen
  real(dp) :: edotout
  real(dp) :: energy
  real(dp) :: mdotos = zero
  real(dp) :: edotos
  real(dp) :: texit
  real(dp) :: dsigng
  real(dp) :: pamb
  real(dp) :: p
  real(dp) :: t
  real(dp) :: mcham
  real(dp) :: echam
  real(dp) :: time = zero
  real(dp) :: propmass = zero
  real(dp) :: drag = zero
  real(dp) :: netthrust = zero
  integer :: nsteps
  integer :: i
  real(dp) :: accel = zero
  real(dp) :: vel = zero
  real(dp) :: altitude = zero
  real(dp) :: rocketmass = zero
  real(dp), allocatable :: output(:,:)
  real(dp) :: den ! air density
end module

module refurbished
  implicit none
contains
  subroutine propwt ! calculate weight of propellent
    use refurbished_mod1
    implicit none

    propmass = pi / 4 * (od**2 - id**2) * length * rhos
    rocketmass = 0.15 * propmass ! assume 85% propellant loading and 15% extra wt of rocket
  end subroutine

  subroutine burnrate
    use refurbished_mod1
    implicit none

    r = rref * (p/pref)**n ! calculate burn rate
    db = db + r*dt ! calculate incremental burn distance
  end subroutine

  subroutine calcsurf
    ! cylinder burning from id outward and from both ends along the length
    use refurbished_mod1
    implicit none

    surf = pi * (id + 2.0d0*db) * (length - 2.0d0*db) + pi * (od**2.0d0 - (id + 2.0*db)**2.0d0) * 0.5

    if (id + 2d0*db .gt. od .or. db.gt.length/2d0) THEN
      surf = 0d0 ! we hit the wall and burned out
      r = 0 ! turn off burn rate so burn distance stops increasing
    endif

    vol = vol + r*surf*dt ! increment the interior volume of the chamber a little
  end subroutine

  subroutine calmdotgen
    use refurbished_mod1
    implicit none

    mdotgen = rhos * r * surf
    edotgen = mdotgen * cp * Tflame
  end subroutine

  subroutine massflow
    USE refurbished_mod1
    implicit none

    REAL(8) :: mdtx
    REAL(8) :: engyx
    REAL(8) :: tx
    REAL(8) :: gx
    REAL(8) :: rx
    REAL(8) :: px
    REAL(8) :: cpx
    REAL(8) :: pcrit
    REAL(8) :: facx
    REAL(8) :: term1
    REAL(8) :: term2
    REAL(8) :: pratio
    REAL(8) :: cstar
    REAL(8) :: ax
    REAL(8) :: hx
    REAL(8) :: p1
    REAL(8) :: p2

    mdotos = 0.
    edotos = 0. ! initially set them to zero prior to running this loop

    p1 = p
    p2 = pamb
    ax = area
    IF (p1.GT.p2) THEN
      dsigng = 1
      tx = t
      gx = g
      rx = rgas
      px = p
      cpx = cp
      hx = cp * t
      pratio = p1 / p2
    else
      dsigng = -1
      tx = tamb
      gx = g
      rx = rgas
      px = pamb
      cpx = cp
      hx = cp * tamb
      pratio = p2 / p1
    end if

    pcrit = (2. / (gx + 1.))**(gx / (gx - 1.))
    IF ((1. / pratio) .LT. pcrit) then
      ! choked flow
      cstar = sqrt((1. / gx) * ((gx + 1.) / 2.)**((gx + 1.) / (gx - 1.)) * rx * tx)
      mdtx= px * ax / cstar
    else
      ! unchoked flow
      facx = pratio**((gx - 1.) / gx)
      term1 = SQRT(gx * rx * tx / facx)
      term2 = SQRT((facx - 1.) / (gx - 1.))
      mdtx = SQRT(2.) * px / pratio / rx / tx * facx * term1 * term2 * ax
    end if
    engyx = mdtx * hx  ! reformulate based on enthalpy of the chamber
    mdotos = mdtx * dsigng ! exiting mass flow (could be negative "dsigng")
    edotos = engyx * dsigng ! exiting enthalpy
  end subroutine

  subroutine addmass
    use refurbished_mod1
    implicit none

    mcham = mcham + (mdotgen - mdotos) * dt
    echam = echam + (edotgen - edotos) * dt
  end subroutine

  subroutine calct
    use refurbished_mod1
    implicit none

    t = echam / mcham / cv
  end subroutine

  subroutine calcp
    use refurbished_mod1
    implicit none

    p = mcham * rgas * t / vol
  end subroutine

  subroutine calcthrust
    use refurbished_mod1
    implicit none

    thrust = (p - pamb) * area * cf ! correction to thrust (actual vs vacuum thrust)
    den = rhob * exp(-gravity * mwair * altitude / RU / tamb)
    drag = -cd * 0.5 * den * vel * abs(vel) * surfrocket
    netthrust=thrust+drag
  end subroutine

  subroutine height
    use refurbished_mod1
    implicit none

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

    use refurbished_mod1
    implicit none

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
    psipa = 6894.76d0 ! pascals per psi (constant)
    pref = 3000d0 * psipa ! reference pressure (constant)

    nsteps=nint(tmax/dt) ! number of time steps

    ! preallocate an output file for simulation infomration
    allocate(output(0:nsteps,11))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

    !! now begin calculating and initializing
    ! gas variables
    rgas = ru / mw
    cv = cp - rgas
    g = cp / cv

    area = pi / 4d0 * dia**2.0d0 ! nozzle area

    pamb = 101325d0 ! atmospheric pressure

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
