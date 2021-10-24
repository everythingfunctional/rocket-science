module refurbished
  use refurbished_mod1, only : pi, dp, tamb, Ru, gravity, rhob, cd, mwair, surfrocket
  implicit none
contains
  subroutine propwt(od, id, length, rhos, propmass, rocketmass)
    !! calculate weight of propellent
    implicit none
    real(dp) id, od, length, rhos, propmass, rocketmass

    propmass = pi / 4 * (od**2 - id**2) * length * rhos
    rocketmass = 0.15 * propmass ! assume 85% propellant loading and 15% extra wt of rocket
  end subroutine

  subroutine burnrate(rref, p, pref, n, dt, r, db)
    implicit none
    real(dp) rref, p, pref, n, dt, r, db

    r = rref * (p/pref)**n ! calculate burn rate
    db = db + r*dt ! calculate incremental burn distance
  end subroutine

  subroutine calcsurf(id, db, length, od, dt, vol, surf, r)
    ! cylinder burning from id outward and from both ends along the length
    implicit none
    real(dp) id, db, length, od, dt, vol, surf, r

    surf = pi * (id + 2.0d0*db) * (length - 2.0d0*db) + pi * (od**2.0d0 - (id + 2.0*db)**2.0d0) * 0.5

    if (id + 2d0*db .gt. od .or. db.gt.length/2d0) THEN
      surf = 0d0 ! we hit the wall and burned out
      r = 0 ! turn off burn rate so burn distance stops increasing
    endif

    vol = vol + r*surf*dt ! increment the interior volume of the chamber a little
  end subroutine

  subroutine calmdotgen(rhos, r, surf, cp, Tflame, mdotgen, edotgen)
    implicit none
    real(dp) rhos, r, surf, cp, Tflame, mdotgen, edotgen

    mdotgen = rhos * r * surf
    edotgen = mdotgen * cp * Tflame
  end subroutine

  subroutine massflow(p, pamb, area, t, g, rgas, cp, dsigng, mdotos, edotos)
    implicit none

    real(dp) p, pamb, area, t, g, rgas, cp, dsigng, mdotos, edotos

    real(dp) :: mdtx, engyx
    real(dp) :: tx, gx, rx, px, cpx, pcrit, facx, term1, term2, pratio, cstar, ax, hx
    real(dp) :: p1, p2

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

  subroutine addmass(mdotos, edotos, mdotgen, edotgen, dt, mcham, echam)
    implicit none
    real(dp) mdotos, edotos, mdotgen, edotgen, dt, mcham, echam

    mcham = mcham + (mdotgen - mdotos) * dt
    echam = echam + (edotgen - edotos) * dt
  end subroutine

  subroutine calct(echam, mcham, cv, t)
    implicit none
     real(dp) echam, mcham, cv, t

    t = echam / mcham / cv
  end subroutine

  subroutine calcp(mcham, rgas, t, vol, p)
    implicit none
    real(dp) mcham, rgas, t, vol, p

    p = mcham * rgas * t / vol
  end subroutine

  subroutine calcthrust(p, pamb, area, cf, altitude, vel, thrust, den, drag, netthrust)
    implicit none
    real(dp) p, pamb, area, cf, altitude, vel, thrust, den, drag, netthrust

    thrust = (p - pamb) * area * cf ! correction to thrust (actual vs vacuum thrust)
    den = rhob * exp(-gravity * mwair * altitude / RU / tamb)
    drag = -cd * 0.5 * den * vel * abs(vel) * surfrocket
    netthrust=thrust+drag
  end subroutine

  subroutine height (mdotgen, netthrust, rocketmass, mcham, gravity, dt, vel, altitude, propmass, accel)
    implicit none
     real(dp) mdotgen, netthrust, rocketmass, mcham, gravity, dt, vel, altitude, propmass, accel

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

    implicit none

    real(dp), intent(in) :: dt_, t_max_
    real(dp), intent(in) :: c_p_, MW_
    real(dp), intent(in) :: temperature_, pressure_
    real(dp), intent(in) :: T_flame_, r_ref_, n_
    real(dp), intent(in) :: id_, od_, length_, rho_solid_
    real(dp), intent(in) :: dia_, C_f_
    real(dp), allocatable :: rocket(:,:)

    real(dp)  area, Cf, Cp, Cv, dia, dt, echam, g, i
    real(dp)  id, length, mcham, mw, n, nsteps, od, p, Pamb, pref, psipa
    real(dp)  Rgas, rhos, rref, T, Tflame, tmax

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
 
    block 
      real(dp) r, surf, mdotgen, edotgen, dsigng, edotos, den
      real(dp) vol, db, thrust, mdotos, time, propmass, drag, netthrust, accel, vel, altitude, rocketmass
      real(dp), allocatable :: output(:,:)
   
      vol = 1.
      db = 0.
      thrust = 0.
      mdotos = 0.
      time = 0.; propmass = 0.; drag = 0.; netthrust = 0.
      accel = 0.; vel = 0.; altitude = 0.; rocketmass = 0.

      ! propellant burn rate information
      psipa = 6894.76d0 ! pascals per psi (constant)
      pref = 3000d0 * psipa ! reference pressure (constant)

      nsteps=nint(tmax/dt) ! number of time steps

      ! preallocate an output file for simulation infomration
      allocate(output(0:nsteps,11))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

      !! now begin calculating and initializing
      ! gas variables
      rgas = Ru / mw
      cv = cp - rgas
      g = cp / cv

      area = pi / 4d0 * dia**2.0d0 ! nozzle area

      pamb = 101325d0 ! atmospheric pressure

      ! calculate initial mass and energy in the chamber
      mcham = p * vol / rgas / t ! use ideal gas law to determine mass in chamber
      echam = mcham * cv * t ! initial internal energy in chamber

      output(0,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]

      call propwt(od, id, length, rhos, propmass, rocketmass)
      do i=1,nsteps
        call burnrate(rref, p, pref, n, dt, r, db)
        call calcsurf(id, db, length, od, dt, vol, surf, r)
        call calmdotgen(rhos, r, surf, cp, Tflame, mdotgen, edotgen) 
          ! [mdot,engy,dsign] = massflow(p1,pamb,t1,tamb,cp,cp,rgas,rgas,g,g,area)
        call massflow(p, pamb, area, t, g, rgas, cp, dsigng, mdotos, edotos)
        call addmass(mdotos, edotos, mdotgen, edotgen, dt, mcham, echam) 
        call calct(echam, mcham, cv, t)
        call calcp(mcham, rgas, t, vol, p)
        call calcthrust(p, pamb, area, cf, altitude, vel, thrust, den, drag, netthrust)
        call height (mdotgen, netthrust, rocketmass, mcham, gravity, dt, vel, altitude, propmass, accel)
        time = time + dt
        output(i,:) = [time, p, t, mdotos, thrust, drag, netthrust, vol, accel, vel, altitude]
      end do
      rocket = output
    end block
  end function
end module
