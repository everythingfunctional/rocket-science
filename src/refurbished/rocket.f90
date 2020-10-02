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
      fuel_density, &

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

  subroutine update_burn_rate_and_depth( &
      time_step_length, &
      burn_rate_exponent, &
      pressure, &
      reference_burn_rate, &

      burn_depth, &

      burn_rate)
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: burn_rate_exponent
    real(dp), intent(in) :: pressure
    real(dp), intent(in) :: reference_burn_rate
    real(dp), intent(inout) :: burn_depth
    real(dp), intent(out) :: burn_rate

    real(dp), parameter :: PASCALS_PER_PSI = 6894.76_dp
    real(dp), parameter :: REFERENCE_PRESSURE = 3000.0_dp * PASCALS_PER_PSI

    associate( &
        rbr => reference_burn_rate, &
        p => pressure, &
        p_r => REFERENCE_PRESSURE, &
        n => burn_rate_exponent)
      burn_rate = rbr * (p/p_r)**n
    end associate
    burn_depth = burn_depth + burn_rate*time_step_length
  end subroutine

  subroutine update_combustion_progress( &
      burn_depth, &
      time_step_length, &
      fuel_inner_diameter, &
      fuel_length, &
      fuel_outer_diameter, &

      burn_rate, &
      chamber_volume, &

      burning_surface_area)
    real(dp), intent(in) :: burn_depth
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: fuel_inner_diameter
    real(dp), intent(in) :: fuel_length
    real(dp), intent(in) :: fuel_outer_diameter
    real(dp), intent(inout) :: burn_rate
    real(dp), intent(inout) :: chamber_volume
    real(dp), intent(out) :: burning_surface_area

    associate( &
        id => fuel_inner_diameter, &
        od => fuel_outer_diameter, &
        l => fuel_length, &
        d => burn_depth)
      associate( &
          new_l => l - 2.0_dp*d, &
          new_id => id + 2.0_dp*d)
        ! cylinder burning from inside and both ends
        burning_surface_area = &
            PI * new_id * new_l &
            + PI * (od**2 - new_id**2) / 4.0_dp * 2.0_dp
        if (new_id > od .or. new_l < 0.0_dp) then
          burning_surface_area = 0.0_dp ! we hit the wall and burned out
          burn_rate = 0.0_dp ! turn off burn rate so burn distance stops increasing
        end if
      end associate
    end associate

    ! increment the interior volume of the chamber
    chamber_volume = &
        chamber_volume &
        + burn_rate*burning_surface_area*time_step_length
  end subroutine

  subroutine calculate_generation_rates( &
      heat_capacity_at_constant_pressure, &
      burn_rate, &
      fuel_density, &
      burning_surface_area, &
      flame_temperature, &

      energy_generation_rate, &
      mass_generation_rate)
    real(dp), intent(in) :: heat_capacity_at_constant_pressure
    real(dp), intent(in) :: burn_rate
    real(dp), intent(in) :: fuel_density
    real(dp), intent(in) :: burning_surface_area
    real(dp), intent(in) :: flame_temperature
    real(dp), intent(out) :: energy_generation_rate
    real(dp), intent(out) :: mass_generation_rate

    mass_generation_rate = fuel_density * burn_rate * burning_surface_area
    energy_generation_rate = &
        mass_generation_rate * heat_capacity_at_constant_pressure * flame_temperature
  end subroutine

  subroutine calculate_flow_rates( &
      flow_area, &
      heat_capacity_at_constant_pressure, &
      heat_capacity_ratio, &
      chamber_pressure, &
      specific_gas_constant, &
      chamber_temperature, &

      energy_outflow_rate, &
      mass_outflow_rate)
    real(dp), intent(in) :: flow_area
    real(dp), intent(in) :: heat_capacity_at_constant_pressure
    real(dp), intent(in) :: heat_capacity_ratio
    real(dp), intent(in) :: chamber_pressure
    real(dp), intent(in) :: specific_gas_constant
    real(dp), intent(in) :: chamber_temperature
    real(dp), intent(out) :: energy_outflow_rate
    real(dp), intent(out) :: mass_outflow_rate

    real(dp) :: critical_pressure_ratio
    real(dp) :: energy_flow_rate
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: f3
    real(dp) :: flow_direction
    real(dp) :: flow_enthalpy
    real(dp) :: flow_pressure
    real(dp) :: flow_speed
    real(dp) :: flow_temperature
    real(dp) :: mass_flow_rate
    real(dp) :: pressure_ratio

    mass_outflow_rate = 0.0_dp
    energy_outflow_rate = 0.0_dp ! initially set them to zero prior to running this loop

    if (chamber_pressure.GT.ATMOSPHERIC_PRESSURE) then
      flow_direction = 1.0_dp
      flow_temperature = chamber_temperature
      flow_pressure = chamber_pressure
      flow_enthalpy = heat_capacity_at_constant_pressure * chamber_temperature
      pressure_ratio = chamber_pressure / ATMOSPHERIC_PRESSURE
    else
      flow_direction = -1.0_dp
      flow_temperature = AMBIENT_TEMPERATURE
      flow_pressure = ATMOSPHERIC_PRESSURE
      flow_enthalpy = heat_capacity_at_constant_pressure * AMBIENT_TEMPERATURE
      pressure_ratio = ATMOSPHERIC_PRESSURE / chamber_pressure
    end if

    associate( &
        cp => critical_pressure_ratio, &
        g => heat_capacity_ratio)
      cp = (2.0_dp / (g + 1.0_dp))**(g / (g - 1.0_dp))
    end associate
    if ((1.0_dp / pressure_ratio) < critical_pressure_ratio) then
      ! choked flow
      associate( &
          g => heat_capacity_ratio, &
          r => specific_gas_constant, &
          t => flow_temperature)
        flow_speed = &
            sqrt( &
                (1.0_dp / g) &
                * ((g + 1.0_dp) / 2.0_dp)**((g + 1.0_dp) / (g - 1.0_dp)) &
                * r * t)
      end associate
      mass_flow_rate = flow_pressure * flow_area / flow_speed
    else
      ! unchoked flow
      f1 = pressure_ratio**((heat_capacity_ratio - 1.0_dp) / heat_capacity_ratio)
      f2 = sqrt(heat_capacity_ratio * specific_gas_constant * flow_temperature / f1)
      f3 = sqrt((f1 - 1.0_dp) / (heat_capacity_ratio - 1.0_dp))
      mass_flow_rate = &
          sqrt(2.0_dp) &
          * flow_pressure &
          / pressure_ratio &
          / specific_gas_constant &
          / flow_temperature &
          * f1 &
          * f2 &
          * f3 &
          * flow_area
    end if
    energy_flow_rate = mass_flow_rate * flow_enthalpy
    mass_outflow_rate = mass_flow_rate * flow_direction
    energy_outflow_rate = energy_flow_rate * flow_direction
  end subroutine

  subroutine update_chamber_contents( &
      time_step_length, &
      energy_generation_rate, &
      energy_outflow_rate, &
      mass_generation_rate, &
      mass_outflow_rate, &

      energy, &
      mass)
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: energy_generation_rate
    real(dp), intent(in) :: energy_outflow_rate
    real(dp), intent(in) :: mass_generation_rate
    real(dp), intent(in) :: mass_outflow_rate
    real(dp), intent(inout) :: energy
    real(dp), intent(inout) :: mass

    mass = mass + (mass_generation_rate - mass_outflow_rate) * time_step_length
    energy = energy + (energy_generation_rate - energy_outflow_rate) * time_step_length
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
      call update_burn_rate_and_depth(dt, n, p, rref,  db,  r)
      call update_combustion_progress(db, dt, id, length, od, r, vol, surf)
      call calculate_generation_rates( &
          cp, r, rhos, surf, Tflame,  edotgen, mdotgen)
      call calculate_flow_rates(area, cp, g, p, rgas, t,  edotos, mdotos)
      call update_chamber_contents( &
          dt, edotgen, edotos, mdotgen, mdotos,  echam, mcham)
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
