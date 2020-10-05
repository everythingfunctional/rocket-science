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
  pure function initial_propellant_mass( &
      propellant_inner_diameter, &
      propellant_length, &
      propellant_outer_diameter, &
      propellant_density)
    real(dp), intent(in) :: propellant_inner_diameter
    real(dp), intent(in) :: propellant_length
    real(dp), intent(in) :: propellant_outer_diameter
    real(dp), intent(in) :: propellant_density
    real(dp) :: initial_propellant_mass

    associate( &
        id => propellant_inner_diameter, &
        od => propellant_outer_diameter, &
        l => propellant_length)
      associate(propellant_volume => PI * (od**2 - id**2) / 4.0_dp * l)
        initial_propellant_mass = propellant_volume * propellant_density
      end associate
    end associate
  end function

  pure function initial_rocket_mass(propellant_mass)
    real(dp), intent(in) :: propellant_mass
    real(dp) :: initial_rocket_mass

    initial_rocket_mass = 0.15_dp * propellant_mass
  end function

  pure function calculate_burn_rate(&
      burn_rate_exponent, &
      pressure, &
      reference_burn_rate) &
      result(burn_rate)
    real(dp), intent(in) :: burn_rate_exponent
    real(dp), intent(in) :: pressure
    real(dp), intent(in) :: reference_burn_rate
    real(dp) :: burn_rate

    real(dp), parameter :: PASCALS_PER_PSI = 6894.76_dp
    real(dp), parameter :: REFERENCE_PRESSURE = 3000.0_dp * PASCALS_PER_PSI

    associate( &
        rbr => reference_burn_rate, &
        p => pressure, &
        p_r => REFERENCE_PRESSURE, &
        n => burn_rate_exponent)
      burn_rate = rbr * (p/p_r)**n
    end associate
  end function

  pure function update_burn_depth( &
      previous_burn_depth, &
      burn_rate, &
      time_step_length) &
      result(new_burn_depth)
    real(dp), intent(in) :: previous_burn_depth
    real(dp), intent(in) :: burn_rate
    real(dp), intent(in) :: time_step_length
    real(dp) :: new_burn_depth

    new_burn_depth = previous_burn_depth + burn_rate*time_step_length
  end function

  pure function calculate_burning_surface_area( &
      burn_depth, &
      propellant_inner_diameter, &
      propellant_length, &
      propellant_outer_diameter) &
      result(burning_surface_area)
    real(dp), intent(in) :: burn_depth
    real(dp), intent(in) :: propellant_inner_diameter
    real(dp), intent(in) :: propellant_length
    real(dp), intent(in) :: propellant_outer_diameter
    real(dp) :: burning_surface_area

    ! cylinder burning from inside and both ends
    associate( &
        id => propellant_inner_diameter, &
        od => propellant_outer_diameter, &
        l => propellant_length, &
        d => burn_depth)
      associate( &
          new_l => l - 2.0_dp*d, &
          new_id => id + 2.0_dp*d)
        if (new_id > od .or. new_l < 0.0_dp) then
          burning_surface_area = 0.0_dp
        else
          burning_surface_area = &
              PI * new_id * new_l &
              + PI * (od**2 - new_id**2) / 4.0_dp * 2.0_dp
        end if
      end associate
    end associate
  end function

  pure function recheck_burn_rate( &
      current_burn_rate, &
      burn_depth, &
      propellant_inner_diameter, &
      propellant_length, &
      propellant_outer_diameter) &
      result(burn_rate)
    real(dp), intent(in) :: current_burn_rate
    real(dp), intent(in) :: burn_depth
    real(dp), intent(in) :: propellant_inner_diameter
    real(dp), intent(in) :: propellant_length
    real(dp), intent(in) :: propellant_outer_diameter
    real(dp) :: burn_rate

    associate( &
        id => propellant_inner_diameter, &
        od => propellant_outer_diameter, &
        l => propellant_length, &
        d => burn_depth)
      associate( &
          new_l => l - 2.0_dp*d, &
          new_id => id + 2.0_dp*d)
        if (new_id > od .or. new_l < 0.0_dp) then
          burn_rate = 0.0_dp
        else
          burn_rate = current_burn_rate
        end if
      end associate
    end associate
  end function

  pure function update_chamber_volume( &
      previous_chamber_volume, &
      burn_rate, &
      burning_surface_area, &
      time_step_length) &
      result(chamber_volume)
    real(dp), intent(in) :: previous_chamber_volume
    real(dp), intent(in) :: burn_rate
    real(dp), intent(in) :: burning_surface_area
    real(dp), intent(in) :: time_step_length
    real(dp) :: chamber_volume

    chamber_volume = &
        previous_chamber_volume &
        + burn_rate*burning_surface_area*time_step_length
  end function

  pure function calculate_mass_generation_rate( &
      propellant_density, &
      burn_rate, &
      burning_surface_area) &
      result(mass_generation_rate)
    real(dp), intent(in) :: propellant_density
    real(dp), intent(in) :: burn_rate
    real(dp), intent(in) :: burning_surface_area
    real(dp) :: mass_generation_rate

    mass_generation_rate = propellant_density * burn_rate * burning_surface_area
  end function

  pure function calculate_energy_generation_rate( &
      mass_generation_rate, &
      heat_capacity_at_constant_pressure, &
      flame_temperature) &
      result(energy_generation_rate)
    real(dp), intent(in) :: mass_generation_rate
    real(dp), intent(in) :: heat_capacity_at_constant_pressure
    real(dp), intent(in) :: flame_temperature
    real(dp) :: energy_generation_rate

    energy_generation_rate = &
        mass_generation_rate &
        * heat_capacity_at_constant_pressure &
        * flame_temperature
  end function

  pure function calculate_mass_flow_rate( &
      flow_area, &
      heat_capacity_ratio, &
      chamber_pressure, &
      specific_gas_constant, &
      chamber_temperature) &
      result(mass_outflow_rate)
    real(dp), intent(in) :: flow_area
    real(dp), intent(in) :: heat_capacity_ratio
    real(dp), intent(in) :: chamber_pressure
    real(dp), intent(in) :: specific_gas_constant
    real(dp), intent(in) :: chamber_temperature
    real(dp) :: mass_outflow_rate

    real(dp) :: critical_pressure_ratio
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: f3
    real(dp) :: flow_direction
    real(dp) :: flow_pressure
    real(dp) :: characteristic_velocity
    real(dp) :: flow_temperature
    real(dp) :: mass_flow_rate
    real(dp) :: pressure_ratio

    if (chamber_pressure.GT.ATMOSPHERIC_PRESSURE) then
      flow_direction = 1.0_dp
      flow_temperature = chamber_temperature
      flow_pressure = chamber_pressure
      pressure_ratio = chamber_pressure / ATMOSPHERIC_PRESSURE
    else
      flow_direction = -1.0_dp
      flow_temperature = AMBIENT_TEMPERATURE
      flow_pressure = ATMOSPHERIC_PRESSURE
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
        characteristic_velocity = &
            sqrt( &
                (1.0_dp / g) &
                * ((g + 1.0_dp) / 2.0_dp)**((g + 1.0_dp) / (g - 1.0_dp)) &
                * r * t)
      end associate
      mass_flow_rate = flow_pressure * flow_area / characteristic_velocity
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
    mass_outflow_rate = mass_flow_rate * flow_direction
  end function

  pure function calculate_energy_flow_rate( &
      heat_capacity_at_constant_pressure, &
      chamber_pressure, &
      chamber_temperature, &
      mass_outflow_rate) &
      result(energy_outflow_rate)
    real(dp), intent(in) :: heat_capacity_at_constant_pressure
    real(dp), intent(in) :: chamber_pressure
    real(dp), intent(in) :: chamber_temperature
    real(dp), intent(in) :: mass_outflow_rate
    real(dp) :: energy_outflow_rate

    real(dp) :: flow_enthalpy

    if (chamber_pressure.GT.ATMOSPHERIC_PRESSURE) then
      flow_enthalpy = heat_capacity_at_constant_pressure * chamber_temperature
    else
      flow_enthalpy = heat_capacity_at_constant_pressure * AMBIENT_TEMPERATURE
    end if

    energy_outflow_rate = mass_outflow_rate * flow_enthalpy
  end function

  subroutine update_chamber_contents( &
      time_step_length, &
      energy_generation_rate, &
      energy_outflow_rate, &
      mass_generation_rate, &
      mass_outflow_rate, &

      chamber_energy, &
      chamber_mass)
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: energy_generation_rate
    real(dp), intent(in) :: energy_outflow_rate
    real(dp), intent(in) :: mass_generation_rate
    real(dp), intent(in) :: mass_outflow_rate
    real(dp), intent(inout) :: chamber_energy
    real(dp), intent(inout) :: chamber_mass

    chamber_mass = &
        chamber_mass &
        + (mass_generation_rate - mass_outflow_rate) * time_step_length
    chamber_energy = &
        chamber_energy &
        + (energy_generation_rate - energy_outflow_rate) * time_step_length
  end subroutine

  subroutine calculate_temperature( &
      heat_capacity_at_constant_volume, &
      chamber_energy, &
      chamber_mass, &

      chamber_temperature)
    real(dp), intent(in) :: heat_capacity_at_constant_volume
    real(dp), intent(in) :: chamber_energy
    real(dp), intent(in) :: chamber_mass
    real(dp), intent(out) :: chamber_temperature

    chamber_temperature = &
        chamber_energy / chamber_mass / heat_capacity_at_constant_volume
  end subroutine

  subroutine calculate_pressure( &
      chamber_mass, &
      specific_gas_constant, &
      chamber_temperature, &
      chamber_volume, &

      chamber_pressure)
    real(dp), intent(in) :: chamber_mass
    real(dp), intent(in) :: specific_gas_constant
    real(dp), intent(in) :: chamber_temperature
    real(dp), intent(in) :: chamber_volume
    real(dp), intent(out) :: chamber_pressure

    chamber_pressure = &
        chamber_mass * specific_gas_constant * chamber_temperature / chamber_volume
  end subroutine

  subroutine calculate_thrust( &
      altitude, &
      flow_area, &
      thrust_correction_factor, &
      chamber_pressure, &
      velocity, &

      drag, &
      net_thrust, &
      thrust)
    real(dp), intent(in) :: altitude
    real(dp), intent(in) :: flow_area
    real(dp), intent(in) :: thrust_correction_factor
    real(dp), intent(in) :: chamber_pressure
    real(dp), intent(in) :: velocity
    real(dp), intent(out) :: drag
    real(dp), intent(out) :: net_thrust
    real(dp), intent(out) :: thrust

    real(dp), parameter :: DRAG_COEFFICIENT = 1.1_dp
    real(dp), parameter :: MOLECULAR_WEIGHT_OF_AIR = 28.96_dp
    real(dp), parameter :: REFERENCE_AIR_DENSITY = 1.225_dp
    real(dp), parameter :: ROCKET_SURFACE_AREA = PI / 4.0_dp
    real(dp) :: air_density

    thrust = &
        (chamber_pressure - ATMOSPHERIC_PRESSURE) &
        * flow_area * thrust_correction_factor
    air_density = &
        REFERENCE_AIR_DENSITY &
        * exp( &
            -GRAVITY &
            * MOLECULAR_WEIGHT_OF_AIR &
            * altitude &
            / UNIVERSAL_GAS_CONSTANT &
            / AMBIENT_TEMPERATURE)
    drag = &
        -DRAG_COEFFICIENT &
        * 0.5_dp &
        * air_density &
        * velocity &
        * abs(velocity) &
        * ROCKET_SURFACE_AREA
    net_thrust = thrust + drag
  end subroutine

  subroutine update_trajectory( &
      time_step_length, &
      chamber_mass, &
      mass_generation_rate, &
      net_thrust, &
      rocket_mass, &

      altitude, &
      propellant_mass, &
      velocity, &

      acceleration)
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: chamber_mass
    real(dp), intent(in) :: mass_generation_rate
    real(dp), intent(in) :: net_thrust
    real(dp), intent(in) :: rocket_mass
    real(dp), intent(inout) :: altitude
    real(dp), intent(inout) :: propellant_mass
    real(dp), intent(inout) :: velocity
    real(dp), intent(out) :: acceleration

    propellant_mass = propellant_mass - mass_generation_rate*time_step_length
    acceleration = net_thrust / (propellant_mass + rocket_mass + chamber_mass) - GRAVITY
    velocity = velocity + acceleration*time_step_length
    altitude = altitude + velocity*time_step_length
  end subroutine

  function rocket( &
      time_step_length, &
      simulation_end_time, &
      heat_capacity_at_constant_pressure, &
      molecular_weight, &
      initial_temperature, &
      initial_pressure, &
      flame_temperature, &
      reference_burn_rate, &
      burn_rate_exponent, &
      propellant_inner_diameter, &
      propellant_outer_diameter, &
      propellant_length, &
      propellant_density, &
      nozzle_diameter, &
      thrust_correction_factor)
    !! this is a basic program of a single stage
    !! rocket motor flowing out of a nozzle, assuming
    !! a thrust coefficient and ignoring the complexities of
    !! what happens to thrust at low pressures, i.e. shock in the nozzle
    real(dp), intent(in) :: time_step_length
    real(dp), intent(in) :: simulation_end_time
    real(dp), intent(in) :: heat_capacity_at_constant_pressure
    real(dp), intent(in) :: molecular_weight
    real(dp), intent(in) :: initial_temperature
    real(dp), intent(in) :: initial_pressure
    real(dp), intent(in) :: flame_temperature
    real(dp), intent(in) :: reference_burn_rate
    real(dp), intent(in) :: burn_rate_exponent
    real(dp), intent(in) :: propellant_inner_diameter
    real(dp), intent(in) :: propellant_outer_diameter
    real(dp), intent(in) :: propellant_length
    real(dp), intent(in) :: propellant_density
    real(dp), intent(in) :: nozzle_diameter
    real(dp), intent(in) :: thrust_correction_factor
    real(dp), allocatable :: rocket(:,:)

    real(dp), parameter :: ONE = 1.0_dp
    real(dp), parameter :: ZERO = 0.0_dp
    real(dp) :: acceleration = ZERO
    real(dp) :: altitude = ZERO
    real(dp) :: burn_depth = ZERO
    real(dp) :: burn_rate
    real(dp) :: burning_surface_area
    real(dp) :: chamber_energy
    real(dp) :: chamber_mass
    real(dp) :: chamber_pressure
    real(dp) :: chamber_temperature
    real(dp) :: chamber_volume = ONE
    real(dp) :: drag = ZERO
    real(dp) :: energy_generation_rate
    real(dp) :: energy_outflow_rate
    real(dp) :: flow_area
    real(dp) :: propellant_mass = ZERO
    real(dp) :: heat_capacity_at_constant_volume
    real(dp) :: heat_capacity_ratio
    integer :: i
    real(dp) :: mass_generation_rate
    real(dp) :: mass_outflow_rate = ZERO
    real(dp) :: net_thrust = ZERO
    integer :: num_time_steps
    real(dp) :: rocket_mass = ZERO
    real(dp) :: specific_gas_constant
    real(dp) :: thrust = ZERO
    real(dp) :: time = ZERO
    real(dp) :: velocity = ZERO

    chamber_temperature = initial_temperature
    chamber_pressure = initial_pressure
    num_time_steps = nint(simulation_end_time / time_step_length)
    allocate(rocket(0:num_time_steps, 11))
    specific_gas_constant = UNIVERSAL_GAS_CONSTANT / molecular_weight
    heat_capacity_at_constant_volume = &
        heat_capacity_at_constant_pressure - specific_gas_constant
    heat_capacity_ratio = &
        heat_capacity_at_constant_pressure / heat_capacity_at_constant_volume
    flow_area = PI / 4.0_dp * nozzle_diameter**2
    chamber_mass = &
        chamber_pressure &
        * chamber_volume &
        / specific_gas_constant &
        / chamber_temperature
    chamber_energy = &
        chamber_mass * heat_capacity_at_constant_volume * chamber_temperature
    rocket(0,:) = [ &
        time, &
        chamber_pressure, &
        chamber_temperature, &
        mass_outflow_rate, &
        thrust, &
        drag, &
        net_thrust, &
        chamber_volume, &
        acceleration, &
        velocity, &
        altitude]

    propellant_mass = initial_propellant_mass( &
        propellant_inner_diameter, &
        propellant_length, &
        propellant_outer_diameter, &
        propellant_density)
    rocket_mass = initial_rocket_mass(propellant_mass)
    do i = 1, num_time_steps
      burn_rate = calculate_burn_rate( &
          burn_rate_exponent, chamber_pressure, reference_burn_rate)
      burn_depth = update_burn_depth(burn_depth, burn_rate, time_step_length)
      burning_surface_area = calculate_burning_surface_area( &
          burn_depth, &
          propellant_inner_diameter, &
          propellant_length, &
          propellant_outer_diameter)
      burn_rate = recheck_burn_rate( &
          burn_rate, &
          burn_depth, &
          propellant_inner_diameter, &
          propellant_length, &
          propellant_outer_diameter)
      chamber_volume = update_chamber_volume( &
          chamber_volume, burn_rate, burning_surface_area, time_step_length)
      mass_generation_rate = calculate_mass_generation_rate( &
          propellant_density, burn_rate, burning_surface_area)
      energy_generation_rate = calculate_energy_generation_rate( &
          mass_generation_rate, &
          heat_capacity_at_constant_pressure, &
          flame_temperature)
      mass_outflow_rate = calculate_mass_flow_rate( &
          flow_area, &
          heat_capacity_ratio, &
          chamber_pressure, &
          specific_gas_constant, &
          chamber_temperature)
      energy_outflow_rate = calculate_energy_flow_rate( &
          heat_capacity_at_constant_pressure, &
          chamber_pressure, &
          chamber_temperature, &
          mass_outflow_rate)
      call update_chamber_contents( &
          time_step_length, &
          energy_generation_rate, &
          energy_outflow_rate, &
          mass_generation_rate, &
          mass_outflow_rate, &

          chamber_energy, &
          chamber_mass)
      call calculate_temperature( &
          heat_capacity_at_constant_volume, &
          chamber_energy, &
          chamber_mass, &

          chamber_temperature)
      call calculate_pressure( &
          chamber_mass, &
          specific_gas_constant, &
          chamber_temperature, &
          chamber_volume, &

          chamber_pressure)
      call calculate_thrust( &
          altitude, &
          flow_area, &
          thrust_correction_factor, &
          chamber_pressure, &
          velocity, &

          drag, &
          net_thrust, &
          thrust)
      call update_trajectory( &
          time_step_length, &
          chamber_mass, &
          mass_generation_rate, &
          net_thrust, &
          rocket_mass, &

          altitude, &
          propellant_mass, &
          velocity, &

          acceleration)
      time = time + time_step_length
      rocket(i,:) = [ &
          time, &
          chamber_pressure, &
          chamber_temperature, &
          mass_outflow_rate, &
          thrust, &
          drag, &
          net_thrust, &
          chamber_volume, &
          acceleration, &
          velocity, &
          altitude]
    end do
  end function
end module
