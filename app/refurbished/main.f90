program main
  use refurbished, only: rocket, dp
  implicit none

  integer :: i
  real(dp), allocatable :: output(:,:)

  integer file_unit

  real(dp) dt_, t_max_
  real(dp) c_p_, MW_
  real(dp) temperature_, pressure_
  real(dp) T_flame_, r_ref_, n_
  real(dp) id_, od_, length_, rho_solid_
  real(dp) dia_, C_f_

  namelist/numerics_list/ dt_, t_max_
  namelist/gas_list/ c_p_, MW_
  namelist/state_list/  temperature_, pressure_
  namelist/combustion_list/ T_flame_, r_ref_, n_
  namelist/grain_list/ id_, od_, length_, rho_solid_
  namelist/nozzle_list/ dia_, C_f_

  open(newunit=file_unit, file="rocket.inp", status="old")

  read(file_unit, nml=numerics_list)
  read(file_unit, nml=gas_list)
  read(file_unit, nml=state_list)
  read(file_unit, nml=combustion_list)
  read(file_unit, nml=grain_list)
  read(file_unit, nml=nozzle_list)

  close(file_unit)

  output = rocket( &
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
  do i = lbound(output, 1), ubound(output, 1)
      print *, output(i,:)
  end do
end program
