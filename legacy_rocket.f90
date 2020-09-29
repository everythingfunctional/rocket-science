module mod1
implicit none
save
integer, parameter :: precision=15, range=307
integer, parameter :: dp = selected_real_kind(precision, range)
real(dp), parameter :: gravity=9.81d0
real(dp), parameter :: pi=3.1415926539
real(dp), parameter :: RU=8314d0
real(dp), parameter :: zero=0._dp, one=1._dp
real(dp), parameter :: cd=1.1, rhob=1.225,tamb=300d0
real(dp), parameter :: mwair=28.96,surfrocket=pi/4  
! assuminng a 1.1 drag coefficient and 

real(dp):: cp,cv,g,rgas,mw,vol=one,dia,cf,id,od,length,rref,rhos,psipa,pref
real(dp):: db=zero,dt,tmax,Tflame
real(dp):: thrust=zero, area, r, n, surf,mdotgen,mdotout,edotgen,edotout,energy
real(dp):: mdotos=zero, edotos, texit, dsigng,pamb,p,t
real(dp):: mcham,echam,time=zero,propmass=zero,drag=zero
integer nsteps,i
real(dp):: accel=zero, vel=zero, altitude=zero, rocketmass=zero
real(dp), allocatable :: output(:,:)
real(dp) den ! air density
end module

subroutine propwt ! calculate weight of propellent
  use mod1
  implicit none
  propmass=pi/4*(od**2-id**2)*length*rhos
  rocketmass=0.15*propmass ! assume 85% propellant loading and 15% extra wt of rocket
end subroutine


subroutine burnrate
  use mod1
  implicit none
  r=rref*(p/pref)**n ! calculate burn rate
  db=db+r*dt  ! calculate incremental burn distance
end subroutine

subroutine calcsurf
  ! cylinder burning from id outward and from both ends along the length
  use mod1
  implicit none

  surf=pi*(id+2.0d0*db)*(length-2.0d0*db)+0.5d0*pi*(od**2.0d0-(id+2.0*db)**2.0d0)

  if(id+2d0*db.gt.od.or.db.gt.length/2d0) THEN
     surf=0d0  ! we hit the wall and burned out
     r=0  ! turn off burn rate so burn distance stops increasing
   endif

vol=vol+r*surf*dt ! increment the interior volume of the chamber a little
end subroutine

subroutine calmdotgen
  use mod1
  implicit none
  mdotgen=rhos*r*surf
  edotgen=mdotgen*cp*Tflame
end subroutine

subroutine massflow
   USE mod1
   implicit none
   REAL (8)::mdtx,engyx
   REAL (8)::tx,gx,rx,px,cpx,pcrit,facx,term1,term2,pratio,cstar,ax,hx
   REAL (8):: p1,p2

   mdotos=0.
   edotos=0.  ! initially set them to zero prior to running this loop

     p1=p
     p2=pamb
     ax=area
     IF(p1.GT.p2) THEN
        dsigng=1
        tx=t
        gx=g
        rx=rgas
        px=p
        cpx=cp
        hx=cp*t
        pratio=p1/p2
     else
        dsigng=-1
        tx=tamb
        gx=g
        rx=rgas
        px=pamb
        cpx=cp
        hx=cp*tamb
        pratio=p2/p1
    end if

    pcrit=(2./(gx+1.))**(gx/(gx-1.))
    IF((1./pratio).LT.pcrit) then
        ! choked flow
        cstar=sqrt((1./gx)*((gx+1.)/2.)**((gx+1.)/(gx-1.))*rx*tx)
        mdtx=px*ax/cstar
    else
        ! unchoked flow
      facx=pratio**((gx-1.)/gx)
      term1=SQRT(gx*rx*tx/facx)
      term2=SQRT((facx-1.)/(gx-1.))
      mdtx=SQRT(2.)*px/pratio/rx/tx*facx*term1*term2*ax
    end if
    engyx=mdtx*hx  ! reformulate based on enthalpy of the chamber
    mdotos=mdtx*dsigng ! exiting mass flow (could be negative "dsigng")
    edotos=engyx*dsigng ! exiting enthalpy
end subroutine

subroutine addmass
    use mod1
    implicit none
    mcham=mcham+(mdotgen-mdotos)*dt
    echam=echam+(edotgen-edotos)*dt
end subroutine

subroutine calct
    use mod1
    implicit none
    t=echam/mcham/cv
end subroutine

subroutine calcp
    use mod1
    implicit none
    p=mcham*rgas*t/vol
end subroutine

subroutine calcthrust
    use mod1
    implicit none
    thrust=(p-pamb)*area*cf ! correction to thrust (actual vs vacuum thrust)
    den=rhob*exp(-gravity*mwair*altitude/RU/tamb)
    drag=-cd*0.5*den*vel*abs(vel)*surfrocket

    thrust=thrust+drag
end subroutine

subroutine height
  use mod1
  implicit none
  propmass=propmass-mdotgen*dt ! incremental change in propellant mass
  accel=thrust/(propmass+rocketmass+mcham)-gravity
  vel=vel+accel*dt
  altitude=altitude+vel*dt
end subroutine



!!  Main program


function legacy_rocket(input_file)
  !! this is a basic program of a single stage
  !! rocket motor flowing out of a nozzle, assuming
  !! a thrust coefficient and ignoring the complexities of
  !! what happens to thrust at low pressures, i.e. shock in the nozzle

use mod1
implicit none

character(len=*), intent(in) :: input_file
real(dp), allocatable :: legacy_rocket(:,:)

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

open(newunit=file_unit, file=input_file, status="old")

read(file_unit, nml=numerics_list)
dt   = dt_
tmax = t_max_

read(file_unit, nml=gas_list)
cp = c_p_
mw = MW_

read(file_unit, nml=state_list)
t = temperature_
p = pressure_

read(file_unit, nml=combustion_list)
Tflame = T_flame_
rref   = r_ref_
n      = n_

read(file_unit, nml=grain_list)
id     = id_
od     = od_
length = length_
rhos   = rho_solid_

read(file_unit, nml=nozzle_list)
dia = dia_
cf  = C_f_

close(file_unit)


!  propellent grain is a cylinder burning radially outward and axially inward from one end.
!  the other end is considered inhibited.
! outer diameter is inhibited because this is a cast propellent: it was poured
! into the tube/chamber and only the inner diameter burns when ignited.

  ! propellant burn rate information
  psipa=6894.76d0 ! pascals per psi (constant)
  pref=3000d0*psipa ! reference pressure (constant)

  nsteps=nint(tmax/dt) ! number of time steps

! preallocate an output file for simulation infomration
  allocate(output(0:nsteps,9))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

!! now begin calculating and initializing
! gas variables
  rgas=ru/mw
  cv=cp-rgas
  g=cp/cv

  area=pi/4d0*dia**2.0d0 ! nozzle area

  pamb=101325d0 ! atmospheric pressure

!  calculate initial mass and energy in the chamber
  mcham=p*vol/rgas/t  ! use ideal gas law to determine mass in chamber
  echam=mcham*cv*t ! initial internal energy in chamber

  output(0,:)=[time,p,t,mdotos,thrust,vol, accel, vel, altitude]
  output(0,:)=0.0
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
   time=time+dt
   output(i,:)=[time,p,t,mdotos,thrust,vol,accel,vel,altitude]

  enddo

  legacy_rocket = output

end function legacy_rocket

program main
  use mod1, only: dp
  implicit none

  interface
    function legacy_rocket(input)
      import dp
      character(len=*), intent(in) :: input
      real(dp), allocatable :: legacy_rocket(:,:)
    end function
  end interface

  integer :: i
  real(dp), allocatable :: output(:,:)

  output = legacy_rocket("rocket.inp")
  do i = 0, size(output, 1)
      print *, output(i,:)
  end do
end program
