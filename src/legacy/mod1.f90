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
real(dp):: mcham,echam,time=zero,propmass=zero,drag=zero,netthrust=zero
integer nsteps,i
real(dp):: accel=zero, vel=zero, altitude=zero, rocketmass=zero
real(dp), allocatable :: output(:,:)
real(dp) den ! air density
end module
