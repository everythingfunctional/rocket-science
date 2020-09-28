module mod1
implicit none
save
real(8):: cp,cv,g,rgas,mw,vol,dia,cf,id,od,length,rref,rhos,psipa,pref,db,dt,tmax
real(8):: thrust, area, r, surf,mdotgen,mdotout,edotgen,edotout,energy
real(8):: pi, mdotos, edotos, texit, dsigng,pamb,p,t
real(8):: mcham,echam,time
real(8):: accel, vel, altitude, rocketmass, propmass, gravity
integer nsteps
real(8), allocatable :: output(:,:)
end module
