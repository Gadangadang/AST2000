import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem


utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

sunmass = system.star_mass
planetmass = system.masses[0]
my = sunmass*planetmass/(sunmass+planetmass)
masstot = sunmass+planetmass

Gr = G_sol
def integrator(planet,sun):

    def gravity_on_sun(x_sun,y_sun,x_planet,y_planet):

        r = np.asarray([x_sun-x_planet, y_sun-y_planet])
        grav = -(Gr*sunmass*planetmass)/(np.linalg.norm(r))**2
        ax = grav*r[0] /(np.linalg.norm(r))
        ay = grav*r[1]/(np.linalg.norm(r))

        return np.asarray([ax,ay],float)

    x0x,x0y,v0x,v0y = planet
    x0xsun,x0ysun,v0xsun,v0ysun = sun

    N = 1000000
    time = 0.02 #time for 20 rotations for our home planet
    #print(time)
    dt = time/N
    t = dt

    x_planet = np.zeros((N,2),float)
    x_sun = np.zeros((N,2),float)

    cm = np.array([0,0])

    x_planet[0,:] = [x0x,x0y]
    x_sun[0,:] = [x0xsun,x0ysun]

    v_planet = np.asarray([v0x,v0y],float)
    v_sun = np.asarray([v0xsun,v0ysun],float)
    F_sun0 = gravity_on_sun(x0xsun,x0ysun,x0x,x0y)
    a_i_sun = F_sun0/sunmass
    a_i_planet = -F_sun0/planetmass
    for i in range(N-1):
        t+=dt
        cm = np.asarray([(sunmass*x_sun[i,0] + planetmass*\
        x_planet[i,0])/(sunmass+planetmass),\
        (sunmass*x_sun[i,1] + planetmass*\
        x_planet[i,1])/(sunmass+planetmass)])

        x_planet[i+1,:] = x_planet[i,:]+(v_planet*dt + 0.5*a_i_planet*dt**2)-cm
        x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm
        F_sun = gravity_on_sun(x_sun[i+1,0],x_sun[i+1,1],x_planet[i+1,0],\
        x_planet[i+1,1])
        a_iplus1_sun = F_sun/sunmass
        a_iplus1_pl = -F_sun/planetmass
        v_planet += 0.5*( a_i_planet + a_iplus1_pl)*dt
        a_i_planet = a_iplus1_pl
        v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
        a_i_sun = a_iplus1_sun

    return x_planet,x_sun , dt, N

Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
system.initial_velocities[0][0],system.initial_velocities[1][0]]

Sun = [0,0,0,0]

planet_orbit, sun_orbit, dt, N = integrator(Planet0,Sun)


"""Defining the intervalls for when the sun is partly or fully infront of star """
Upper =  system.radii[0]*1000/AU + system.star_radius*1000/AU
x_planet = np.copy(planet_orbit)
dt = dt*yr/3600
"""
Creating bolenean arrays with true for the timesteps where planet is in intervall
"""
condition1 = np.logical_and(np.greater_equal(x_planet[:,1],0),np.greater_equal(x_planet[:,0],0))
condition2 = np.logical_and(np.less_equal(x_planet[:,1],Upper),np.greater_equal(x_planet[:,0],0))

dt_in_intervall1 = np.logical_and(condition1,condition2).astype(np.float32)
time_in_int1 = sum(dt_in_intervall1)*dt

""" Defining the intervall where the planet is only partly infront of star"""
Upper = system.star_radius*1000/AU + system.radii[0]*1000/AU
Under = system.star_radius*1000/AU - system.radii[0]*1000/AU

condition3 = np.logical_and(np.greater_equal(x_planet[:,1],Under),np.greater_equal(x_planet[:,0],0))
condition4 = np.logical_and(np.less_equal(x_planet[:,1],Upper),np.greater_equal(x_planet[:,0],0))

dt_in_intervall2 = np.logical_and(condition3,condition4).astype(np.float32)
time_in_int2 = sum(dt_in_intervall2*dt)
"""
dt_in_intervall3 will be only where the planet is fully infront of star
"""
dt_in_intervall3 = dt_in_intervall1-dt_in_intervall2

Areal_sun = 2*np.pi*(system.star_radius*1000/AU)
Areal_planet = 2*np.pi*(system.radii[0]*1000/AU)

Flux_drop = (Areal_planet/Areal_sun)
v = (Areal_planet/Areal_sun)/time_in_int2

time = np.linspace(0,N*dt,N)
curve = np.zeros(N)
"""
Creating baseline for when planet is fully infront of star
"""
curve = np.where(dt_in_intervall1 == 1, dt_in_intervall1-Flux_drop, dt_in_intervall1)
"""
Creating reverse array to mirror first array.
"""
curve1 = curve[::-1]
dt_in_intervall2_2 = dt_in_intervall2[::-1]

for i in range(len(dt_in_intervall2)):
    if dt_in_intervall2[i] == 1:
        curve[i] = curve[i-1] + v*dt
"""
Making sure that when planet is not infront of star, flux=1
"""
curve =  np.where(curve==0,curve+1,curve)
curve1 =  np.where(curve1==0,curve1+1,curve1)
for i in range(len(dt_in_intervall2)):
    if dt_in_intervall2_2[i] == 1:
        curve1[i] = curve1[i-1] - v*dt


curve = np.append(curve1,curve)
time = np.linspace(0,2*N*dt,2*N)

sigma = (Areal_planet/Areal_sun)*0.2
noise = np.random.normal(0,sigma, size = int(2*N))
curve +=noise

"""
Length of array that is plottet was choosen basen on apperence of plot.
"""
"""
This is the write code
"""
"""infile = open("Light_curve_data.txt", "w")
for i in range(len(curve[int(140/dt):int(210/dt)])):
    infile.write("{} {}\n".format(time[int(140/dt):int(210/dt)][i],\
    curve[int(140/dt):int(210/dt)][i]))
infile.close()"""
