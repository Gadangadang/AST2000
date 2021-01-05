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


sigma = 0.0000126*1/5
N = 152000
noise = np.random.normal(0,sigma, size = int(N))
V_rad = [noise[0]]
Vofsystem = 46 # AU/year

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

    N = 152000
    time = 7.6*2 #time for 20 rotations for our home planet
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

        v_inp = np.linalg.norm(v_sun)*x_sun[i+1,1]/np.linalg.norm(x_sun[i+1])
        V_rad.append(v_inp+noise[i+1])
    return x_planet,x_sun, V_rad, dt

Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
system.initial_velocities[0][0],system.initial_velocities[1][0]]

Sun = [0,0,0,0]

planet_orbit, sun_orbit, V_rad, dt = integrator(Planet0,Sun)

V_rad = np.asarray(V_rad) + Vofsystem


infile = open("Radial_vel.txt","w")
for i in range(len(V_rad)):
    infile.write("{} {}\n".format(i*dt,V_rad[i]))
infile.close()
