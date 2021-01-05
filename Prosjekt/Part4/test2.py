import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts




utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
codes = [8208] # The codes recieved from the group teacher (in this case just one)
shortcuts = SpaceMissionShortcuts(mission, codes)

shortcuts.place_spacecraft_on_escape_trajectory(0.4,30000,pi/24, 3000)

#Interpolate a coordinate at a given time
def interpol(a,b,t):
    t = t
    f = interpolate.interp1d(a,b,axis=0)
    return f(t)

#Find position of each object, with given radii of circle for the reference points
def find(x,y):

    x1,y1 = (system.initial_positions[0][3],system.initial_positions[1][3])
    x2,y2 = (system.initial_positions[0][4],system.initial_positions[1][4])
    x3,y3 = (0,0)
    r1 = np.sqrt((x-x1)**2+(y-y1)**2)
    r2 = np.sqrt((x-x2)**2+(y-y2)**2)
    r3 = np.sqrt((x-x3)**2+(y-y3)**2)
    print("r1 = {:g}, r2 = {:g}, r3 = {:g}".format(r1,r2,r3))
    return x1,y1,x2,y2,x3,y3,r1,r2,r3

#Trilaterate position
def position(x,y):

    x1,y1,x2,y2,x3,y3,r1,r2,r3 = find(x,y)
    A = 2*x2 - 2*x1
    B = 2*y2 - 2*y1
    C = r1**2 - r2**2 - x1**2 + x2**2 - \
    y1**2 + y2**2

    D = 2*x3 - 2*x2
    E = 2*y3 - 2*y2
    F = r2**2 - r3**2 - x2**2 + x3**2 - \
    y2**2 + y3**2

    x = (C*E - F*B) / (E*A - B*D)
    y = (C*D - A*F) / (B*D - A*E)


    return x,y,x1,y1,x2,y2,x3,y3,r1,r2,r3

"""Loops through each planet and makes array with neccesary info"""
Planets = np.zeros((8,5))
for i in range(len(Planets)):
    Planet0_x0 = [system.initial_positions[0,i],system.initial_positions[1,i]]
    Planet0_v0 = [system.initial_velocities[0][i],system.initial_velocities[1][i]]
    Planets[i,0] = Planet0_x0[0]
    Planets[i,1] = Planet0_x0[1]
    Planets[i,2] = Planet0_v0[0]
    Planets[i,3] = Planet0_v0[1]
    Planets[i,4] = i

Gr = G_sol
sunmass = system.star_mass

"""Numerical integrator that calculates orbit with leapfrog method"""
def integrator(Planet):
    def gravity(x,y):
        grav = -(Gr*sunmass)/(np.sqrt(x**2+y**2))**2
        ax = grav*x /(np.sqrt(x**2+y**2))
        ay = grav*y/(np.sqrt(x**2+y**2))
        return np.asarray([ax,ay])
    x0x,x0y,v0x,v0y,planet_index = Planet

    #N = 1400000

    time = 7.6*40 #time for 20 rotations for our home planet

    N=1400000
    dt = time/N
    t = dt
    x = np.zeros((N,2),float)
    x[0,:] = [x0x,x0y]
    v = np.asarray([v0x,v0y])
    a_i = gravity(x[0,0],x[0,1])
    for i in range(N-1):
        t+=dt
        x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
        a_iplus1 = gravity(x[i+1,0],x[i+1,1])
        v += 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
    return x, planet_index


"""
Generating lists of objects with their positions
"""
planetorb = []

N=1400000
sunpos = np.zeros((N,2),float)
time = 2
p=0

planet_orbit,planet_index = integrator(Planets[3])
planetorb.append(planet_orbit)
plt.plot(planet_orbit[:,0],planet_orbit[:,1],label="Planet {} orbit".format(planet_index))


planet_orbit,planet_index = integrator(Planets[4])
planetorb.append(planet_orbit)
plt.plot(planet_orbit[:,0],planet_orbit[:,1],label="Planet {} orbit".format(planet_index))


planetorb.append(sunpos)

planetorb = np.asarray(planetorb)
time = 7.6*40 #time for 20 rotations for our home planet
N=1400000
dt = time/N
t = np.zeros((N),float)
t[0] = dt
for i in range(N-1):
    t[i+1] = t[i] + dt

#Find position
tidspunkt = 0.4
#test values
x_test=6
y_test=6


x,y,x1,y1,x2,y2,x3,y3,r1,r2,r3 = position(x_test,y_test)
theta = np.linspace(0, 2*np.pi, 100)
print("X: {} Y:{} Xerror: {} Yerror: {}".format(x,y,abs(x-x_test),abs(y-y_test)))

s1 = x1 + r1*np.cos(theta)
s2 = y1 + r1*np.sin(theta)
s3 = x2 + r2*np.cos(theta)
s4 = y2 + r2*np.sin(theta)
s5 = x3 + r3*np.cos(theta)
s6 = y3 + r3*np.sin(theta)


plt.plot(s1,s2,"b--")
plt.plot(s3,s4,"r--")
plt.plot(s5,s6,"g--")
plt.plot(x,y,"yd", label="Pos [{:g}, {:g}]".format(x,y))
plt.plot(x_test,y_test,"rx", label="Orig pos [{:g}, {:g}]".format(x_test,y_test))
plt.plot(x1,y1,"bo")
plt.plot(x2,y2,"go")
plt.plot(x3,y3,"ro",label="Sun")
plt.plot([x,x1],[y,y1],"k--")
plt.plot([x,x2],[y,y2],"k--")
plt.plot([x,x3],[y,y3],"k--")

plt.axis("equal")
plt.xlabel("X in AU")
plt.ylabel("Y in AU")
plt.title("PathFinder Measuments")
plt.legend()
plt.savefig("Trilateration_simulated_test_3_4.jpeg")
plt.show()
