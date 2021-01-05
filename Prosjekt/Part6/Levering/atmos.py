import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
#import part_3
#import test_part6

utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

"""

"""
def t0(avg_distanse):
    temp = (system.star_temperature**4/4*((system.star_radius*1000/AU)/(avg_distanse))**2)**(1/4)
    return temp

avg_distanse_planet1 = np.sqrt(system.initial_positions[0][1]**2+system.initial_positions[1][1])
t0 = t0(avg_distanse_planet1)


gamma = 1.4
R = system.radii[1]*1000
K = k_B
g = 11
m_H = m_p
mu = (18.01528/(N_A*1000))/(3*m_H) + (44.01/(N_A*1000))/(3*m_H) + (16/(N_A*1000))/(3*m_H)
rho_R = system.atmospheric_densities[1]
print("Pressure at sealevel is {}".format(rho_R*K*t0/(mu*m_H)))



def T(r):
    temp = t0 - g*(r-R)*mu*m_H/K*(gamma-1)/gamma
    if temp >= t0/2:#Adiabatic condition
        return temp
    else:#Isothermal condition
        return t0/2


def rho_adia(r):
    rho = (rho_R**(gamma-1)-g*mu*m_H*(r-R)/(K*t0)*(gamma-1)/gamma*rho_R**(gamma-1))**(1/(gamma-1))
    return rho


def rho_iso(r,rho_m,rm):
    rho = rho_m*np.exp((-2*g*mu*m_H*(r-rm))/(K*t0))
    return rho


r = np.linspace(system.radii[1]*1000,system.radii[1]*1000+15000000,15000001)
densities = np.zeros(15000001)
temperatures = np.zeros(15000001)
rho_func = lambda r: rho_adia(r)


p = True
for i,j in enumerate(r):
    temp_adia = t0 - g*(j-R)*mu*m_H/K*(gamma-1)/gamma
    if temp_adia < t0/2 and p:
        rho_m = rho_adia(j)
        rm = j
        rho_func = lambda r: rho_iso(r,rho_m,rm)
        p=False

    densities[i] = rho_func(j)
    temperatures[i] = T(j)

densities = np.asarray(densities)
np.save("densities.npy",densities, allow_pickle=True)

temperatures = np.asarray(temperatures)

c_temp = temperatures.copy()
c_temp = c_temp - 274.15

plt.plot(r/1000,densities,label="Density")
plt.legend()
plt.xlabel("Distanse from surface of planet[km]")
plt.ylabel("Density [kg/m^3]")
plt.title("Density of planet")
plt.savefig("Atmosphere_prifile_density.jpeg")
plt.show()

plt.plot(r/1000,temperatures, label="T[k]")
plt.plot(r/1000,c_temp, label="T[C]")
plt.legend()
plt.xlabel("Distanse from surface of planet")
plt.ylabel("Temperature")
plt.title("Temperature of planet")
plt.savefig("Atmosphere_profile_temp.jpeg")
plt.show()
