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


#print(system.print_info())

r = system.star_radius*1000
T= system.star_temperature
L = sigma*T**4*4*pi*r**2
lum_sun = 3.828E26
print("L star = {}, L/Lsun = {}".format(L, L/lum_sun))

#Luminosity/Mass relation

m_sun_ea = m_sun

LMstar = L/(system.star_mass*m_sun)
LMsun = lum_sun/m_sun
print("L/M difference: {}".format(LMstar/LMsun ))

#Temperature/Mass relation
Tstar = system.star_temperature
Mstar = system.star_mass*m_sun
Tsun = 1.57E7
Msun = m_sun
TMstar = Tstar/Mstar
TMsun = Tsun/Msun

print("T/M difference: {:e}".format(TMstar/TMsun))

# Største R:


K = k_B
m_H = m_p
mu = 3*(1.00784/(N_A*1000))/(4*m_H) + (4.002602/(N_A*1000))/(4*m_H)
print(mu)
T = 10 #K

totm = system.star_mass*m_sun*1.5
Rj = (G*mu*m_H*3)/(15*K*T)*totm
print("r < {:e}AU".format(Rj/AU))



L = sigma*T**4*4*pi*Rj**2
lum_sun = 3.828E26
print("L star = {}, L/Lsun = {}".format(L, L/lum_sun))

#New
mu = (1.00784/(N_A*1000))/(m_H)
star_volume= ((4/3)*pi*(system.star_radius*1000)**3)
rho_0 = system.star_mass*m_sun/star_volume
T_core = (system.star_temperature + 2*(pi/3)*G*rho_0*mu*m_H/K)*(system.star_radius*1000)**2
T_core = G*system.star_mass*m_sun*m_H/(K*system.star_radius*1000)
print("T = {:g}K".format(T_core))
#Find new luminosity and compare

Epp = 1.08*10**(-12)*rho_0*(0.745)**2*(37.2)**2
Ecno = 8.24*10**(-31)*0.745*0.002*(37.2)**2
M_core = system.star_mass*m_sun*(0.2)**3
L_rev = (Epp+Ecno)*M_core

print("L/Lsun with new method = {:g}".format(L_rev/lum_sun))
