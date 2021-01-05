import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)


print("Planet mass = {:.4e} kg".format(system.masses[1]*m_sun))
print("Planet radii = {:.3e} km".format(system.radii[1]))
print("Planet volume = {:.3e} km^3".format(4/3*np.pi*system.radii[1]**3))
print("Lander mass = {}kg".format(mission.lander_mass))
