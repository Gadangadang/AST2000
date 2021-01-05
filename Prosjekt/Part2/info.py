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

print(system.print_info())
"""

for i in range(len(system.masses)):
    affect = system.masses[i]/system.radii[i]**2
    print("The affect of planet {} is {}".format(i,affect))
"""
