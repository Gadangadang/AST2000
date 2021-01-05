import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import seaborn as sns
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem


utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

print(system.rotational_periods[0])
print(system.initial_positions)
print(day)
print(yr)
print(system.radii)
print(system.initial_positions)
print(system.initial_velocities)
print(system.rotational_periods)
