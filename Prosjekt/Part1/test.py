import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import seaborn as sns
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
import Parameter as par

seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = mission.system
print(system.print_info())
