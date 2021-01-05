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

sigma = sigma


"""
(base) 1x-193-157-182-74:Part3 Sakarias$ python solarpanel.py
Area of solar panel is 0.209958m^2, with length and width 0.458212m on planet 1
"""
