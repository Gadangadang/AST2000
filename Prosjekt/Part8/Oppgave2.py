import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

from ast2000tools.relativity import RelativityExperiments

utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
experiments = RelativityExperiments(seed)
system = SolarSystem(seed)

planet_idx = 1 # I want to perform the experiment near planet 4
experiments.spaceship_race(planet_idx)
