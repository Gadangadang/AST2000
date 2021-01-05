import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
#import part_3
#import test_part7
utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

densities = np.load("densities.npy",allow_pickle=True)
print(densities)

r = np.linspace(system.radii[1]*1000,system.radii[1]*1000+15000000,15000001)
print(len(densities),len(r))

def interpol(a,b):

    f = interpolate.interp1d(a,b,axis=0)
    return f

density_func = interpol(r,densities)
