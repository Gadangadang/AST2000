import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
"""
Her bruker vi bare geometriske sammenhenger og argumentasjonen fra bloggen.
"""

def coor_change(init_pos,time_elapsed):
    omega = 2*np.pi/system.rotational_periods[1]*day
    v_phi = 2*np.pi/(system.rotational_periods[1]*day)*system.radii[1]*1000
    x,y,z = init_pos
    r,phi, theta = np.array([np.linalg.norm(init_pos),np.arctan(y/x),np.arccos(z/\
    np.linalg.norm(init_pos))])
    new_phi = phi + v_phi*time_elapsed
    new_pos = np.array([r*np.sin(theta)*np.cos(new_phi),r*np.sin(theta)*np.sin(new_phi),
    r*np.cos(theta)])
    return new_pos
