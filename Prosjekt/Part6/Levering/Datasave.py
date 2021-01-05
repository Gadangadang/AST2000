import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

infile = open("sigma_noise.txt", "r")
noise_data = []
lambda_data_noise = []
for line in infile:
    word = line.split()
    noise_data.append(float(word[0]))
    lambda_data_noise.append(float(word[1]))

infile = open("spectrum_seed49_600nm_3000nm.txt", "r")
flux_data = []
lambda_data = []
for line in infile:
    word = line.split()
    flux_data.append(float(word[0]))
    lambda_data.append(float(word[1]))
np.save('flux_data',np.asarray(flux_data))
np.save('noise_data',np.asarray(noise_data))
np.save('lambda_data',np.asarray(lambda_data))
