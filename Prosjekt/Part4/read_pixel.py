import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

img = Image.open("sample0000.png") # Open existing png

pixels = np.array(img) # png into numpy array
width = len(pixels[0, :])
length = len(pixels[:,1])


alpha = 70 *pi/180
"""
Finding the max- and min-value for x and y.
"""
x_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
x_min = -x_max

y_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
y_min = -x_max
print(x_max,y_max)

"""
Creating a grid.
"""
X = np.linspace(x_min,x_max, length)
Y = np.linspace(y_min,y_max, width)
x,y = np.meshgrid(X,Y)

theta_0 = np.pi/2
phi_0 = 0

def thetha_convert(x,y):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return theta_0 - np.arcsin(np.cos(beta)*np.cos(theta_0)\
     +(y/rho)*np.sin(beta)*np.sin(theta_0))
def phi_convert(x,y):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return phi_0 + np.arctan((x*np.sin(beta))/(rho*np.sin(theta_0)*np.cos(beta) -\
    y*np.cos(theta_0)*np.sin(beta)))

theta_array= thetha_convert(x,y)
phi_array = phi_convert(x,y)
pixels = np.arraylike(theta_array)





"""
pixles = mission.get_sky_image_pixel(theta_array,phi_array)
img2 = Image.fromarray(pixels)
img2.save('del4A.png')
"""
