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

col = np.load('himmelkule.npy')
img = Image.open("sample0000.png") # Open existing png

pixels = np.array(img) # png into numpy array
width, length = img.size

alpha = 70 *pi/180

x_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
x_min = -x_max

y_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
y_min = -y_max

X = np.linspace(x_min,x_max, width)
Y = np.linspace(y_min,y_max, length)
xg,yg = np.meshgrid(X,Y)

theta_0 = np.pi/2
phi_0 = pi/6
def thetha_convert(x,y):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return theta_0 - np.arcsin(np.cos(beta)*np.cos(theta_0)\
     +(y/rho)*np.sin(beta)*np.sin(theta_0))
def phi_convert(phi_0,x,y):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return phi_0 + np.arctan((x*np.sin(beta))/(rho*np.sin(theta_0)*np.cos(beta) -\
    y*np.cos(theta_0)*np.sin(beta)))
theta= thetha_convert(xg,yg)
phi = phi_convert(phi_0,xg,yg)

pic = np.zeros((length, width, 3), dtype ='uint8')

for i in range(length):
    for j in range(width):
        index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
        pic[i,j] = np.array([col[index][2],col[index][3],col[index][4]])

def new_phi(img):
    pixels = img
    width = len(pixels[0, :])
    length = len(pixels[:,1])
    alpha = 70 *pi/180
    x_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
    x_min = -x_max

    y_max = (2*np.sin(alpha/2))/(1+np.cos(alpha/2))
    y_min = -x_max

    X = np.linspace(x_min,x_max, width)
    Y = np.linspace(y_min,y_max, length)
    xg,yg = np.meshgrid(X,Y)
    phi_0 = 0
    error = 1000000
    epsilon = 4000
    sigma = pi

    while error>epsilon:
        error_array = np.array([])
        origo = phi_0
        sigma = sigma/5
        phi_0_array = [origo + sigma*i for i in [-5,-5/2,0,5/2,5]]
        for i in phi_0_array:
            theta_0 = np.pi/2
            phi_0 = i
            theta = thetha_convert(xg,yg)
            phi = phi_convert(phi_0,xg,yg)
            pic = np.zeros((length, width, 3), dtype ='uint8')
            for i in range(length):
                for j in range(width):
                    index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
                    pic[i,j] =  np.array([col[index][2],col[index][3],col[index][4]])
            error_array = np.append(error_array,np.sum((pixels-pic)**2))
        error = np.min(error_array)
        print(error)
        print(np.where(error_array == error)[0])
        phi_0 = phi_0_array[int(np.where(error_array == error)[0])]
        print(phi_0)

    return phi_0

print(new_phi(pic))
