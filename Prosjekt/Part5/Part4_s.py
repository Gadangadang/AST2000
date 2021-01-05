import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts
from PIL import Image
#utils.check_for_newer_version()
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

def thetha_convert(x,y,theta_0):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return theta_0 - np.arcsin(np.cos(beta)*np.cos(theta_0)\
     +(y/rho)*np.sin(beta)*np.sin(theta_0))
def phi_convert(phi_0,x,y,theta_0):
    rho = np.sqrt(x**2 + y**2)
    beta = 2*np.arctan(rho/2)
    return phi_0 + np.arctan((x*np.sin(beta))/(rho*np.sin(theta_0)*np.cos(beta) -\
    y*np.cos(theta_0)*np.sin(beta)))

def predict_spacecraft_pointing_phi(captured_image_filename):
    """
    Here you can implement the angle prediction for challenge B
    of Part 4 of the project.
    """

    img = Image.open(captured_image_filename)
    col = np.load('himmelkule.npy')
    pixels = np.array(img)
    alphaphi = 70*pi/180
    alphatheta = alphaphi
    width, length = img.size
    x_max = (2*np.sin(alphaphi/2))/(1+np.cos(alphaphi/2))
    x_min = -x_max
    y_max = (2*np.sin(alphatheta/2))/(1+np.cos(alphatheta/2))
    y_min = -y_max
    X = np.linspace(x_min,x_max, width)
    Y = np.linspace(y_min,y_max, length)
    xg,yg = np.meshgrid(X,Y)

    error = 10000000
    epsilon = 500000
    sigma = 2*pi
    phi_0 = sigma/2
    origo = 0
    for i in range(5):
        error_array = np.array([])
        if error == 1000000:
            sigma = sigma/60
            phi_0_array = [origo + sigma*i for i in range(60)]
        else:
            sigma = sigma/4
            origo = phi_0 - 2*sigma
            phi_0_array = [origo + sigma*i for i in [0,1,2,3,4]]
        for i in phi_0_array:
            theta_0 = pi/2
            theta = thetha_convert(xg,yg,theta_0)
            phi = phi_convert(i,xg,yg,theta_0)
            pic = np.zeros((length, width, 3), dtype ='uint8')
            for i in range(length):
                for j in range(width):
                    index = mission.get_sky_image_pixel(theta[i,j],phi[i,j])
                    pic[i,j] =  np.array([col[index][2],col[index][3],col[index][4]])
            error_array = np.append(error_array,np.sum((pixels-pic)**2))
        error = np.min(error_array)
        if len(np.where(error_array == error)[0])>1:
            phi_0 = phi_0_array[int(np.where(error_array == error)[0][0])]
        else:
            phi_0 = phi_0_array[int(np.where(error_array == error)[0])]
    return phi_0


def predict_spacecraft_velocity(mission,
                                measured_star_1_doppler_shift,
                                measured_star_2_doppler_shift):

    m_lamda1 = measured_star_1_doppler_shift
    m_lamda2 = measured_star_2_doppler_shift
    lamda1, lamda2 = mission.star_doppler_shifts_at_sun[0],mission.star_doppler_shifts_at_sun[1]
    theta1,theta2 = mission.star_direction_angles[0],mission.star_direction_angles[1]
    def radi(theta):
        return theta*pi/180

    theta1,theta2 = radi(theta1),radi(theta2)
    def dp_shift_3(lamda,m_lamda):
        r = c_AU_pr_s*yr*((lamda-m_lamda))/mission.reference_wavelength
        return r

    def unit_convert(phi1,phi2,vel):
        A = np.array([[np.sin(phi2),-np.sin(phi1)],[-np.cos(phi2),np.cos(phi1)]])
        vel_vec = np.array([vel[0],vel[1]])
        new_vel = 1/np.sin(phi2-phi1)*np.dot(A,vel_vec)
        return new_vel

    r1 = dp_shift_3(lamda1,m_lamda1)
    r2 = dp_shift_3(lamda2,m_lamda2)
    d_vel = [r1,r2]
    d_vel = tuple(d_vel)
    vx,vy = unit_convert(theta1,theta2,d_vel)[0],unit_convert(theta1,theta2,d_vel)[1]

    return np.asarray([vx,vy])


def predict_spacecraft_position(time_of_measurement, measured_distances):
    """
    Here you can implement the position prediction for challenge D
    of Part 4 of the project.
    """

    # Load the exact planet trajectories from Part 2
    times, planet_positions = np.load('Planet_orbits.npy',allow_pickle=True)

    #Interpolate a coordinate at a given time
    def interpol(a,b,t):
        t = t
        f = interpolate.interp1d(a,b,axis=0)
        return f(t)

    #Find position of each object, with given radii of circle for the reference points
    def find(t,tidspunkt):
        p0x = planet_positions[:][0][0]
        p0y = planet_positions[:][1][0]
        p1x = planet_positions[:][0][4]
        p1y = planet_positions[:][1][4]
        s0x = np.zeros(len(p0x))
        s0y = np.zeros(len(p0x))
        x1 = interpol(t,p0x,tidspunkt)
        y1 = interpol(t,p0y,tidspunkt)
        x2 = interpol(t,p1x,tidspunkt)
        y2 = interpol(t,p1y,tidspunkt)
        x3 = interpol(t,s0x,tidspunkt)
        y3 = interpol(t,s0y,tidspunkt)
        r1 = measured_distances[0]
        r2 = measured_distances[4]
        r3 = measured_distances[-1]
        return x1,y1,x2,y2,x3,y3,r1,r2,r3

    #Trilaterate position
    def position(t,tidspunkt):
        x1,y1,x2,y2,x3,y3,r1,r2,r3 = find(t,tidspunkt)
        A = 2*x2 - 2*x1
        B = 2*y2 - 2*y1
        C = r1**2 - r2**2 - x1**2 + x2**2 - \
        y1**2 + y2**2

        D = 2*x3 - 2*x2
        E = 2*y3 - 2*y2
        F = r2**2 - r3**2 - x2**2 + x3**2 - \
        y2**2 + y3**2

        x = (C*E - F*B) / (E*A - B*D)
        y = (C*D - A*F) / (B*D - A*E)

        return x,y
    x,y = position(times,time_of_measurement)

    return np.asarray([x,y])


def launch_and_orient_spacecraft(mission,
                                 launch_direction,
                                 time_of_launch):
    """
    This function performs a rocket launch before predicting and
    verifying the orientation of the spacecraft.
    """

    codes = [8208] # The codes recieved from the group teacher (in this case just one)
    shortcuts = SpaceMissionShortcuts(mission, codes)
    shortcuts.place_spacecraft_on_escape_trajectory(150000,15,9E-10,20000000,pi/3, 400000)

    # Predict pointing angle after launch
    captured_image_filename = 'sky_picture.png'
    mission.take_picture(filename=captured_image_filename)

    predicted_phi = predict_spacecraft_pointing_phi(captured_image_filename)

    # Predict velocity after launch
    measured_star_1_doppler_shift, \
    measured_star_2_doppler_shift  \
      = mission.measure_star_doppler_shifts()

    predicted_velocity = predict_spacecraft_velocity(mission,
                                                     measured_star_1_doppler_shift,
                                                     measured_star_2_doppler_shift)

    # Predict position after launch
    time_of_measurement = mission.time_after_launch
    measured_distances = mission.measure_distances()
    predicted_position = predict_spacecraft_position(time_of_measurement,
                                                     measured_distances)

    # Verify manually inferred orientation values
    mission.verify_manual_orientation(predicted_position,
                                      predicted_velocity,
                                      predicted_phi)


# Prevent the following code from executing when calling `import part_4`
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available


    launch_direction = 0
    time_of_launch = 15.5


    # Perform launch and orientation
    launch_and_orient_spacecraft(mission,
                                 launch_direction,
                                 time_of_launch)
