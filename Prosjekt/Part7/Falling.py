import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
import atmos_part7
#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

rho = atmos_part7.density_func

#Kalkuler luftmotsand
def Fd(r,v,A):
    h = np.linalg.norm(r) #finn avstand

    """Hvis vi er lengre ut enn den interpolerte lufttettheten,
    så settes den til 0"""
    if h> system.radii[1]*1000 +15000000:
        return np.asarray([0,0,0])

    #Finn vår relative hastighet relativt til atmosfæren
    rel_vel = v -(h*2*np.pi/(system.rotational_periods[1]*day))\
     * np.asarray([-r[1]/h,r[0]/h,0])
    a = -1/2 * rho(h)*A*np.linalg.norm(rel_vel)*rel_vel

    #Hvis akselerasjonen er for stor stopper koden, skal tilsvare at fallskjermen ryker
    assert np.linalg.norm(a)<250000,'a={}, height = {}'.format(np.linalg.norm(a),\
    np.linalg.norm(r)-system.radii[1]*1000)
    return a

#Gravitasjon fra planeten
def grav_planet(r,m_ship):
    a = -G*m_ship*system.masses[1]*m_sun*(r/np.linalg.norm(r)**3)
    return a
#Lander booster, hvis raketten er over eller under grensa vi satt returnerer den
# enten 0 eller det vi velger i radiell retning fra planeten til raketten
def thrusters(r,d,thruster_boost,v):
    if np.linalg.norm(r)> system.radii[1]*1000+d:
        return np.asarray([0,0,0])
    if np.linalg.norm(r)< system.radii[1]*1000+d:
        return thruster_boost*(r/np.linalg.norm(r))

#Simuler selve landingsprosessen
def falling(initial_pos,initial_vel,dt_length,tot_crossec,mass,\
    d,nr,thruster_boost):
    #Sett parametre, konstanter og init-verdier
    m_ship = mass#kg
    A = tot_crossec
    dt = dt_length
    x = np.zeros((1,3),float)
    x[0,:] = initial_pos
    v = initial_vel
    t = 0
    initial_time=0
    a_i = (grav_planet(x[0,:],mass) + Fd(x[0,:],v,A))/m_ship
    j = 0

    #Loop over posisjon til raketten med leapfrog
    while (np.linalg.norm(x[j])-(system.radii[1]*1000)>0):
        t+=dt
        #Dersom falling brukes flere ganger gjennom landingen legges posisjons
        #arrayene sammen til en array
        x = np.append(x,[x[j,:]+v*dt + 0.5*a_i*dt**2],axis=0)

        #Hvis posisjonen til raketten er mindre enn 0, returner siste posisjon,
        #med hastighet og tid
        if ((np.linalg.norm(x[j+1])-system.radii[1]*1000)<0):
            return x[:-1],v,t

        #Nr er posisjon til fallskjermutslippen, settes til 0 etter første gang
        #den brukes
        if np.linalg.norm(x[j+1])<(system.radii[1]*1000+nr) and nr>0:
            return x[:-1],v,t
        a_iplus1 = (grav_planet(x[j+1],mass) + Fd(x[j+1],v,A) + thrusters(x[j+1],d,thruster_boost,v))\
        /m_ship
        #print(a_iplus1)
        v = v + 0.5*( a_i + a_iplus1 )*dt
        a_i = a_iplus1
        j +=1

        #print avstand og radiell hastighet til raketten
        distanse = np.linalg.norm(x[j,:])-system.radii[1]*1000
        vr = np.dot(v,-x[j,:]/np.linalg.norm(x[j,:]))

        print("x = {:.5e} moh vr = {:.5f} m/s" .format(distanse,vr))
    return x[:-1],v,t
if __name__ == '__main__':

    """ mission = mission.load('mission_after_part_6.pickle')

    # Initiate the landing sequence
    landing_sequence = mission.begin_landing_sequence()
    time0, pos_act, vel = landing_sequence.orient()""" #s,m,m/s
    x,v,t = falling(np.asarray([9472707.88605198,6375196.26580251,0]),\
    np.asarray([-2398.29145663,3563.54792079,0.]),20000,1,0.3,mission.lander_mass,20)
    fig, ax = plt.subplots()
    ax.plot(x[:,0],x[:,1],label='lander')
    r = system.radii[1]*1000 +1
    circle = plt.Circle((0,0),system.radii[1]*1000, color='peru')
    ax.add_artist(circle)

    plt.show()

    
