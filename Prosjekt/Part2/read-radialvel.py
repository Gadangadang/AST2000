import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

#utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)

L = 10
sigma = 0.00001/5
sigma_P = 7.5/5
sigma_t0 = 5.35/5
max_val_v = 0.000011+sigma
min_val_v = 0.000011-sigma
max_val_P = 7.5+sigma_P
min_val_P = 7.5-sigma_P
max_val_t0 = 5.35+sigma_t0
min_val_t0 = 5.35-sigma_t0

val_of_cm = 46


infile = open("Radial_vel.txt", "r")
infile.readline()

time = []
vel = []
for line in infile:
    word = line.split()
    time.append(float(word[0]))
    vel.append(float(word[1]))


time = np.asarray(time)
vel_imp = np.asarray(vel)



for i in range(len(vel_imp)):
    vel_imp[i] -= val_of_cm
N = len(time)



def minste_kvadrat(P,v_star, t_0, vel_imp):
    value = 0
    for i in range(len(time)):
        a = (vel_imp[i]-v_star*np.cos((2*np.pi/P)*(time[i]-t_0)))**2
        value += a
    return P,v_star, t_0,value

"""
d_value_P = (max_val_P-min_val_P)/(L-1)
d_value_v = (max_val_v-min_val_v)/(L-1)
d_value_t0 = (max_val_v-min_val_v)/(L-1)

P_values = np.asarray([min_val_P + d_value_P*i for i in range(L)])
v_star_values = np.asarray([min_val_v + d_value_v*i for i in range(L)])
t_0_values = np.asarray([min_val_t0 + d_value_t0*i for i in range(L)])

minste_kvadrat_liste = []

for n in range(len(P_values+1)):
    for g in range(len(P_values+1)):
        for j in range(len(P_values+1)):
            b = minste_kvadrat(P_values[n],v_star_values[g],t_0_values[j],vel_imp)
            minste_kvadrat_liste.append(b)

minste_kvadrat_liste = np.asarray(minste_kvadrat_liste)

for i in range(len(minste_kvadrat_liste)):
    if minste_kvadrat_liste[i,3] == min(minste_kvadrat_liste[:,3]):
        index = np.where(minste_kvadrat_liste[i,3])



print(minste_kvadrat_liste[index[0]])
P,v_star,t_0,value =minste_kvadrat_liste[index[0]][0]
print("Period = {},star velocity = {}, t_0 = {}".format(P,v_star,t_0,value))
def V(t):
    return v_star*np.cos((2*np.pi*(t-t_0)/P))
"""
plt.plot(time,vel_imp,label="Data")
#plt.plot(time,V(time),label="Estimated curve")
plt.plot()

plt.xlabel("Years")
plt.ylabel("AU/Years")
plt.legend()
#plt.savefig("read-vals.jpeg")
plt.show()
