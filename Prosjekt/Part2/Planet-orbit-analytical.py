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

# Extract associated SolarSystem object

e = [i for i in system.eccentricities]
a = [i for i in system.semi_major_axes]
print(len(a),len(e))
name = [i for i in range(8)]
x= zip(e,a,name)

def orbit_analytical(e,a,f):
    r = (a*(1-e**2))/(1+e*np.cos(f))
    return r

angle = np.linspace(0,2*np.pi,1000)

for i in x:
    r = orbit_analytical(i[0],i[1],angle)
    x = r*np.cos(angle)
    y = r*np.sin(angle)
    plt.plot(x, y, label="{}".format(i[2]))

#plt.axis('equal')
plt.grid()
plt.legend()
plt.xlabel("X-position (AU)")
plt.ylabel("Y-position (AU)")
plt.savefig("Analytical-orgbit.jpeg")
plt.show()
