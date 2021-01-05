import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
print(seed)

rocket_info = np.load("rocket_info.npy",allow_pickle=True)



semimajor_axis = np.max(abs(rocket_info[1]))
semiminor_axis = np.max(abs(rocket_info[2]))

print(semiminor_axis/semimajor_axis)
if semimajor_axis > semiminor_axis:
    e = np.sqrt(1-(semimajor_axis/semiminor_axis)**2)
else:
    e = np.sqrt(1-(semiminor_axis/semimajor_axis)**2)
print("a = {}, b = {}, e = {}".format(semimajor_axis,semiminor_axis,e))


vel_der = np.gradient(rocket_info[3])
vel_der_change = np.where(np.diff(np.sign(vel_der)))[0]
print("P/yr = {}".format(len(vel_der_change)/2/7))

print(rocket_info[2])

plt.plot(rocket_info[1],rocket_info[2],"b",label="Rocket pos")
plt.title('Rocket position')
plt.xlabel('yr')
plt.ylabel('AU/yr')
plt.legend()
plt.savefig('Position of rocket.jpeg')
plt.show()


plt.plot(rocket_info[0],rocket_info[3],"b",label="Radial velocity")
plt.title('Rocket radial velocity')
plt.xlabel('yr')
plt.ylabel('AU/yr')
plt.legend()
plt.savefig('Radial velocity of rocket.jpeg')
plt.show()

plt.plot(rocket_info[0],vel_der,"b",label="Radial acceleration")
plt.title('Rocket radial velocity')
plt.xlabel('yr')
plt.ylabel('AU/yr^2')
plt.legend()
plt.savefig('Radial acceleration of rocket.jpeg')
plt.show()
