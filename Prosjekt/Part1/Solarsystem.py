import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem

seed = utils.get_seed('Sgfrette')
system = SolarSystem(seed)


print(system.print_info())
