import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem
from scipy import interpolate
from ast2000tools.shortcuts import SpaceMissionShortcuts

flux_data = np.load('lambda_data.npy')
noise_data = np.load('noise_data.npy')
lambda_data = np.load('flux_data.npy')

dt = 2400/len(lambda_data)
c /= 1E-9
k_B /= 1E-18
"""
Kode for å regne ut fluks verdiene for en lambda
"""
def F(lambdai,width,lambda_min,lambda_0):
    return 1 + (lambda_min-1)*np.exp(-(lambdai-lambda_0)**2/(2*width**2))
"""
Gjorde om temperatur til sigma
"""
def width_func(lambda_0,T,mass):
    return ((lambda_0)/c * np.sqrt(k_B*T/mass))

def calculate_fault(lambda_0_og,mass):
    """
    Definerer intervaller for parametrene jeg vil teste for.
    """
    T = np.linspace(135,475,45)
    lambda_min = np.linspace(0.6,0.99,30)
    doppler_max = lambda_0_og*1E13/(c)
    lambda_0_list = np.linspace(lambda_0_og-doppler_max,lambda_0_og+doppler_max,70)
    """
    Liste for indexene til de forskjellige parametrene.
    """
    index_list0,error_list, index_listT, index_listmin = np.array([]), np.array([]),\
     np.array([]), np.array([])
    """
    Lager kurver av alle kombinasjonene og finner kurven som gir minst feil.
    """
    for k in range(len(lambda_0_list)):
        lambda_0 = lambda_0_list[k]
        for i in range(len(T)):
            width = width_func(lambda_0_og,T[i],mass)
            min,max = int((-600+lambda_0-doppler_max-width)/dt),int((-600+lambda_0+doppler_max+width)/dt)
            for j in range(len(lambda_min)):
                error =  np.sum(((flux_data[min:max]-F(lambda_data[min:max],width,\
                lambda_min[j],lambda_0))/noise_data[min:max])**2)
                error_list = np.append(error_list,error)
                index_listT = np.append(index_listT,[i])
                index_listmin = np.append(index_listmin,[j])
                index_list0 = np.append(index_list0,[k])
    """
    Bruker at index listene til alle parametrene har like mange elemeneter som antall kurver.
    Altså, dersom kurve nr j gir best tilnærming, vil index_listT[j] gi parametrene
    brukt i å lage kurve nr j.
    """
    smallest_error = np.amin(error_list)
    index = int(np.where(error_list==smallest_error)[0][0])
    lambda_0_calc = lambda_0_list[int(index_list0[index])]
    T_calc = T[int(index_listT[index])]
    Fmin_calc = lambda_min[int(index_listmin[index])]
    v_rad_calc = c*1E-9*(lambda_0_calc-lambda_0_og)/lambda_0_og
    return Fmin_calc,lambda_0_calc,T_calc,abs(v_rad_calc)

proton_mass= 1.6726219E-27 #kg
neutron_mass=1.6749274E-27#kg

hydrogen_mass = proton_mass
carbon_mass = 6*proton_mass + 6*neutron_mass
nitrogen_mass = 7*proton_mass + 7*neutron_mass
oksygen_mass = 8*proton_mass + 8*neutron_mass

"""
Finner massene til alle molekylene.
"""
O2 = np.array([632,690,760])
H2O = np.array([720,820,940])
CO2 = np.array([1400,1600])
CH4 = np.array([1660,2200])
CO = np.array([2340])
N2O = np.array([2870])
Element = [O2,H2O,CO2,CH4,CO,N2O]
Element_names = ['O2','H2O','CO2','CH4','CO','N2O']
masses = [oksygen_mass*2,oksygen_mass+hydrogen_mass*2,carbon_mass+oksygen_mass*2,\
carbon_mass+hydrogen_mass*4,carbon_mass+oksygen_mass,nitrogen_mass*2+oksygen_mass]
i = 0
j = 1
print('Element: Fmin, lambda0, T, v_rad')
"""
Finner og plotter alle kurvene som gir best tilnærming.
"""
for elem in Element:
    mass = masses[i]
    name = Element_names[i]
    i+=1
    for lambda_0_og in elem:
        Fmin_calc,lambda_0_calc,T_calc,v_rad_calc = calculate_fault(lambda_0_og,mass)
        print('{}:    {:.2f}--{:.2f}--{:.2f}--{:.2f}'.format(str(name),Fmin_calc,lambda_0_calc,T_calc,v_rad_calc))
        print('----------------------------------')
        width = width_func(lambda_0_og,T_calc,mass)
        doppler_max = lambda_0_og*1E13/c
        min,max = int((-600+lambda_0_og-2*doppler_max-width)/dt),int((-600+lambda_0_og+2*doppler_max+width)/dt)
        plt.subplot(2,6,j)
        plt.plot(lambda_data[min:max], flux_data[min:max])
        plt.plot(lambda_data[min:max],F(lambda_data[min:max],width,Fmin_calc,lambda_0_calc), label =str(name))
        plt.legend()
        j+=1
plt.savefig('Part6b.jpeg')
plt.show()
