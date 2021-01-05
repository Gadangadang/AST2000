import numpy as np
import matplotlib.pyplot as plt

N = 10**5 #Antall partikler
T = 3000 #Temp i Kelvin
a = 0 #m/s
b = 3*10**4 #m/s
n = 100000 #Antall iterasjoner
k = 1.38*10**(-23) #Boltzmanns konstant
m = 1.67*10**(-27) #Masse H2 gass

sigma = np.sqrt((k*T)/m)

def P(x):
    B = (1/((2*np.pi*sigma**2)**(3/2)))*np.exp(-x**2/(2*sigma**2))
    return B

g = np.linspace(a,b,n)

def trapezoidal(f, a, b, n):
    h = float(b-a)/n
    result = 0.5*f(a) + 0.5*f(b)
    for i in range(1, n):
        result += f(a + i*h)
    result *= h
    return result

print(trapezoidal(P,a,b,n)*N)

plt.plot(g,P(g))
plt.xlabel("Velocity m/s")
plt.ylabel("Probability distribution")
plt.grid()
plt.show()
