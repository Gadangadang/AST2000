from scipy.stats import norm
from math import sqrt, pi
import numpy as np

#Funksjon for gaussiske fordelingen av noe
def f(x):
    p = (1/(sqrt(2*pi)*sigma))*np.exp((-1/2)*((x-mu)/sigma)**2)
    return p

#Numerisk integrasjonsløkke
def trapezoidal(f, a, b, n):
    h = float(b-a)/n
    result = 0.5*f(a) + 0.5*f(b)
    for i in range(1, n):
        result += f(a + i*h)
    result *= h
    return result

#Variabler for løkka
sigma = 2
mu =10
n = 1000
a = mu-sigma
a2 = mu-2*sigma
a3 = mu-3*sigma
b = mu+sigma
b2 = mu+2*sigma
b3 = mu+3*sigma


P1 = trapezoidal(f,a,b,n)*100
P2 = trapezoidal(f,a2,b2,n)*100
P3 = trapezoidal(f,a3,b3,n)*100
print("Difference of standard deviation: {:.2f}% {:.2f}% {:.2f}%"\
      .format(P1,P2,P3))
