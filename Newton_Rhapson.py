#Premiere partie sur Newton-Rhapson

from scipy import misc
import numpy as np
import sympy

#Fonction de calcul de la jacobienne
def compute_jacobian(vectorial, x):
    h = 1.0e-4
    n = len(x)
    jacobian = np.zeros((n, n))
    f0 = vectorial(x)
    for i in range(n):
        temp = x[i]
        x[i] = temp + h
        f1 = vectorial(x)
        x[i] = temp
        jacobian[:,i] = (f1 - f0)/h
    return jacobian, f0



x1, x2, x3 = sympy.symbols('x1 x2 x3', Real=True)

compute_jacobian()