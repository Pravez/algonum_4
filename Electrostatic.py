# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linspace, special
import Newton_Rhapson as NR # Newton-Raphson method

               ################# Part 4 #####################

               ########## Electrostatic equilibrium #########

N = 3

######### Question 1 ##########


def Jacobian(X):

    J = np.zeros([N, N])

    for i in range(0, N):
        for j in range(0, N):
            if i != j:
                J[i, j] = 1 / ((X[i] - X[j])**2)
            
            else: #elements of diagonal
                J[i, j] = (-1 / ((X[i] + 1)**2) + (-1 / ((X[i] - 1)**2)))

    return J


######### Question 2 ##########


def sum_deriv(X, i):

    "Compute the sum derived from E"

    s = 0.0
    
    for j in range(0, N):
        if i != j:
            s += 1.0 / (X[i] - X[j])

    return s
