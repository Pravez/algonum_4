# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linspace, special
import Newton_Rhapson as NR # Newton-Raphson method

               ################# Part 4 #####################

               ########## Electrostatic equilibrium #########

N = 3 # All part

######### Question 1 ##########

def sum_diagonal(X, i):

	"Return 1 / ((X[i] - X[k])**2)"

	s = 0.0
	
	for k in range(0, N):
		if i != k:
			s += 1.0 / ((X[i] - X[k])**2)
			
	return s


def Jacobian(X):

    "Return the Jacobian (here, N = 3)"

    J = np.zeros([N, N])

    for i in range(0, N):
        for j in range(0, N):
            if i != j:
                J[i, j] = 1 / ((X[i] - X[j])**2)
            
            else: #elements of diagonal
                J[i, j] = (-1 / ((X[i] + 1)**2) + (-1 / ((X[i] - 1)**2)) - sum_diagonal(X, i))

    return J


######### Question 2 ##########


def sum_deriv(X, i):

    "Return the sum derived from E"

    s = 0.0
    
    for j in range(0, N):
        if i != j:
            s += 1.0 / (X[i] - X[j])

    return s


def sum_log(X, i):

	"Return log(|X[i] - X[j]|) with j different to i in (O, N)"

	s = 0.0
	
	for j in range(0, N):
		if i != j:
			s += math.log(math.fabs(X[i] - X[j]))

	return s
	

def Laplacian_E(X):

	"Return a vector corresponding to the Laplacian E"

	E = np.zeros(N)
	
	for i in range(0, N):
		E[i] = (1.0 / (X[i] + 1)) + (1.0 / (X[i] - 1)) + sum_deriv(X, i)
		
	return E
	

def Energy(X):

	"Return E(X1, ..., XN)"

	s = 0.0
	
	for i in range(0, N):
		s += math.log(math.fabs(X[i] + 1)) + math.log(math.fabs(X[i] - 1)) + ((1.0 / 2.0) * sum_log(X, i))
		
	return s

#def P1(x):
    
    #   "Dérivée du second polynome de Legendre"
    #return 3*x

#def P2(x):
        
    #   "Dérivée du troisième polynome de Legendre"
    #return 0.5*(15.*x*x-3)

#def P3(x):
    
    #   "Dérivée du quatrième polynome de Legendre"
    #return 0.125*(140.*x*x*x-60.*x)

#def P4(x):
        
    #"Dérivée du cinquième polynome de Legendre"
    #return 0.125*(315.*x**4-210.*x**2+15.)


def Polynomials_Legendre():

    "Plot the Legendre Polynomials"
	
    x = linspace(-1,1,100)
    
    for n in range(0, 5):
		legendre = special.legendre(n)
		y = legendre(x)
		plt.plot(x, y)
	
    
    
    #  g = plt.plot(x, P1(x), 'm', linewidth=1.0)
    #  h = plt.plot(x, P2(x), 'g', linewidth=1.0)
    #  i = plt.plot(x, P3(x), 'r', linewidth=1.0)
    #   j = plt.plot(x, P4(x), 'b', linewidth=1.0)
                        
    #   U = np.asmatrix([-0.9]).T
    #   Res = NR.Newton_Raphson_Backtrack(Energy(1), Jacobian(1), U, 10, epsilon)
    #   plt.plot(Res[0],0,'o')
                                        
    #   U = np.asmatrix([-0.9,0.9]).T
    #   Res = NR.Newton_Raphson_Backtrack(Energy(2), Jacobian(2), U, 10, epsilon)
    #   plt.plot(Res[0],0,'o')
    #   plt.plot(Res[1],0,'o')
                                                        
#       U = np.asmatrix([-0.9,0.,0.9]).T
#       Res = NR.Newton_Raphson_Backtrack(Energy(3), Jacobian(3), U, 10, epsilon)
#       plt.plot(Res[0],0,'o')
#       plt.plot(Res[1],0,'o')
#       plt.plot(Res[2],0,'o')
                                                                    
#       U = np.asmatrix([-0.8,-0.1,0.1,0.8]).T
#       Res = NR.Newton_Raphson_Backtrack(Energy(4), Jacobian(4), U, 10, epsilon)
#       plt.plot(Res[0],0,'o')
#       plt.plot(Res[1],0,'o')
#       plt.plot(Res[2],0,'o')
#       plt.plot(Res[3],0,'o')
                                                                                                    
#       plt.axis([-1,1,-5,5])


    
    
    
    
    
        plt.ylabel('Pn(X)')
        plt.xlabel("X")
        plt.title("Legendre Polynomials")
	plt.legend((g,h,i,j),('P0(X) = 0', 'P1(X) = 1', 'P2(X) = 2', 'P3(X) = 3', 'P4(X) = 4', 'P5(X) = 5'))
	plt.show()

