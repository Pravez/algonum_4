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
	

def Polynomials_Legendre():

	"Plot the Legendre Polynomials"
	
	x = linspace(-1,1,100)

	for n in range(0, 5):
		legendre = special.legendre(n)
		y = legendre(x)
		plt.plot(x, y)
	
	plt.legend(('P0(X) = 0', 'P1(X) = 1', 'P2(X) = 2', 'P3(X) = 3', 'P4(X) = 4', 'P5(X) = 5'))
	plt.show()
	

	
######## Test part ###########


X_test = np.array([0.2, 0.4, 0.5]) # For all test part

print "Test of Jacobian :"
print ""
print Jacobian(X_test)
print ""
print ""

print "Test of Energy's Laplacian :"
print ""
print Laplacian_E(X_test)
print ""
print ""

print "Resolution of the system with Newton-Raphson method (backtracking) :"
print ""
res = NR.Newton_Raphson_Back(Laplacian_E, Jacobian, X_test, 50, 0.0001)
print res
print ""
print ""

print "Plotting the Legendre polynomials :"
print ""
Polynomials_Legendre()
	


############ Question 3 ############


a = np.array([0.8, 0.6, 0.4])
b = np.array([0.5, 0.6, -0.7])
c = np.array([-0.1, 0.8, 0.6])

print "Energy of {}:" . format(res)
print ""
print Energy(res)
print ""
print ""

print "Energy of {}:" . format(a)
print ""
print Energy(a)
print ""
print ""

print "Energy of {}:" . format(b)
print ""
print Energy(b)
print ""
print ""

print "Energy of {}:" . format(c)
print ""
print Energy(c)
print ""
print ""

print "Jacobian of {}:" . format(res)
print ""
print Jacobian(res)
print ""
print ""

print "Eigenvalues ​​of the Jacobian:"
print ""
print np.linalg.eigvals(Jacobian(res))
print ""
print ""

print "We obtain negative values ​​so you get maximum energy."
print ""				
	

# Evolution of Energy

x = np.linspace(-0.9, 0.9, N)
y = res

print x
print y
print ""

plt.plot(x, y)
plt.show()	
	
	
	
	
	
	
	
	
	
	
	
	
	
			
		
		


