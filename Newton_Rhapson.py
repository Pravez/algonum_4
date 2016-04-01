# -*- coding: utf-8 -*-

#Premiere partie sur Newton-Rhapson

from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

################### Newton-Raphson ######################

def Newton_Raphson(f, J, U0, N, epsilon):
    """Solves non-linear equations using Newton-Raphson method"""
    i = 0
    U = U0
    
    while (i < N and np.linalg.norm(f(U)) > epsilon):
        V = np.linalg.lstsq(J(U), -f(U))[0]
        U = U + V
        i = i + 1
    
    return U


########### Newton-Raphson with Backtracking ##########

def Newton_Raphson_Back(f, J, U0, N, epsilon):
    """Solves non-linear equations using Newton-Raphson method
        with backtracking"""
    i = 0
    U = U0
    norm_f = np.linalg.norm(f(U))
    
    while (i < N and norm_f > epsilon):
        V = np.linalg.lstsq(J(U), -f(U))[0]
        
        while np.linalg.norm(f(U + V)) > norm_f:
            V = V*(2/3)
        
        U = U + V
        i = i + 1
    
    return U

###### Test part #######

def f(x):
    return np.array([x[0]*x[0]-2, x[1]*x[1]-5])


def df(x):
    return np.array([[2*x[0], 0], [0, 2*x[1]]])

def newton_raphson_convergence(f, J, U0, N, epsilon):
    """Solves non-linear equations using Newton-Raphson method and
        returns a list containing the norm of f(U)"""
    i = 0
    U = U0
    nfU = []
    while (i < N and np.linalg.norm(f(U)) > epsilon):
        nfU.append(np.linalg.norm(f(U)))
        V = np.linalg.lstsq(J(U), -f(U))[0]
        U = U + V
        i = i + 1
    return nfU


def plot_convergence_speed(f, df, eps):
    y = newton_raphson_convergence(f, df, np.array([2, 5]), 50, eps)
    plt.plot(y)
    plt.xlabel('Iterations of Newton-Raphson method')
    plt.ylabel('$||f(U)||$')
    plt.title('Convergence speed of Newton-Raphson method without backtracking')
    plt.show()


def test_all():
    """Makes a bunch of tests only if the
        file is 'executed' directly (not imported)"""
    if __name__ == "__main__":
        print "Test on function: square root with 2 and 5"
        print "------------------------------------------"
        
        print "With numpy:\tsqrt(2) = {}\tsqrt(5) = {}" . format(np.sqrt(2), np.sqrt(5))
        res1 = Newton_Raphson(f, df, np.array([2, 5]), 50, 0.00000001)
        print "Without backtracking:\t{}" . format(res1)
        res2 = Newton_Raphson_Back(f, df, np.array([2, 5]), 50, 0.00000001)
        print "With backtracking:\t{}\n" . format(res2)
        
        print "Plotting convergence speed of Newton-Raphson method without backtracking"
        plot_convergence_speed(f, df, 0.0001)
        print "... Done"
    else:
        print "> File which contains implementation of Newton-Raphson method was imported successfully.\n\n"


test_all()
