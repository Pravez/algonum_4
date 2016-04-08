# -*- coding: utf-8 -*-

#Premiere partie sur Newton-Rhapson
import numpy as np

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

