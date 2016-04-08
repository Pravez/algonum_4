import numpy as np
import numdifftools as nd
from Newton_Raphson.py import *


def elastic_force(dx,k):
    return -dx*k

def centrifugal_force(x0, y0, x, y, k):
    return (k*(x-x0), k*(y-y0))

def gravitational_force(x0 ,y0, x, y, k):
    return (-k*(x-x0)/((x-x0)**2 + (y-y0)**2)**(3/2),
            -k*(y-y0)/((x-x0)**2 + (y-y0)**2)**(3/2))

def f(x, y):
    b=(0.01/1.01, 0)

    fx= gravitational_force(0, 0, x, y, 1)[0] +
    gravitational_force(1, 0, x, y, 0.01)[0] +
    centrifugal_force(b[0], b[1], x, y, 1)[0]

    fy= gravitational_force(0, 0, x, y, 1)[1] +
    gravitational_force(1, 0, x, y, 0.01)[1] +
    centrifugal_force(b[0], b[1], x, y, 1)[1]

    return (fx,fy)

def J(


def equilibrium_points(Ux0 ,Uy0):
