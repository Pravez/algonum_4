import numpy as np
from Newton_Raphson.py import *

def elastic_force(dx,k):
    return -dx*k

def centrifugal_force(x0, y0, x, y, k):
    return (k*(x-x0), k*(y-y0))

def gravitational_force(x0 ,y0, x, y, k):
    return (-k*(x-x0)/((x-x0)**2 + (y-y0)**2)**(3/2),
            -k*(y-y0)/((x-x0)**2 + (y-y0)**2)**(3/2))

def f(U):
    return np.array([f1(U),f2(U)])

def forces_sum(U):
    f1=gravitational_force(0,0,x,y,1)
    f2=gravitational_force(1,0,x,y,0.01)
    m1=1
    m2=0.01
    f3=centrifugal_force((m1*0+m2*1)/(m1+m2),(m1*0+m2*0)/(m1+m2),1)
    return f1+f2+f3

def deriv_x(f,x0,y0,k,U):
    h=10**-10
    T=np.array([h,0])
    return ((f(U+T,x0,y0,k) - f(U-T,x0,y0,k))/(2*h))

def deriv_y(f,x0,y0,k,U):
    h=10**-10
    T=np.array([0,h])
    return ((f(U+T,x0,y0,k) - f(U-T,x0,y0,k))/(2*h))

def gravit_jac(U,x0,y0,k):
    J=np.array([deriv_x(gravitational_force,x0,y0,k,U),deriv_y(gravitational_for
ce,x0,y0,k,U)])
    return J

def centrif_jac(U,x0,y0,k):
    J=np.array([deriv_x(centrifugal_force,x0,y0,k,U),deriv_y(centrifugal_force,x
0,y0,k,U)])
    return J

def Jacobian_matrix(U):
    f1=gravit_jac(U,0,0,1)
    f2=gravit_jac(U,1,0,0.01)
    m1=1
    m2=0.01
    x0=(m1*0+m2*1)/(m1+m2)
    y0=(m1*0+m2*0)/(m1+m2)
    f3=centrif_jac(U,x0,y0,1)
    return f1+f2+f3

