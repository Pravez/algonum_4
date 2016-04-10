import numpy as np

def Newton_Raphson(f, J, U0, N, epsilon):
    """Solves non-linear equations using Newton-Raphson method"""
    i = 0
    U = U0
    
    while (i < N and np.linalg.norm(f(U)) > epsilon):
        V = np.linalg.lstsq(J(U), -f(U))[0]
        U = U + V
        i = i + 1
        
    return U

#----------Computation of the Lagrangian points----------

##
#The function elastic_force returns a vector containing -k*U.
##
def elastic_force(U,k):
    return -k*U

##
#The function centrifugal_force returns a vector containing the centrifugal force application on U.
##
def centrifugal_force(U,x0,y0,k):
    return np.array([k*(U[0] - x0), k*(U[1] - y0)])

##
#The function gravitational_force returns a vector containing the gravitational force application on U.
##
def gravitational_force(U,x0,y0,k):
    n = np.sqrt(((U[0] - x0)**2 + ((U[1] - y0)**2))**3)
    return np.array([-k*(U[0] - x0)/n, ((-k)*(U[1] - y0))/n])

##
#The function f returns the value of a function f valued at U.
##
def f(U):
    return np.array([f1(U),f2(U)])

##
#The function forces_sum returns the sum of the forces.
##
def forces_sum(U):
    f1=gravitational_force(U,0,0,1)
    f2=gravitational_force(U,1,0,0.01)
    m1=1
    m2=0.01
    f3=centrifugal_force(U,(m1*0+m2*1)/(m1+m2),(m1*0+m2*0)/(m1+m2),1)
    return f1+f2+f3

#----------Other important functions----------

##
#The function deriv_x returns the derivative of f(U) on x with parameters x0, y0 and k.
##
def deriv_x(f,x0,y0,k,U):
    h=10**-10
    T=np.array([h,0])
    return ((f(U+T,x0,y0,k) - f(U-T,x0,y0,k))/(2*h))

##
#The function deriv_y returns the derivative of f(U) on y with parameters x0, y0 and k.
##
def deriv_y(f,x0,y0,k,U):
    h=10**-10
    T=np.array([0,h])
    return ((f(U+T,x0,y0,k) - f(U-T,x0,y0,k))/(2*h))

##
#The function gravit_jac returns the jacobian matrix of the gravitational force, using parameters x0, y0 and k.
##
def gravit_jac(U,x0,y0,k):
    J=np.array([deriv_x(gravitational_force,x0,y0,k,U),deriv_y(gravitational_force,x0,y0,k,U)])
    return J

##
#The function centrif_jac returns the jacobian matrix of the centrifugal force, using parameters x0, y0 and k.
##
def centrif_jac(U,x0,y0,k):
    J=np.array([deriv_x(centrifugal_force,x0,y0,k,U),deriv_y(centrifugal_force,x0,y0,k,U)])
    return J

##
#The function jacobian_matrix returns the jacobian matrix of the function forces_sum.
##
def jacobian_matrix(U):
    f1=gravit_jac(U,0,0,1)
    f2=gravit_jac(U,1,0,0.01)
    m1=1
    m2=0.01
    x0=(m1*0+m2*1)/(m1+m2)
    y0=(m1*0+m2*0)/(m1+m2)
    f3=centrif_jac(U,x0,y0,1)
    return f1+f2+f3

#----------Tests----------

def forces_display():
    U=np.array([1,1])
    print ("The elastic force is: ")
    print (elastic_force(U, 1))
    print ("The centrifugal force is: ")
    print (centrifugal_force (U,0,0,1))
    print ("The gravitational force is: ")
    print (gravitational_force (U,0,0,1))
    return

def test__forces_sum():
    U=np.array([1.5,0])
    print ("\n Example from the project: \n")
    print ("The total force applied on the two solids is : ")
    print (forces_sum(U))
    print ("Its Jacobian matrix is: ")
    print (jacobian_matrix(U))
    print ("By applying the Newton Raphson method, the result is: ")
    print (Newton_Raphson(forces_sum,jacobian_matrix,U,50,10**(-11))[0])
    return

def test__gravit_jac():
    print("\n Tests about the generation of the gravitational jacobian matrix: \n")
    cpt=0
    for i in range(0,100):
        gen=gravit_jac(np.array([np.random.randint(1),np.random.randint(1)]),np.random.randint(1),np.random.randint(1),1)
        assert(gen[0][1] == gen[1][0])
        cpt=cpt+1
        print("For the gravitational force :\n%s successful test(s) out of 100.\n" % cpt)
        return

def test__centrif_jac():
    print("\n Tests about the generation of the centrifugal jacobian matrix: \n")
    cpt=0
    for i in range(0,100):
        gen=centrif_jac(np.array([np.random.randint(1),np.random.randint(1)]),np.random.randint(1),np.random.randint(1),1)
        assert(gen[0][1] == gen[1][0])
        cpt=cpt+1
        print("For the centrifugal force :\n%s successful test(s) out of 100.\n" % cpt)
        print("End of tests about the computation of Lagrangian points.\n")
        return

def test__Lagrange():
    forces_display()
    test__forces_sum()
    test__gravit_jac()
    test__centrif_jac()
    return

test__Lagrange()
