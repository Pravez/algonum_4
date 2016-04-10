# -*- coding: utf-8 -*-
from Lagrange import *
import unittest

display_forces()

class Tests_Electrostatic(unittest.TestCase):

    def test__forces_sum(self):
        U = np.array([1.5, 0])
        print ("\n Example from the project: \n")
        print ("The total force applied on the two solids is : ")
        print (forces_sum(U))
        print ("Its Jacobian matrix is: ")
        print (jacobian_matrix(U))
        print ("By applying the Newton Raphson method, the result is: ")
        print (Newton_Raphson(forces_sum, jacobian_matrix, U, 50, 10 ** (-11))[0])

    def test__gravit_jac(self):
        print("Tests about the generation of the gravitational jacobian matrix")
        for i in range(0, 100):
            gen = gravit_jac(np.array([np.random.randint(1), np.random.randint(1)]), np.random.randint(1),
                             np.random.randint(1), 1)
            assert (gen[0][1] == gen[1][0])

    def test__centrif_jac(self):
        print("Tests about the generation of the centrifugal jacobian matrix")
        for i in range(0, 100):
            gen = centrif_jac(np.array([np.random.randint(1), np.random.randint(1)]), np.random.randint(1),
                              np.random.randint(1), 1)
            assert (gen[0][1] == gen[1][0])


if __name__ == '__main__':
    unittest.main()