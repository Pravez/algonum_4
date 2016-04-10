# -*- coding: utf-8 -*-
from Newton_Rhapson import *
import unittest
import matplotlib.pyplot as plt
import tkMessageBox
import Tkinter as tk

#Suppression de la fenetre en trop de tkinter
root = tk.Tk()
root.withdraw()

epsilon = 10**-10
#Definition de la fonction et de la jacobienne pour ce cas
def f(x):
    return np.array([x[0]*x[0]-2, x[1]*x[1]-5])

def df(x):
    return np.array([[2*x[0], 0], [0, 2*x[1]]])

def speed_convergence():
    # recuperation de la liste de valeurs pour la convergence
    y = newton_raphson_convergence(f, df, np.array([2, 5]), 50, epsilon)
    # mise en place du plot
    plt.plot(y)
    plt.xlabel('Iterations of Newton-Raphson method')
    plt.ylabel('$||f(U)||$')
    plt.title('Convergence speed of Newton-Raphson method without backtracking')
    plt.show()


class Test_Newton_Raphson(unittest.TestCase):

    print('Tests sur Newton-Raphson')

    def test_NR_without_backtracking(self):
        print('Tests sans backtracking...')
        NR_result = Newton_Raphson(f, df, np.array([2, 5]), 50, epsilon)

        #Testing results with np.sqrt
        assert(np.abs(NR_result[0]-np.sqrt(2)) < epsilon)
        assert(np.abs(NR_result[1]-np.sqrt(5)) < epsilon)
        print('Fait')

    def test_NR_with_backtracking(self):
        print('Tests avec backtracking...')
        NR_result_back = Newton_Raphson_Back(f, df, np.array([2, 5]), 50, epsilon)

        #Testing results with np.sqrt once again
        assert(np.abs(NR_result_back[0] - np.sqrt(2)) < epsilon)
        assert(np.abs(NR_result_back[1] - np.sqrt(5)) < epsilon)
        print('Fait')

    def test_speed_convergence(self):
        if(tkMessageBox.askyesno("Vitesse de convergence","Tester la vitesse de convergence ? (une figure sera produite)")):
            root.quit()
            root.destroy()
            speed_convergence()


if __name__ == '__main__':
    unittest.main()