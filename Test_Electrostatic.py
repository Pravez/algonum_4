# -*- coding: utf-8 -*-
from Electrostatic import *
import unittest
import matplotlib.pyplot as plt
import tkMessageBox
import Tkinter as tk

# Suppression de la fenetre en trop de tkinter
root = tk.Tk()
root.withdraw()

X_test = np.array([0.2, 0.4, 0.5])

class Tests_Electrostatic(unittest.TestCase):

    def test_energy_borders(self):
        result = NR.Newton_Raphson_Back(Laplacian_E, Jacobian, X_test, 50, 0.0001)

        a = np.array([0.8, 0.6, 0.4])
        b = np.array([0.5, 0.6, -0.7])
        c = np.array([-0.1, 0.8, 0.6])

        assert(Energy(result) < 0)
        assert(Energy(a) < 0)
        assert(Energy(b) < 0)
        assert(Energy(c) < 0)

        jacobian = Jacobian(result)
        eigvals = np.linalg.eigvals(jacobian)
        for i in range(0, eigvals.size):
            assert(eigvals[i] < 0)

        if(tkMessageBox.askyesno("Polynômes de legendre et évolution de l'énergie","Faire un plot des polynômes de legendre et de l'évolution de l'énergie ?")):
            root.quit()
            root.destroy()
            Polynomials_Legendre()

            x = np.linspace(-0.9, 0.9, N)
            y = result

            plt.plot(x, y)
            plt.ylabel('Y (Values of resolution system)')
            plt.xlabel("X")
            plt.title("Evolution of Energy")
            plt.show()




if __name__ == '__main__':
    unittest.main()
