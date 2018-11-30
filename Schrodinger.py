"""
Code to Solve Schrodinger's Equation for continuous and well behaved potentials.
Purpose - PHY106 Project Work
Author - Kartikeya Rambhatla

"""

from decimal import *
from math import *
import numpy as np
import Eigenvalue
from matplotlib import pyplot as plt

getcontext().prec = 100


def matrix(n):
    L = []
    for i in range(n):
        l = []
        for i in range(n):
            l.append(0)
        L.append(l)
    return L


def differential(n):
    l = matrix(n)
    for i in range(n):
        for j in range(n):

            if i == j:
                l[i][j] = -2

                if j + 1 <= n - 1:
                    l[i][j + 1] = 1

                if j - 1 >= 0:
                    l[i][j - 1] = 1
    return l


def c_matrix(c, A):
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] *= c
    return A


def add_matrix(A, B):
    C = matrix(len(A))
    for i in range(len(A)):
        for j in range(len(B)):
            C[i][j] = A[i][j] + B[i][j]
    return C


def potential_matrix(V):
    n = len(V)
    l = matrix(n)
    for i in range(n):
        for j in range(n):

            if i == j:
                l[i][j] = V[i]
    return l


def H(V, d, m, h):
    n = len(V)
    V = potential_matrix(V)

    f = -(h * h) / float(2 * m * d * d)
    D = c_matrix(f, differential(n))

    A = add_matrix(D, V)
    return A


def generate_potential(function, limit, increment, m):
    V = []
    X = []
    x = -(limit)
    while x <= limit:

        try:  # If no error is found

            y = eval(function)
            V += [y]
            X += [x]
            # print("V(",x,")  = ", y)
            x += increment


        except:
            V += [0]
            X += [x]
            x += increment
            pass

    return [V, X]


def plot(x, y, title, xlabel, ylabel):
    plt.title(str(title))
    plt.xlabel(str(xlabel))
    plt.ylabel(str(ylabel))
    plt.plot(x, y)
    plt.grid()
    plt.figtext(0.05, 0.02, "Press Ctrl+W to close the window")
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()


def numpy_solve(H_matrix, h):
    """ The following code uses built in numpy function to compute n eigenvalues.
        Note that this is a very fast method for computing the eigenvalues."""

    H_matrix = np.array(H_matrix)
    E, Psi = np.linalg.eig(H_matrix)

    idx = E.argsort()[::-1]
    E = E[idx]
    Psi = Psi[:, idx]

    E = (E.tolist())
    print("Energy values", E)
    Psi = Psi.tolist()
    return Psi


def inverse_power_solve(H_matrix, h):
    """ The following commented code uses inverse power method
        to find out only one of the energy eigenvalues that is the smallest one.
        This plots a really good graph for the wavefuntion for the ground state energy.
        Also it should be noted that it is very slow."""

    E = Eigenvalue.eigenvalues_inverse_power_method(H_matrix)

    Energy = E[0]
    Psi = E[1]

    print("Energy values", Energy)

    return Psi


def execute(H_matrix, method, h, x):
    if method == 0:

        mode = int(input(
            "Enter n if you want nth Energy-Wavefunction to be plotted (experimntal but fun) or -1 if you want the complete wavefunction plot (perfected)"))

        if mode == -1:
            Psi = numpy_solve(H_matrix, h)
            plot(x, Psi, "Psi vs x", "x", "Psi")
            print(Psi[0])
            plt.show()

        else:
            Psi = numpy_solve(H_matrix, h)
            plot(x, Psi[mode], "Psi vs x", "x", "Psi")
            print(Psi[0])
            plt.show()

    elif method == 1:
        Psi = inverse_power_solve(H_matrix, h)
        plot(x, Psi, "Psi (Ground state) vs x", "x", "Psi")
        plt.show()

    elif method == 2:
        Psi1 = numpy_solve(H_matrix, h)
        Psi2 = inverse_power_solve(H_matrix, h)
        plt.subplot(1, 2, 1)
        plot(x, Psi1, "Psi vs x [Numpy method]", "x", "Psi")
        plt.subplot(1, 2, 2)
        plot(x, Psi2, "Psi (Ground state) vs x [Inverse Power method]", "x", "Psi")
        plt.show()

    else:
        main__()


def main__():
    q = Decimal(1.60217662 * (10 ** (-19)))
    q = float(q)

    k = Decimal(9 * (10 ** (9)))
    k = float(k)

    m = Decimal(1.6737236 * (10 ** (-27)))
    m = float(m)

    h = Decimal(1.054571 * (10 ** (-34)))
    h = float(h)

    d = float(Decimal(10 ** -20))

    potential = str(input("Input your continuous potential function ((q*q*k)/x or x*x for easy visualization) : "))
    limits = float(input("Enter the highest number for the range of the graph to be plotted : "))
    increment = float(input("Enter the interval size : "))
    method = float(input("Enter 0 for numpy module solution\nEnter 1 for inverse power method solution\nEnter 2 for "
                         "both"))

    L = generate_potential(potential, limits, increment, m)
    V, x = L[0], L[1]

    H_matrix = H(V, d, m, h)

    execute(H_matrix, method, h, x)


main__()
