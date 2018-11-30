"""
Code to Solve Eigenvalue Problems using Power and Inverse Power methods.
Purpose - PHY106 Project Work
Author - Kartikeya Rambhatla

"""

import random
import numpy as np


# A = [[100, 0, 0], [0, -1, 0], [0, 0, 600]]


def multiply_vector(A, v):
    z = []
    for i in range(len(A)):
        x = 0.0
        for j in range(len(A)):
            x += A[i][j] * v[j]
        z += [x]
    return z


def norm(v):
    n = 0.0
    for i in range(len(v)):
        n += v[i] * v[i]

    return n ** 0.5


def normalize(v):
    return scale(v, 1 / norm(v))


def scale(v, k):
    x = []
    for i in range(len(v)):
        x += [k * v[i]]
    return x


def random_normalized_vector(dimension, o=2):
    """ Generate a random vector v0 """
    v = []
    for i in range(dimension):
        v += [random.uniform(-o, o)]

    ''' Normalize it'''
    v = normalize(v)
    return v


def mod(v):
    x = []
    for i in range(len(v)):
        x += [abs(v[i])]
    return x


def eigenvalues_power_method(A, o=2):
    v = random_normalized_vector(len(A), o)
    v1 = []
    l = []
    l += [v1]
    l += [v]
    print(l)
    count = 0
    while mod(l[0]) != mod((l[1])):
        print(l[1], " end")
        z = multiply_vector(A, l[1])  # Operate matrix on v
        print(z)
        print(norm(z))
        x = normalize(z)
        print(z, "2")
        l[0] = l[1]
        l[1] = x
        count += 1

        print("Number of iterations : ", count)

    if count % 2 != 0:
        print("Eigenvalue : ", -1 * norm(z))

    else:
        print("Eigenvalue : ", norm(z))

    return l[1]


def eigenvalues_inverse_power_method(A, o=1):
    v = random_normalized_vector(len(A), o)
    v1 = []
    l = []
    l += [v1]
    l += [v]
    count = 0
    while mod(l[0]) != mod((l[1])):

        A = np.array(A)
        b = np.array(l[1])
        z = np.linalg.solve(A, b)
        if np.allclose(np.dot(A, z), b):
            z = z.tolist()
        else:
            print("Cannot solve the problem")
            break
        x = normalize(z)
        l[0] = l[1]
        l[1] = x
        count += 1

    # print("Number of iterations : ", count)

    eigenvector = l[1]

    ''' Find the sign of the eigenvalue'''
    Av = np.dot(np.array(A), np.array(eigenvector))
    lambdav = np.dot(1.0 / norm(z), np.array(eigenvector))

    for i in range(len(Av)):
        if Av[i] != 0:
            sign = np.sign(Av[i] * lambdav[i])
            break

    eigenvalue = sign * 1.0 / norm(z)

    # print("Eigenvalue : ", eigenvalue)

    return [eigenvalue, eigenvector]

# print(eigenvalues_inverse_power_method(A))
