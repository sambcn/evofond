#!/usr/bin/env python3

import numpy as np
import scipy.optimize as op

def f(listOfArgs, ):
    """
    Conditions de fonctionnements :
    - même nombre de sortie que d'entrée
    - chaque sortie est de
    """
    a, b = np.array(listOfArgs[0]), np.array(listOfArgs[1])
    print("test")
    print(listOfArgs)
    x = []
    x.append(a*(b+np.ones(2))-4*np.ones(2))
    x.append((a-b)**2-np.ones(2))
    return x

def func(x):

    return [x[0] * np.cos(x[1]) - 4, x[1] * x[0] - x[1] - 5]


print(f([2, 1]))

root = op.fsolve(f, [np.array([1., 1., 1.]), np.array([2., 2., 2.])])
print(root)
print(f(root))
print(np.isclose(f(root), [0.0, 0.0]))