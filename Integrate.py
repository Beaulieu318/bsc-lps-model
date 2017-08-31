import numpy as np
import random

def array(X, Y):
    array = np.zeros((X, Y))
    for y in range(Y):
        for x in range(X):
            array[y][x] = random.random()  
    return array

def trapezium(y, increment = 1):
    summation = 0
    if len(y) != 1:
        for i in range(1, len(y)-1):
            summation += y[i]
        inte = 0.5*increment*((y[0]+y[-1]) + 2*summation)
        return inte
    else:
        return 0

def x_integrate(array):
    array_out = np.zeros((len(array), len(array[0])))
    for y in range(0, len(array)):
        for x in range(0, len(array[0])):
            array_out[y][x] = trapezium(array[y][x:])
    return array_out