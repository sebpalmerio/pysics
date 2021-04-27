"""
   This module is the hub for the library
   For details and applications of these functions, check documentation on GitHub
"""

from numpy import arange, empty, array
from math import floor

def trapezoid(f, a, b, N):
    """
    Implements the trapezoid integration scheme. 
    Arguments:
        f: function that returns an integer value
        a: left integration bound
        b: right integration bound
        N: number of partitions to divide the interval into (i.e. how many slices are we seperating our interval into)

    trapezoid: function, float, float, int -> float

    """
    h = (b-a)/N
    sum = (f(a)+f(b))/2
    for i in range(1,N):
        sum += f(a+i*h)
    
    return sum

def simpsons(f, a, b, N):
    """
    Implements Simpson's  integration scheme. 
    Arguments:
        f: function that returns an integer value
        a: left integration bound
        b: right integration bound
        N: number of partitions to divide the interval into (i.e. how many slices are we seperating our interval into)

    simpsons: function, float, float, int -> float

    """    
    h = (b-a)/N
    sum = (f(a)+f(b))*(h/3)
    for i in range(1,N):
        sum += 4 * f(a+i*h)
        sum += 2 * f(a+i*h)

    return sum

def RK4(f, ti, tf, h, yinit):
    '''
    Implements the fourth order Runge Kutta method for solving an initial value
      problem for a system of n ordinary differential equations. f is the vector
      function, ti and tf are the initial and final times, respectively, h is the
      step size, and yinit are the initial conditions.
    
    RK4: Func Float Float Float List(Float(s)) -> Array
    requires: ti < tf

    Example
    -------
    dy/dt = 5*y(t), y(0) = 0
    def f(t, y):
        return array([5*y[0]])

    RK4(f, 0.0, 1.0, 0.5, [0.0]) returns an array:
     index 0 is time values
     index [1][0] are values of y(t)
     index [1][1] are values of y'(t)
    '''  
    N_steps = int((tf-ti)/h) # number of steps
    t = h * arange(N_steps) # grid for the independent variable   
    y = empty([len(yinit), len(t)]) #(n, N_step) array for the dependent variable
    ys = array(yinit) # n-dim buffer array, set to the initial conditions

    for i, ts in enumerate(t): # RK4 integrator
        y[:,i] = ys
        k1 = h*f(ts, ys)
        k2 = h*f(ts + 0.5*h, ys + 0.5*k1)
        k3 = h*f(ts + 0.5*h, ys + 0.5*k2)
        k4 = h*f(ts + h, ys + k3)
        ys += (k1 + 2*k2 + 2*k3 + k4)/6
    return t, y