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
        sum += f(a+k*h)
    
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
        sum += 4 * f(a+k*h)
        sum += 2 * f(a+k*h)

    return sum

"""
   This module is the hub for the library

   For details and applications of these functions, check Documentation on GitHub
"""
