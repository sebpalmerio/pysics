class Integration:

    """
    A Class for creating instances of functions we wish to integrate. The class is made up of a function attribute f which stores a function we wish to integrate.
    Attributes a and b are for our bounds of integration on an interval (a,b) (i.e. a<b). The current implementation also takes in a parameter N for the number of partitions for a given integration scheme (i.e. how many slices/terms in the sumation). 
    The added attribute h is the step size used when calculating the area of our differential slices. 

    Integration: function, float, float, int

    Methods:

        trapezoid:
            Implements the trapezoid integration scheme as a method returning the value of the integral. 

            Arguments:
                self: takes in itself and returns the value of the integral 

            trapezoid: self -> float


        simpsons:
            Implements the Simpson's integration scheme as a method returning the value of the integral. 

            Arguments:
                self: takes in itself and returns the value of the integral  
       
            simpsons: self -> float
    """

    #def __init__(self, f, a, b, N):
     #   self.a = a
      #  self.b = b
       # self.N = N
        #self.h = (self.b-self.a)/self.N
        #self.f = f
    
    def trapezoid(self):
    
        total = (self.f(self.a)+self.f(self.b))/2
        for i in range(1,self.N):
            total += self.f(self.a+i*self.h)
    
        return total

    def simpsons(self):
      
        total = (self.f(self.a)+self.f(self.b))*(self.h/3)
        for i in range(1,self.N):
            total += 4 * self.f(self.a+i*self.h)
            total += 2 * self.f(self.a+i*self.h)

        return total

class Ode:
    #__init__(self, ti, tf, h, yinit):
     #   self.ti,self.tf, self.h, self.yinit = ti, tf, h, yinit

    def RK4(self):
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
    N_steps = int((self.tf-self.ti)/self.h) # number of steps
    t = self.h * arange(N_steps) # grid for the independent variable   
    y = empty([len(self.yinit), len(t)]) #(n, N_step) array for the dependent variable
    ys = array(self.yinit) # n-dim buffer array, set to the initial conditions

    for i, ts in enumerate(t): # RK4 integrator
        y[:,i] = ys
        k1 = self.h*f(ts, ys)
        k2 = self.h*f(ts + 0.5*self.h, ys + 0.5*k1)
        k3 = self.h*f(ts + 0.5*self.h, ys + 0.5*k2)
        k4 = self.h*f(ts + self.h, ys + k3)
        ys += (k1 + 2*k2 + 2*k3 + k4)/6
    return t, y

class Function(Ode, Integration):
    __init__(self, yinit, a = 0, b = 10, N = 1000, ti = 0.0, tf = 20.0, h = 0.001):
        self.yinit, self.a, self.b, self.ti, self.tf = yinit, a, b, ti, tf
        self.N, self.h = N, h 



        