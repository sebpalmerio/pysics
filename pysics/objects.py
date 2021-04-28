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

    def __init__(self, f, a, b, N):
        self.a = a
        self.b = b
        self.N = N
        self.h = (self.b-self.a)/self.N
        self.f = f
    
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