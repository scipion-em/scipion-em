import math
import numpy as np

class FibonacciSphere:
    def __init__(self, half_points=2500):

        self.half_points= half_points
        self.num_points = self.half_points*2

        self.sphX, self.sphY , self.sphZ, self.ga = FibonacciSphere.fibonacci_sphere(self.num_points)

        self.searchRange = round(180. /self.ga)+1

        self.weights = np.zeros(len(self.sphX))

    def append(self, x,y,z):

        # z is a linear space array from -1 to 1
        # In a simple case for num_points= 100
        # z = [1/100-1, ... 1/100+1]

        # z index
        zIndex = round((z*self.half_points)-1 +(self.half_points)) # + 1 
        numPoints = self.num_points-1
        if zIndex < 0:
            zIndex = 0
        elif zIndex > numPoints:
            zIndex = numPoints

        zMin = max(zIndex-self.searchRange, 0)
        zMax = min(zIndex + self.searchRange, numPoints)

        distance = 4 # way maximun distance
        finalIndex = -1

        # x, y , go through the sphere
        for index in range(zMin,zMax):
            xSph = self.sphX[index]
            ySph = self.sphY[index]

            newDist = (xSph-x)**2 + (ySph-y)**2

            if newDist < distance:
                finalIndex = index
                distance = newDist

        # Add one to the final weight
        self.weights[finalIndex] = self.weights[finalIndex]+1

    def cleanWeights(self):

        X, Y, Z, W = [],[],[],[]

        for x,y,z,w in zip(self.sphX, self.sphY, self.sphZ, self.weights):
            if w != 0:
                X.append(x)
                Y.append(y)
                Z.append(z)
                W.append(w)

        self.sphX = X
        self.sphY = Y
        self.sphZ = Z
        self.weights = W

    @staticmethod
    def fibonacci_sphere(num_points: int = 5000):
        """
        1000 OK small graph
        5000 OK large graph
        more than 10000 does not make sense

        Computes a set of x,y,z equidistant points on the surface of a sphere of radius 1 using the fibonacci method

        :param num_points: number of points to compute.
        :returns x,y,z as lists

        """
        ga = (3 - np.sqrt(5)) * np.pi  # golden angle

        # Create a list of golden angle increments along tha range of number of points
        theta = ga * np.arange(num_points)

        # Z is a split into a range of -1 to 1 in order to create a unit circle
        z = np.linspace(1 / num_points - 1, 1 - 1 / num_points, num_points)

        # a list of the radii at each height step of the unit circle
        radius = np.sqrt(1 - z * z)

        # Determine where xy fall on the sphere, given the azimuthal and polar angles
        y = radius * np.sin(theta)
        x = radius * np.cos(theta)
        return x, y, z, ga


class TrigonometricMemoization(object):

    """
    Stores sines of 0 to 90 degrees in a list.
    These are used as the basis for calculating and
    returning sines, cosines and tangents of any angle.

    With the default, errors are in the 4th decimal point
    """
    def __init__(self, factor=100):

        """
        Calculates sines of 0 to 90 degrees and appends them to a list.

        :param factor: default(100). List size = 90 * factor. The higher the factor the higher the precision.

        """
        self.memoized_sines = []
        degrees_in_radian = (factor * 180.0) / math.pi
        self.factor = factor
        maximum = 90 * factor + 1
        for d in range(0, maximum):  # samble 1/10 degree
                                 # we may try 1/100 -> 9001
            self.memoized_sines.append(math.sin(d / degrees_in_radian))


    def sine(self, degrees):
        """
        Calculates sine of any angle from values stored in list.
        """
        sine_negative = False
        # get an angle between -360 and 360 using modulus operator
        degrees_shifted = degrees % 360

        # move negative angles by 360 into the positive range
        if(degrees_shifted < 0):
            degrees_shifted += 360

        # if angle is in the 180-360 range shift it to the 0-180 range
        # and record that we need to negate the sine
        if(degrees_shifted > 180):
            degrees_shifted-= 180
            sine_negative = True

        # if degrees_shifted is 90-180 we need to "mirror" it into the 0-90 range
        if(degrees_shifted > 90 and degrees_shifted <= 180):
            degrees_shifted = 90 - (degrees_shifted - 90)

        # now get the sine from the memoized lookup table
        sine = self.memoized_sines[round(degrees_shifted * self.factor)]

        # negate if angle was 180-360
        if(sine_negative):
            sine*= -1

        return sine


    def cosine(self, degrees):
        """
        Cosine can easily be calculated from sine by shifting 90 degrees.
        """
        return self.sine(degrees + 90)

    def tangent(self, degrees):
        """
        Tangent is calculated using the sine and cosine from the other two methods.
        """
        try:
            sine = self.sine(degrees)
            cosine = self.cosine(degrees)
            return sine / cosine
        except:
            return math.nan

