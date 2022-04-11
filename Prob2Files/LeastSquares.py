import numpy as np
import matplotlib.pyplot as pyplot
from math import *

class LeastSquaresFit_Class():
    def __init__(self, xdata=None, ydata=None):
        self.x=xdata if xdata is not None else np.array([])
        self.y=ydata if ydata is not None else np.array([])
        self.coeffs=np.array([])

    def RSquared(self, a):
        '''
        To calculate the R**2 value for a set of x,y data and a LeastSquares fit with polynomial having coefficients a
        :param x: array of actual x data
        :param y: array of actual y data
        :param a:  the coefficients of the polynomial fit
        :return:
        '''
        AvgY=np.mean(self.y) #calculates the average value of y
        SSTot=0
        SSRes=0
        for i in range(len(self.y)):
            SSTot+=(self.y[i]-AvgY)**2
            SSRes+=(self.y[i]-self.Poly(self.x[i],a))**2
        RSq=1-SSRes/SSTot
        return RSq

    def Poly(self,xval, a):
        """
        calculates the value for a polynomial given a value for x and the coefficients of the polynomial.
        f(x)=y=a[0]*x**(n-1)+a[1]*x**(n-2)+...+a[n], where n=len(a)
        :param a:  the coefficients for an n-1 order polynomial it
        :return: the value of the fit at xval
        """
        p=np.poly1d(a)
        return p(xval)

    def LeastSquares(self, power):
        """
        Uses polyfit from numpy.  See that documentation.
        :param power:
        :return:
        """
        self.coeffs=np.polyfit(self.x, self.y, power)
        return self.coeffs

    def GetCoeffsString(self):
        """
        Get output of coefficients as a formatted string.
        :return: a formatted string of the coefficients that is comma delimited.
        """
        s=''
        n=0
        for c in self.coeffs:
            s += ('' if n == 0 else ', ')+"{:0.4f}".format(c)
            n += 1
        return s

    def GetPlotInfo(self, power, npoints=500):
        Xmin = min(self.x)
        Xmax = max(self.x)
        Ymin = min(self.y)
        Ymax = max(self.y)
        dX = 1.0 * (Xmax - Xmin) / npoints

        a = self.LeastSquares(power)

        xvals = []
        yvals = []
        for i in range(npoints):
            xvals.append(Xmin + i * dX)
            yvals.append(self.Poly(xvals[i], a))
        RSq = self.RSquared(a)
        return xvals,yvals,RSq
