#!/usr/bin/env python
#
# File: Corana.py
# Corana's parabola -- multimodal paraboloid function
#
# Usage: Corana.py <parm1> .. <parm4>
import sys, math
from DeDriver import *

class corana:
    _LARGEVAL = 1.0E20
    def __init__(self, fname='de_fcn.dat'):
        """ Setup parms for the funciton        
        """
        self.D = 4
        self.x = [corana._LARGEVAL]*self.D
        self.outfname = fname
        self.f = corana._LARGEVAL


    def setVect(self, x):
        """ Set the parm vector from the command line args
        """
        if len(x) != self.D:
            return False
        for j in range(self.D):
            self.x[j] = float(x[j])
        return True


    def sgn(self, arg):
        """ Return the sign of the argument
        """
        if arg == 0.0:
            return 0.0
        if arg < 0.0:
            return -1.0
        if arg > 0.0:
            return 1.0
        

    def zj(self, xj):
        return math.floor(abs(xj / 0.2) + 0.49999) * self.sgn(xj) * 0.2


    def calc(self):
##     def calc(self, x):
        """ Calculate the cost value
        """
        #--->
##         if not self.setVect(x):
##             self.f = corana._LARGEVAL
##             return self.f
        self.f = corana._LARGEVAL
##         self.x = x
        try:
            #---^
            d = (1, 1000, 10, 100)
            sum = 0.0
            for j in range(self.D):
                zj = self.zj(self.x[j])
                if abs(self.x[j] - zj) < 0.05:
                    sum += 0.15 * (zj - 0.05 * self.sgn(zj))**2 * d[j]
                else:
                    sum += d[j] * self.x[j]**2
            self.f = sum
        except:
            pass
        #--->
        return [self.f]


    def __call__(self, x, args=None):
        if self.setVect(x):
            return self.calc()
        else:
            return corana._LARGEVAL

    def saveOutput(self):
        """ Save the output to a file.
	"""
	try:
	    file = open(self.outfname,'w')
        except:
	    return False
        else:
	    file.write(`self.f`)
            file.close()
	    return True



def deMinDrive():
    c = corana()
    ir = [(0, -1e4, 1e4), (1, -1e4, 1e4),
          (2, -1e4, 1e4), (3, -1e4, 1e4)]
    dim = 4
    fdim = 1
    #np = dim*6
    np = dim*10
    softConstr = []
    hardConstr = ['-x[0]-1e4', 'x[0]-1e4',
                  '-x[1]-1e4', 'x[1]-1e4',
                  '-x[2]-1e4', 'x[2]-1e4',
                  '-x[3]-1e4', 'x[3]-1e4']
    # cxi = 0.25
    cxi = 0.15
    # cxi = 0.10
    #
    # cmt = 0.05
    cmt = 0.10
    gmin = 5
##     gmax = 5000
    gmax = 200
    tol = 1e-6
    dtol = 5e-7
    m = 5
    truncFrac = 0.4
    de = DeDriver(dim, fdim, np)
    de.local = True
    de.costFcn = c
    de.setConstraints(softConstr, hardConstr)
    de.setDEStrategy(xi=cxi, mt=cmt)
    de.setConvergenceStrategy('fractol', gmin, gmax, tol, dtol, m, truncFrac)
    de.ctrlPop.mt = 0.05
    # de.ctrlPop.mt = 0.1
    de.initialize(ir)
    print `de`
    print
    # de.sade = False
    de.sade = True
    # de.zeta = 0.5
    de.zeta = 0.2
##     de.zeta = 0.7
    minCost = de.minimize()
    print 'prevMGenStats:'
    for i in range(len(de.prevMGenStats)):
        print
        for j in range(len(de.prevMGenStats[0])):
            print '%.6E ' % de.prevMGenStats[i][j],
    print 'maxDeltaFracTol = %.6E' % de.maxDeltaStats[4]
    print 
    print `de`


if __name__ == "__main__":
##     c = corana('de_fcn.dat')
##     if not c.setVect(sys.argv[1:]):
##         sys.exit(1)
##     c.calc()
##     if not c.saveOutput():
##         sys.exit(1)
    deMinDrive()
    

"""
Test Results:

Without SA (max number of generations required for convergence):
>>> gmaxNoSA = [103, 99, 102, 106, 107, 92, 98, 108, 99, 106]
>>> mean(gmaxNoSA)
102.0
>>> std(gmaxNoSA)
4.7749345545253288


WITH SA:
xi = 0.1
mt = 0.10
zeta = 0.5
>>> gmax = [68, 75, 82, 86, 74, 105, 78, 78, 86, 88]
>>> stats.mean(gmax)
82.0
>>> std(gmax)
9.6850400102426022


** 19.6% relative gain in the number of generations required for convergence.

xi = 0.1
mt = 0.05
zeta = 0.5
>>> gmaxSAmt5 = [60, 64, 90, 93, 117, 67, 83, 95, 64, 75]
>>> mean(gmaxSAmt5)
80.799999999999997
>>> std(gmaxSAmt5)
17.238329385413195

Removing the outlier (117):
>>> s = [60, 64, 90, 93, 67, 83, 95, 64, 75]
>>> mean(s)
76.777777777777771
>>> std(s)
12.97671228502794


xi = 0.1
mt = 0.05
zeta = 0.2
>>> gz02 = [101, 94, 75, 82, 132, 139, 84, 101, 139, 74]
>>> mean(gz02)
102.09999999999999
>>> std(gz02)
24.373961516339524


xi = 0.1
mt = 0.05
zeta = 0.7
>>> gz07 = [91, 111, 62, 82, 75, 87, 73, 63, 99, 75]
>>> mean(gz07)
81.799999999999997
>>> std(gz07)
14.749915253993834

xi = 0.25
mt = 0.05
zeta = 0.5
>>> gxi25 = [95, 89, 116, 124, 77, 87, 113, 121, 104, 90]
>>> mean(gxi25)
101.59999999999999
>>> std(gxi25)
15.415576538034509

** Changed ctrlPop mutation scale factor to sample from a Breit-Wigner
** pdf (mu=0.5, sigma=0.1) rather than Gaussian (mu=0.5, sigma=0.2)

>>> g
[81, 91, 96, 87, 105, 74, 75, 94, 85, 100]
>>> mean(g)
88.799999999999997
>>> std(g)
9.7959175170067656
>>> sigma
[0.00189449325507, 0.0071933742942199997, 0.039728248135599997, 0.016281945507699999, 0.01031007687, 0.0094201534354699992, 0.012768007369799999, 0.00543563193602, 0.012605486114400001, 0.0074389628909599997]
>>> mean(sigma)
0.012307637980923998
>>> ftol
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

"""
