#!/usr/bin/env python
#
# File: Katsuura.py
# Author: J. Galyardt
# Description:
#   This is an implimentation of Katsuura's function
#   for use in testing the DeDriver.py minimization tool.
#   GLOBAL MINIMUM: f(0) = 1.0
#
# Usage: Kasuura.py <D> <beta> <parm 1> <parm 2> ... <parm D>

import sys, math
from DeDriver import *

#from DeMin import *


class katsuura:
    _LARGEVAL = 1.0E20
    def __init__(self, D=10, beta=32, fname='out.dat'):
        """ Setup the parameters for the function.
        """
        self.D = D
        self.beta = beta
        self.x = [katsuura._LARGEVAL]*self.D
        self.outfname = fname
        self.f = katsuura._LARGEVAL


    def setVect(self, x):
        """ Set the parameter vector from the command-line args.
        """
        if len(x) != self.D:
            return False
        for j in range(self.D):
            self.x[j] = float(x[j])
        return True


    def calc(self, x, args=None):
        """ Calculate Katsuura's function.
        """
        self.f = katsuura._LARGEVAL
        self.x = x
        try:
            f = 1.0
            for j in range(self.D):
                sum = 0.0
                for k in range(self.beta+1):
                    sum += (abs(2**k * self.x[j] - round(2.0**k * self.x[j])) *
                            2**(-k))
                f *= 1.0 + (j + 1.0) * sum
            self.f = f
        except:
            pass
        return self.f
    

    def saveOutput(self):
        """ Save the output to a file.
        """
        try:
            file = open(self.outfname, 'w')
        except:
            return False
        else:
            file.write(`self.f`)
            return True


def driveDeMin(D=10, beta=32):
    k = katsuura(D, beta)
    dim = D
    np = dim * 6
    lo = -10.0
    hi = 10.0
    ir = []
    softConstr = []
    hardConstr = []
    for i in range(dim):
        ir.append((i, lo, hi))
        hardConstr.append('-x['+`i`+']+'+`lo`)
        hardConstr.append('x['+`i`+']-'+`hi`)
    xi = 0.6
    gmin = 5
    gmax = 10000
    tol = 1e-6
    dtol = 5e-7
    m = 5
    #truncFrac = 0.4
    truncFrac = 0.2
    # de = DeMin(k.calc, dim, np, -1, ir)
    fdim = 1
    de = DeDriver(dim, fdim, np)
    de.local = True
    de.sade = True
    de.costFcn = k.calc
    de.zeta = 0.5
    # de.zeta = 0.2
    # de.zeta = 0.7
    de.setConstraints(softConstr, hardConstr)
    de.setDEStrategy(xi)
    de.setConvergenceStrategy('fractol', gmin, gmax, tol, dtol, m, truncFrac)
    de.ctrlPop.mt = 0.05
    #de.initialize()
    de.initialize(ir)
    print `de`
    print
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
    # kat.D should be the dimensionality of the DE vectors
    # (Member objects)
    if len(sys.argv) != 3:
        print 'Usage: Katsurra.py <D> <beta>'
        sys.exit(1)
    D = int(sys.argv[1])
    if D < 0:
        sys.exit(1)
    beta = int(sys.argv[2])
    if beta < 0:
        sys.exit(1)
##     kat = katsuura(D, beta, 'de_fcn.dat')
##     if not kat.setVect(sys.argv[3:]):
##         sys.exit(1)
##     kat.calc()
##     if not kat.saveOutput():
##         sys.exit(1)
    print 'Katsuura.py: using D = %d, beta = %d\n' % (D, beta)
    driveDeMin(D, beta)
    
