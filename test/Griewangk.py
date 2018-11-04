#!/usr/bin/env python
#
# Griewangk.py
#
# Griewangk's function:
# Sum(j=1,D) x_j**2/4000 - Prod(j=1,D) cos(x_j/sqrt(j)) + 1
#
import sys, math
from DeDriver import *

class griewangk:
    _LARGEVAL=1e20
    def __init__(self, D=100):
        self.D = D
        self.x = [griewangk._LARGEVAL]*self.D
        self.f = [griewangk._LARGEVAL]


    def calc(self, x, args=None):
        self.x = x
        self.f = [griewangk._LARGEVAL]
        try:
            summ = 0.0
            prod = 1.0
            for j in xrange(self.D):
                summ += self.x[j]**2/4000.0
##                 prod *= math.cos(self.x[j]/math.sqrt(j+1.0)) + 1.0
                prod *= math.cos(self.x[j]/math.sqrt(j+1.0))
##             self.f = summ - prod 
            self.f = [summ - prod + 1.0]
        except:
            pass
        return self.f


def driveDeMin(D):
    g = griewangk(D)
    dim = D
    # np = dim * 4
    np = dim * 3
    # np = dim * 6
##     lo = -600.0
##     hi = 600.0
    lo = -400.0
    hi = 400.0
    ir = []
    softConstr = []
    hardConstr = []
    for i in range(dim):
        ir.append((i, lo, hi))
        hardConstr.append('-x['+`i`+']+'+`lo`)
        hardConstr.append('x['+`i`+']-'+`hi`)
    # cxi = 0.1
    cxi = 0.15
    # cxi = 0.2
##     cxi = 0.3
    gmin = 5
    gmax = 1000
    tol = 1e-7
    dtol = 5e-8
    m = 5
    truncFrac = 0.25
##     de = DeMin(g.calc, dim, np, -1, ir)
    fdim = 1
    de = DeDriver(dim, fdim, np)
    de.local = True
    de.costFcn = g.calc
    de.setConstraints(softConstr, hardConstr)
##     de.setDEStrategy(xi=cxi, mt=0.05, ph2HDCr=0.5, ph2F=0.6)
##     de.setDEStrategy(xi=cxi, mt=0.1, ph1HDCr=0.1, ph1F=0.5, ph2HDCr=0.1, ph2F=0.5)
    de.setDEStrategy(xi=cxi, mt=0.05, ph2StratStr='best-trig1')
    de.setConvergenceStrategy('fractol', gmin, gmax, tol, dtol, m, truncFrac)
    de.sade = True
    de.ctrlPop.mt = 0.05
##     de.sade = False
    # de.zeta = 0.2
    de.ctrlEps = 5E-2
    # de.zeta = 0.3
    de.zeta = 0.1
##     de.zeta = 0.5
    de.initialize(ir)
    print `de`
    print
    minCost = de.minimize()
    print 'maxDeltaFracTol = %.6E' % de.maxDeltaStats[4]
    print 
    print `de`
    

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: Griewangk.py <D>'
        sys.exit(1)
    D = int(sys.argv[1])
    if D <= 0:
        sys.exit(1)
    driveDeMin(D)
    

"""
D = 10
SADE:
>>> gSA
[495, 640, 403, 758]
>>> sigmaSA
[3.7187044869299999e-11, 1.03508919563e-16, 1.09553641431e-16, 7.8252601022399998e-17]
>>> ftolSA
[7.7907509016499994e-08, 0.0, 0.0, 0.0]
>>> saTrials
10

Notes: saTrials := number of trials attempted.  Difference btwn
saTrials and len(gSA) is the number of trials which failed to converge
to the global minimum.

Standard DeDriver:
>>> g
[393]
>>> sigma
[7.5372405479399995e-17]
>>> ftol
[0.0]
>>> trials
10


SADE with ctrlPop scale factor randomly generated from a cauchy pdf (peak: 0.5, fwhm: 0.1)
>>> g
[433, 376, 641, 669, 632, 488]
>>> mean(g)
539.83333333333337
>>> std(g)
112.80871223250249
>>> sigma
[8.3337864601700004e-17, 1.60886601221e-16, 0.0035763011385000001, 8.2957758334900004e-17, 8.2623730916499998e-17, 6.5561273426600005e-17]
>>> mean(sigma)
0.00059605018975007924
>>> ftol
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


"""
