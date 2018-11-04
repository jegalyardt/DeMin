#!/usr/bin/env python
#
# Benchmark problems for optimization
#
#~~~~~
import numpy as np
import scipy.optimize as sopt
from DeDriver import DeDriver

_LARGEVAL = 1e100

# Shekel data, from ICEO1 competition
shekel_a = np.array([[9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020],
                     [9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374],
                     [8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982],
                     [2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426],
                     [8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567],
                     [7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208],
                     [1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448],
                     [8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762],
                     [0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637],
                     [7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247],
                     [0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016],
                     [2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789],
                     [8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109],
                     [2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564],
                     [4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670],
                     [8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826],
                     [8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591],
                     [4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740],
                     [2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675],
                     [6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258],
                     [0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070],
                     [5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234],
                     [3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027],
                     [8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064],
                     [1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224],
                     [0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644],
                     [0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229],
                     [4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506],
                     [9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732],
                     [4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500]])

shekel_c = np.array([0.806,
                     0.517,
                     0.1,
                     0.908,
                     0.965,
                     0.669,
                     0.524,
                     0.902,
                     0.531,
                     0.876,
                     0.462,
                     0.491,
                     0.463,
                     0.714,
                     0.352,
                     0.869,
                     0.813,
                     0.811,
                     0.828,
                     0.964,
                     0.789,
                     0.360,
                     0.369,
                     0.992,
                     0.332,
                     0.817,
                     0.632,
                     0.883,
                     0.608,
                     0.326])

langerman_a = np.array([[9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020],
	[9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374],
	[8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982],
	[2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426],
	[8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567],
	[7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208],
	[1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448],
	[8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762],
	[0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637],
	[7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247],
	[0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016],
	[2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789],
	[8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109],
	[2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564],
	[4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670],
	[8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826],
	[8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591],
	[4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740],
	[2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675],
	[6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258],
	[0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070],
	[5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234],
	[3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027],
	[8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064],
	[1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224],
	[0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644],
	[0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229],
	[4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506],
	[9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732],
	[4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500]])

langerman_c = np.array([0.806,
                        0.517,
                        1.5,
                        0.908,
                        0.965,
                        0.669,
                        0.524,
                        0.902,
                        0.531,
                        0.876,
                        0.462,
                        0.491,
                        0.463,
                        0.714,
                        0.352,
                        0.869,
                        0.813,
                        0.811,
                        0.828,
                        0.964,
                        0.789,
                        0.360,
                        0.369,
                        0.992,
                        0.332,
                        0.817,
                        0.632,
                        0.883,
                        0.608,
                        0.326])

"""
valuestoreach:
Here are the putative global minima for the ICEO1 test functions.
For Langerman's functions, I could only retrieve the value to reach
specified by the contest. If someone has the global minimum, 
please let me know at neum@cma.univie.ac.at

Shekel's Fox Holes (5 Dimensions)
fmin = -10.40395206000838

Shekel's Fox Holes (10 Dimensions)
fmin = -10.20787684027058

Michalewicz's Function (5 Dimensions)
fmin = -4.68765817908815

Michalewicz's Function (10 Dimensions)
fmin = -9.66015167547059

Langerman's Function (5 Dimensions)
Value to reach in the ICEO contest: -1.4

Langerman's Function (10 Dimensions)
Value to reach in the ICEO contest: -1.4
"""

def Sphere(x, n):
    s = 0.0
    for i in xrange(n):
        s += (x[i] - 1.0) * (x[i] - 1.0)
    #return s
    return [s]

griewangk_D = 4000.0
def Griewangk(x, n):
    f = _LARGEVAL
    # try:
    summ = 0.0
    prod = 1.0
    for j in xrange(n):
        summ += (x[j] - 100.0) * (x[j] - 100.0)
        prod *= np.cos(x[j]/np.sqrt(j+1.0))
    f = summ / griewangk_D - prod + 1.0
    # except:
    #     pass
    # return f
    return [f]


def Shekel(x, n):
    """ Shekel's Foxholes, n <= 10.
    """
    sp = 0.0
    result = 0.0
    h = 0.0
    for i in xrange(30):
        sp = 0.0
        for j in xrange(n):
            h = x[j] - shekel_a[i,j]
            sp += h * h
        result += 1.0 / (sp + shekel_c[i])
    # return -result
    return [-result]

micha_m = 10.0

def Micha(x, n): #, m=10.0):
    """ Michalewitz's function.
    """
    u = 0.0
    for i in xrange(n):
        u += np.sin(x[i]) * (np.sin((i+1.0)*x[i]*x[i]/np.pi))**(2.0*micha_m)
    # return -u
    return [-u]


def SqrDst(x1, x2, n):
    """ For use with Langerman's function
    """
    dist = 0.0
    d = 0.0
    for i in xrange(n):
        d = x1[i] - x2[i]
        dist += d*d
    return dist


def Langerman(x,n): # mn=500.0
    """ Langerman's function
    """
    s = 0.0
    for i in xrange(n):
        dist = SqrDst(x, langerman_a[i], n)
        s -= langerman_c[i] * ( np.exp(-dist/np.pi) * np.cos(np.pi * dist) )
    # return s
    return [s]


def PittyParabola(x,n):
    return ( np.cos(14.5 * x[0] - 0.3) + (x[0] + 0.2) * x[0] 
             + (x[1] + 0.2) * x[1] )


problems = {0:{'fcn':Sphere, 'dim':5, 'ir':(-10.0,10.0), 'vtr':0.0},
            1:{'fcn':Sphere, 'dim':10, 'ir':(-10.0,10.0), 'vtr':0.0},
            2:{'fcn':Griewangk, 'dim':5, 'ir':(-400.0,400.0), 'vtr':0.0},
            3:{'fcn':Griewangk, 'dim':10, 'ir':(-400.0,400.0), 'vtr':0.0},
            4:{'fcn':Shekel, 'dim':5, 'ir':(-100.0,100.0), 'vtr':-10.40395206000838},
            5:{'fcn':Shekel, 'dim':10, 'ir':(-100.0,100.0), 'vtr':-10.20787684027058},
            6:{'fcn':Micha, 'dim':5, 'ir':(-100.0,100.0), 'vtr':-4.68765817908815},
            7:{'fcn':Micha, 'dim':10, 'ir':(-100.0,100.0), 'vtr':-9.66015167547059},
            8:{'fcn':Langerman, 'dim':5, 'ir':(-100.0,100.0), 'vtr':-1.4},
            9:{'fcn':Langerman, 'dim':10, 'ir':(-100.0,100.0), 'vtr':-1.4},
            10:{'fcn':PittyParabola, 'dim':2, 'ir':(-100.0,100.0), 'vtr':-999.0}}


def OptBasinHop(problemID, niter=200, x0=np.zeros(2), disp=False, niter_success=10):
    """ Optimize the benchmark functions via the Basin Hopping algorithm.
    """
    args = {"method": "Powell", "args":(problems[problemID]['dim'],)}
    #x0 = 
    ret = sopt.basinhopping(problems[problemID]['fcn'], x0, niter=niter, disp=disp,
                            minimizer_kwargs=args, niter_success=niter_success)
    print 'global minimum: x = [%.4f, %.4f], f(x0) = %.4f' % (ret.x[0],
                                                              ret.x[1],
                                                              ret.fun)
    print '  Number of function evals = %d' % (ret.nfev)
    return None


def BenchDeDriver(problemID, evalBudgetFac=100, npfac=4, xi=0.2, mt=0.05, 
                  ctrl_mt=0.05, zeta=0.15):
    print '\n=========================================================\n'
    print 'Problem ID: %d' % (problemID)
    print '\n=========================================================\n'
    dim = problems[problemID]['dim']
    evalBudget = dim * evalBudgetFac
    # np = dim * 4
    np = dim * npfac
    print 'np = %d' % (np)
    #gmax = 500
    gmax = int(evalBudget / float(np))
    print 'gmax = %d' % (gmax)
    gmin = 5
    # Optimization strategy transition point (exploration first, 
    # then exploitation);
    # 'xi' is a fraction of gmax at which point the strategy switches.
    #cxi = 0.20
    #
    # Initial range for each parameter under optimization
    ir = []

    # Hard constraints: each parameter should lie in [0.0, 1.0]
    lo = problems[problemID]['ir'][0]
    hi = problems[problemID]['ir'][1]
    softConstr = []
    hardConstr = []
    for i in xrange(dim):
            ir.append((i, lo, hi))
            hardConstr.append('-x[%d]+%.1f' % (i, lo))
            hardConstr.append('x[%d]-%.1f' % (i, hi))
    # tolerance for convergence test
    # tol = 1e-10
    # tol = 1e-20
    # tol = 1e-100
    tol = 1e-300
    # delta tolerance btwn the best truncFrac*np members
    dtol = 0.5 * tol
    # Must meet convergence criteria for 'm' generations
    m = gmax
    truncFrac = 0.25
    fdim = 1
    de = DeDriver(dim, fdim, np)
    de.costArgs = dim
    de.local = True
    de.costFcn = problems[problemID]['fcn']
    de.setConstraints(softConstr, hardConstr)
    #de.setDEStrategy(xi=xi, mt=mt, ph2StratStr='best-trig1')
    #de.setDEStrategy(xi=xi, mt=mt, ph2StratStr='best2')
    de.setDEStrategy(xi=xi, mt=mt, ph1StratStr='rand1', ph2StratStr='best-trig1')
    #de.setDEStrategy(xi=xi, mt=mt, ph1StratStr='rand1', ph2StratStr='best1')
    #de.setDEStrategy(xi=xi, mt=mt, ph1StratStr='rand1', ph2StratStr='best2')
    de.setConvergenceStrategy('fractol', gmin, gmax, tol, dtol, m, truncFrac)
    # Use self-adaptive algorithm
    # de.saveInitPop = True
    de.sade = True
    # Control population evolution parameter
    de.ctrlPop.mt = ctrl_mt
    #de.ctrlPop.strategy = ('best-trig',1)
    de.ctrlPop.strategy = ('rand-trig',1)
    #de.ctrlPop.strategy = ('rand',1)
    # zeta - percentage of Pareto front to use when selecting DE control 
    #        parameter vector (F, Cr).
    #de.zeta = 0.15
    de.zeta = zeta
    #de.ctrlPop.strategy = ('rand',1)
    de.initialize(ir)
    print `de`
    print
    minCost = de.minimize()
    print 'maxDeltaFracTol = %.6E' % de.maxDeltaStats[4]
    print 
    print `de`
    print
    print 'Value_to_reach - Best cost = %.4e' % (problems[problemID]['vtr']-de.population[de.population.ibest].y[0])
    return None



if __name__ == "__main__":
    # probID = 5 # PittyParabola
    # niter = 200
    # #x0 = np.zeros(2)
    # x0 = np.random.random(2) * 100.0

    # probID = 10
    # niter = 200
    # x0 = np.random.random(problems[probID]['dim'])
    # OptBasinHop(probID, niter, x0, disp=False)

    # probID = 2
    # xi = 0.2
    # #xi = 0.10
    # # mt = 0.05
    # mt = 0.2
    # # zeta = 0.4
    # zeta = 0.1
    # #zeta = 1.0
    # # zeta = 0.75
    # #ctrl_mt = 0.05
    # #ctrl_mt = 0.5
    # ctrl_mt = 0.3
    # npfac = 4
    # #budgetFac = 1e4
    # budgetFac = 100

    probID = 4
    # xi = 0.10
    # xi = 0.2
    # xi = 0.3
    # xi = 0.4
    # xi = 0.5
    # xi = 0.6
    xi = 0.7
    # xi = 0.8
    # mt = 0.05
    mt = 0.2
    # mt = 0.3
    # mt = 0.5
    zeta = 0.05
    #zeta = 0.1
    #zeta = 0.2
    # zeta = 0.3
    # zeta = 0.4
    #zeta = 0.2
    #zeta = 1.0
    # zeta = 0.75
    #ctrl_mt = 0.05
    #ctrl_mt = 0.10
    ctrl_mt = 0.20
    #ctrl_mt = 0.4
    #ctrl_mt = 0.3
    #npfac = 4
    #npfac = 5
    npfac = 6
    #npfac = 10
    budgetFac = 1e4
    # budgetFac = 500

    # probID = 5
    # # xi = 0.10
    # # xi = 0.2
    # # xi = 0.3
    # # xi = 0.4
    # # xi = 0.5
    # # xi = 0.7
    # xi = 0.8
    # # mt = 0.05
    # mt = 0.2
    # #zeta = 0.1
    # zeta = 0.2
    # #ctrl_mt = 0.05
    # ctrl_mt = 0.10
    # #npfac = 4
    # #npfac = 5
    # npfac = 6
    # # budgetFac = 1e4
    # budgetFac = 500

    # probID = 6
    # # xi = 0.10
    # # xi = 0.2
    # # xi = 0.3
    # # xi = 0.4
    # xi = 0.5
    # # xi = 0.7
    # # xi = 0.8
    # # mt = 0.05
    # mt = 0.2
    # # zeta = 0.1
    # # zeta = 0.2
    # zeta = 0.3
    # #ctrl_mt = 0.05
    # ctrl_mt = 0.10
    # #npfac = 4
    # #npfac = 5
    # npfac = 6
    # #budgetFac = 1e4
    # budgetFac = 500

    # probID = 7
    # # xi = 0.10
    # # xi = 0.2
    # # xi = 0.3
    # # xi = 0.4
    # # xi = 0.5
    # xi = 0.6
    # # xi = 0.7
    # # xi = 0.8
    # # mt = 0.05
    # mt = 0.2
    # # zeta = 0.1
    # # zeta = 0.2
    # zeta = 0.3
    # #ctrl_mt = 0.05
    # ctrl_mt = 0.10
    # #ctrl_mt = 0.20
    # #ctrl_mt = 0.50
    # #npfac = 4
    # #npfac = 5
    # npfac = 6
    # #budgetFac = 1e4
    # #budgetFac = 100
    # budgetFac = 500

    # probID = 8
    # # xi = 0.10
    # # xi = 0.2
    # # xi = 0.3
    # # xi = 0.4
    # # xi = 0.5
    # xi = 0.7
    # # xi = 0.8
    # # mt = 0.05
    # mt = 0.2
    # # zeta = 0.1
    # # zeta = 0.2
    # zeta = 0.3
    # #ctrl_mt = 0.05
    # ctrl_mt = 0.10
    # #npfac = 4
    # #npfac = 5
    # npfac = 6
    # budgetFac = 1e4
    # #budgetFac = 100

    # probID = 9
    # # xi = 0.10
    # # xi = 0.2
    # # xi = 0.3
    # # xi = 0.4
    # xi = 0.5
    # # xi = 0.7
    # # xi = 0.8
    # # mt = 0.05
    # mt = 0.2
    # # zeta = 0.1
    # # zeta = 0.2
    # zeta = 0.3
    # #ctrl_mt = 0.05
    # ctrl_mt = 0.10
    # #npfac = 4
    # #npfac = 5
    # npfac = 6
    # budgetFac = 1e4
    # #budgetFac = 100

    BenchDeDriver(probID, budgetFac, npfac=npfac, xi=xi, zeta=zeta, 
                  mt=mt, ctrl_mt=ctrl_mt)

