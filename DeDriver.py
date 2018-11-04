#!/usr/bin/env python
#
# TO-DO:
# * Cleanup command-line parser code -- use callback functions!
#   -> esp. for the save/load initial file options.
# * Resume the DeDriver alg after a pause from a specified generation
#   (especially useful if errors occur and no CAF output is available)
# * Clean up error messages, log file messages
# * Improve Local running such that downloading the function eval
#   tarball is only done for generation 0, then subsequently just
#   copied to the run directory.
# * Improve Error handling!
#   - Esp. w.r.t. exception data (reason exception
#     was raised)
#   - CAF infrastructure error handling (e.g. job submission)
# * Cleanup the implementation of the DE strategies! ('rand-trig',1) is
#   an ugly implementation.... and not really useful as a tuple.
#
# SADE NOTES:
# * Store of ctrlPop members in same pickle file as main DE pop
#   members?  If so, cPickle.dump(f) will alternate between
#   main DE pop and ctrlPop members
# -> YES.
# * Need to check for valid ctrlPop members before use, or
#   is validity ensured through selection?
# * Need to avoid zero division errors in gainFcn, densityFcn
#
import sys, os, shutil, popen2, cPickle, shelve, math, time, gzip
import random, copy, re, urllib, xml.parsers.expat
from optparse import OptionParser

__version__ = '2.4.4'

farmdict = {}
curFarm = ''

def Pause(lockFName='lockme.loc'):
    """ Pause the algorithm by creating a lock file.
    """
    loc = open(lockFName, 'w')
    loc.close()



def SendMail(address, subject, message):
    """ Send mail to 'address'.
    """
    child = popen2.Popen4('mail -s \"%s\" %s' % (subject, address))
    child.tochild.write(message)
    child.tochild.close()
    retVal = child.wait()
    child.fromchild.close()
    return retVal


def WriteLog(fname, message):
    """ Write 'message' to file 'fname'.
    """
    try:            
        log = open(fname, 'a')
        log.write(message)
        log.close()
    except:
        print message


def XmlStartElement(attrName, attrs):
    """ Start element handler for the Expat XML parser.  Used in parsing
    the output of the XML-based CAF best submit CGI script.  We explicitly
    exclude the RedHat 7.3 account group.
    """
    global farmdict
    global curFarm
    if attrName == 'farm':
        curFarm = attrs['name']
        # clear this farm's list of account groups
        farmdict[curFarm] = []
    elif (attrName == 'acctgroup' and
          attrs['name'] != 'group_rh73'):
        farmdict[curFarm].append(attrs)



def StdDevScaleFac(F, cr, np):
    """ Calculate the theoretical scale factor which acts upon the
    mean standard deviation of trial solution elements in the current
    generation to give the mean std dev. of trial solution elements
    for the next generation (this std. dev. is an average over the
    population and the trial vector elements).

    See Zaharie, D., 'Critical Values for the Control Parameters of
    Differential Evolution Algorithms.'
    """
    return math.sqrt(2. * F**2 * cr + (cr**2 - 2. * cr) / np + 1.0)



def OptMutationScaleFac(stdDevScale, cr, np):
    """ Calculate the theoretical mutation scale factor ('F') for the
    DE rand/1/bin algorithm, given a desired std. dev. scale factor.
    
    See Zaharie, D., 'Critical Values for the Control Parameters of
    Differential Evolution Algorithms.'
    """
    return math.sqrt(0.5*((stdDevScale**2 - 1.) / cr + (2. - cr) / np ))



def OptCrossProb(stdDevScale, F, np):
    """ Calculate the theoretically optimal crossover probability
    ('Cr') for the DE rand/1/bin algorithm, given a desired
    std. dev. scale factor.
    
    See Zaharie, D., 'Critical Values for the Control Parameters of
    Differential Evolution Algorithms.'
    """
    try:
        cr = -(np*F**2 - 1) + math.sqrt((np*F**2 - 1)**2 - np*(1 - stdDevScale**2))
    except ValueError, data:
        print data
        cr = 0.2
        print 'Setting cross prob to the eminently reasonable value of %f' % cr
    return cr



class Member:
    """Population member; holds one trial vector + (tuple of) fcn eval(s)"""
    _LARGEVAL = 1.0E20
    _defLoBnd = 0.0
    _defHiBnd = 1.0E3
    _constrEvalMax = 500
    _constrCallMax = 10

    def __init__(self, dim=1, fdim=1, gen=-1, ipop=-1):
        self.nDim = dim
        self.fDim = fdim
        self.generation = gen     # Negative gen --> trial vect
        self.popIndex = ipop      # index in the population (< 0 --> trial)
        self.x = [0.0] * self.nDim
        self.y = [Member._LARGEVAL]*self.fDim
        self.parmBounds = []
        self.deltaConstr = 0.0
        self.isValid = True
        for i in range(self.nDim):
            self.parmBounds.append((Member._defLoBnd, Member._defHiBnd))


    def __repr__(self):
        sep = '\n'
        a = '  Generation: %s \t Population Index: %s' \
            % (self.generation,self.popIndex)
        b = '  Cost: %s \t\t Vector: %s' \
            % (self.y, self.x)
        return sep.join((a,b))


    def isTrial(self):
        return (self.generation < 0 or self.popIndex < 0)


    def getKey(self):
        # generate a string object from self.generation and
        # self.popIndex which can be used as a
        # dictionary key.
        return self.popIndex


    def __getitem__(self, key):
        return self.x[key]


    def __setitem__(self, key, value):
        if type(value) != type(1.0) and type(value) != type(1):
            raise TypeError, 'Only floats and ints can be stored in Member.x[].'
        self.x[key] = value


    def __len__(self):
        return self.x.__len__()


    def __cmp__(self, other):
        # Overload the comparison function -- useful for sorting
        # Member objects
        # May want to define the rich comparison operators
        # (e.g. __lt__() for '<') instead of __cmp__().
        if not isinstance(other, Member):
            raise TypeError, ('Comparison with type %s not supported!' 
                  % (other.__class__))
        if len(self.y) > 1:
            raise ValueError, 'Pareto DE not yet implemented.'
        else:
            if self.y[0] > other.y[0]:
                return 1
            elif self.y[0] == other.y[0]:
                return 0
            else:
                return -1


    def makeTrial(self):
        # Make a copy of self, but tag it as a trial;
        # return the copy
        trialMember = copy.deepcopy(self)
        trialMember.generation = -1
        trialMember.popIndex = -1
        trialMember.y = [Member._LARGEVAL]*self.fDim
        return trialMember

    
    def __add__(self, other):
        """ Overload the '+' operator; add the parameter values for
        each index rather than concatenating the underlying list objects. """
##         if not isinstance(other, Member):
        if not isinstance(other, self.__class__):
            raise TypeError, (
                  'Right-hand operand (type %s) must be a %s object!' 
                  % (other.__class__, self.__class__))
        if other.nDim != self.nDim:
            raise ValueError, (
                  'Member dimension mismatch! left.nDim = %s, right.nDim = %s' 
                  % (trialMember.nDim, other.nDim))
        trialMember = self.makeTrial()
        trialMember.x = map((lambda a,b: a+b), trialMember.x, other.x)
        return trialMember


    def __iadd__(self, other):
        """ Overload the '+=' operator;   """
##         if not isinstance(other, Member):
        if not isinstance(other, self.__class__):
            raise TypeError, (
                  'Right-hand operand (type %s) must be a %s object!' 
                  % (other.__class__, self.__class__))
        if other.nDim != self.nDim:
            raise ValueError, (
                  'Member dimension mismatch! left.nDim = %s, right.nDim = %s' 
                  % (self.nDim, other.nDim))
        self.x = map((lambda a,b: a+b), self.x, other.x)
        return self


    def __sub__(self, other):
        """ Overload the '-' operator; subtract the parameter values for
        each index. """
##         if not isinstance(other, Member):
        if not isinstance(other, self.__class__):
            raise TypeError, (
                  'Right-hand operand (type %s) must be a %s object!' 
                  % (other.__class__, self.__class__))
        if other.nDim != self.nDim:
            raise ValueError, (
                  'Member dimension mismatch! left.nDim = %s, right.nDim = %s' 
                  % (trialMember.nDim, other.nDim))
        trialMember = self.makeTrial()
        trialMember.x = map((lambda a,b: a-b), trialMember.x, other.x)
        return trialMember

    
    def __isub__(self,other):
        """ Overload the '-=' operator; change self.x in place. """
##         if not isinstance(other, Member):
        if not isinstance(other, self.__class__):
            raise TypeError, (
                  'Right-hand operand (type %s) must be a %s object!' 
                  % (other.__class__, self.__class__))
        if other.nDim != self.nDim:
            raise ValueError, (
                  'Member dimension mismatch! left.nDim = %s, right.nDim = %s' 
                  % (self.nDim, other.nDim))
        self.x = map((lambda a,b: a-b), self.x, other.x)
        return self

    
    def __mul__(self, scale):
        """ Overload the '*' operator; scales each
        element of self.x by factor 'scale'. """
        if type(scale) != type(1.0) and type(scale) != type(1):
            raise TypeError, (
                  'Right-hand operand (type %s) must be an int or float!' 
                  % type(scale))
        trialMember = self.makeTrial()
        trialMember.x = map((lambda a,b=scale: a*b), trialMember.x)
        return trialMember


    def __imul__(self, scale):
        """ Overload the '*=' operator; change self.x in place. """
        if type(scale) != type(1.0) and type(scale) != type(1):
            raise TypeError, (
                  'Right-hand operand (type %s) must be an int or float!' 
                  % type(scale))
        self.x = map((lambda a,b=scale: a*b), self.x)
        return self

    
    def __rmul__(self, scale):
        """ Overload the '*' operator for the case when 'scale'
        is the left-hand operand. """
        if type(scale) != type(1.0) and type(scale) != type(1):
            raise TypeError, (
                  'Left-hand operand (type %s) must be an int or float!' 
                  % type(scale))
        trialMember = self.makeTrial()
        trialMember.x = map((lambda a,b=scale: a*b), trialMember.x)
        return trialMember

         
    def __div__(self, scale):
        """ Overload the '/' operator. """
        if type(scale) != type(1.0) and type(scale) != type(1):
            raise TypeError, (
                  'Right-hand operand (type %s) must be an int or float!' 
                  % type(scale))
        if scale == 0:
            raise ZeroDivisionError
        trialMember = self.makeTrial()
        trialMember.x = map((lambda a,b=scale: a/b), trialMember.x)
        return trialMember

    
    def __idiv__(self, scale):
        """ Overload the '/=' operator; change self.x in place """
        if type(scale) != type(1.0) and type(scale) != type(1):
            raise TypeError, (
                  'Right-hand operand (type %s) must be an int or float!' 
                  % type(scale))
        if scale == 0:
            raise ZeroDivisionError 
        self.x = map((lambda a,b=scale: a/b), self.x)
        return self


    def setParmBounds(self, i, loBnd=0.0, hiBnd=1.0E2):
        """ Set the lower and upper bounds for the ith parameter in this
        Member object.
        """
        if type(i) != type(1) or (type(loBnd) != type(1)
                                  and type(loBnd) != type(1.0)) or \
                                  (type(hiBnd) != type(1) and
                                   type(hiBnd) != type(1.0)):
            raise TypeError        
        if i < 0 or i >= self.nDim:
            raise IndexError, 'Invalid index, %s.' % i
        self.parmBounds[i] = (loBnd, hiBnd)


    def crossOver(self, crossProb, mutant):
        """ Perform the crossover operation with 'self' as the parent
        Member object using the binomial distribution; return a new Member object.
        """
        if type(crossProb) != type(1.0):
            raise TypeError
        if crossProb < 0.0 or crossProb > 1.0:
            raise ValueError, 'crossProb should be in [0.0,1.0).'
##         if not isinstance(mutant, Member):
        if not isinstance(mutant, self.__class__):
            raise TypeError, 'crossOver is only defined for two Member objects.'
        trialMember = self.makeTrial()
        trialMember.popIndex = self.popIndex
        # Randomly select one index for crossover
        irand = random.randrange(0, self.nDim)
        for i in range(self.nDim):
            if (i == irand or
                random.random() < crossProb):                
                trialMember[i] = mutant[i]
        return trialMember
    

    def repairHardConstr(self, hardConstraints, callCtr=0):
        """ Repair rule for hard constraint violation.
        """
        evalCtr = 0
        # hardConstraints should be a sequence of strings
        # each string should be of the form 'x[0] - 10.0'
        x = copy.deepcopy(self.x)
        # Create a compiled RE object
        patt = re.compile("\[([0-9]+)\]")
        for constr in hardConstraints:
            # re.compile('pattern').findall('string') returns a list of
            # matching substrings
            delta = 0.0
            j = 0
            indexlist = []
            strlist = patt.findall(constr)            
            if len(strlist) == 0:
                # Malformed constraint! There are no substrings like 'x[0]'.
                # Skip this constraint.
                continue
            for i in strlist:
                # Convert the list indices found to integers
                indexlist.append(int(i))
            # Set j to be the first list index found in the constraint string
            i = 0
            j = indexlist[i]
            try:
                delta = eval(constr)
            except(ZeroDivisionError, ValueError, AssertionError,
                   FloatingPointError, OverflowError):
                delta += Member._LARGEVAL
            except (IndexError, KeyError, NameError, TypeError,
                    SyntaxError, AttributeError):
                # Malformed constraint; skip it. (KeyError will probably
                # only be raised if the underlying parameter container
                # for the Member class is changed to a dictionary)
                continue
            except:
                continue
            while delta > 0.0:
##                 if evalCtr < Member._constrEvalMax:
                if evalCtr < self._constrEvalMax:
                    evalCtr += 1
                else:
                    break
                try:
                    # Throw a random number within the bounds of
                    # parameter j.
                    x[j] = random.uniform(self.parmBounds[j][0],
                                          self.parmBounds[j][1])
                    # evaluate the string constr as a Python expression
                    delta = eval(constr)
                except (ZeroDivisionError, ValueError, AssertionError,
                        FloatingPointError, OverflowError):
                    delta += Member._LARGEVAL
                    # this should not be an infinite loop, given
                    # that we've already successfully eval'd the
                    # constraint once.
                    continue 
                except (IndexError, KeyError):
                    # Fixup j and continue?
                    break
                # Cycle through the list indices found in this
                # constraint -- attempt to systematically vary
                # parameters when constraints involve more than
                # one parameter.  For example, if the constraint
                # string is 'x[0] + 2.0 * x[1]', and it evaluates
                # to a positive definite value (i.e. constraint is
                # violated), first throw a random number for x[0]
                # and test the constraint again; if it is still violated,
                # throw a random number for x[1] and re-eval.
                ######
                # There's probably a better way to code this, e.g. using
                # iterators, but this should work for now...
                if i in range(len(indexlist)-1):
                    i += 1
                    # max value of i is (len(indexlist) - 1)
                else:
                    i = 0
                j = indexlist[i]
        if evalCtr > 0:
            # At least one constraint was violated at least once
##             if callCtr < Member._constrCallMax:
##                 if evalCtr < Member._constrEvalMax:
            if callCtr < self._constrCallMax:
                if evalCtr < self._constrEvalMax:
                    # update the parm vector
                    self.x = x
                    # Call this function recursively to verify that
                    # the new parm vector satisfies the constraints (this
                    # is necessary for a set of constraints which is 
                    # interwoven -- e.g. x[0] occurs in more than one
                    # constraint.
                    self.repairHardConstr(hardConstraints, callCtr+1)
                    return
                else:                    
                    # Having trouble getting a viable parm vector,
                    # so try starting over.
                    self.repairHardConstr(hardConstraints, callCtr+1)
                    return
            else:
                # We've reached the recursive call limit; give up on
                # this Member object and ensure that it will not be
                # passed to the objective function.
                self.isValid = False
                return
        else:
            # Only evaluated each constraint once (exactly), and all constraints
            # were satisfied.  We have a viable parm vector, so return.
            self.isValid = True
            return
            

    def repairSoftConstr(self, softConstraints):
        """ Repair rule for soft constraint violation.
        """
        if self.isValid:
            # softConstraints should be a sequence of strings,
            x = copy.deepcopy(self.x)
            sumDelta = 0.0
            # sum up the total magnitude of the constraint violation
            for constr in softConstraints:
                delta = 0.0
                try:
                    delta = eval(constr)
                except (ZeroDivisionError, ValueError, AssertionError,
                        FloatingPointError, OverflowError):
                    sumDelta += Member._LARGEVAL
                    continue
                except (IndexError, KeyError, NameError, TypeError,
                        SyntaxError, AttributeError):
                    # This constraint seems to be malformed;
                    # ignore it.
                    continue
                if delta > 0.0:
                    sumDelta += delta
            self.deltaConstr = sumDelta


    def dist(self, other, xmin, xmax):
        """ Calculate the normalized Euclidean distance between
        vectors self.x and other.x.  Note: normalization is over the
        parameter x_min[i] and x_max[j] as derived from the current
        population (passed to this method as args).
        """
        # SADE
        # make sure to check for constraint violation before using
        # this distance in an objective fcn.
        # Should have no ZeroDivision problems, as the parmBounds
        # are not user-settable.
        if not isinstance(other, Member):
            raise TypeError, 'Member.dist() only valid btwn two Member objects.'
        d = 0.0
        for i in xrange(self.nDim):
##             d += ((self.x[i] - other.x[i]) / (self.parmBounds[i][1] -
##                                               self.parmBounds[i][0]))**2
            delta = xmax[i] - xmin[i] + 1.0E-6
            d += ((self.x[i] - other.x[i]) / delta)**2
        return math.sqrt(d)


#--------------------------------------------------------------------------


class SAMember(Member):
    """ A subclass of Member which implements a history of function values,
    with the average over the history taken as the reported function value.
    """
    _constrEvalMax = 500
    _constrCallMax = 10
    _TINY = 1.0E-15

    def __init__(self, dim=1, fdim=1, gen=-1, ipop=-1, histlen=10):
        self.nDim = dim
        self.fDim = fdim
        self.histLen = histlen
        self.generation = gen     # Negative gen --> trial vect
        self.popIndex = ipop      # index in the population (< 0 --> trial)
        self.x = [0.0] * self.nDim        
##         self.yHist = [[Member._LARGEVAL]*self.fDim for i in range(self.histLen)]
        self.yHist = [[0.0]*self.fDim for i in range(self.histLen)]
        self.y = [Member._LARGEVAL]*self.fDim
        self.parmBounds = []
        self.deltaConstr = 0.0
        self.isValid = True
        for i in range(self.nDim):
            self.parmBounds.append((Member._defLoBnd, Member._defHiBnd))


    def dominates(self, other):
        """ Return values:
        -1 if self dominates other
        1 if other dominates self
        0 if both self and other are non-dominated
        """
        flag1 = False
        flag2 = False
        for i in range(self.fDim):
            if self.y[i] < other.y[i]:
                flag1 = True
            elif self.y[i] > other.y[i]:
                flag2 = True
        if flag1 and not flag2:
            return -1
        elif not flag1 and flag2:
            return 1
        else:
            return 0
            

    def __cmp__(self, other):
        # Overload the comparison function -- useful for sorting
        # Member objects
        # May want to define the rich comparison operators
        # (e.g. __lt__() for '<') instead of __cmp__().
        if not isinstance(other, SAMember):
            raise TypeError, ('Comparison with type %s not supported!' 
                  % (other.__class__))
        # SADE
        return self.dominates(other)


    def isUnique(self, other):
        """ Determine whether this is a unique individual based upon
        decision vector.
        """
        # calc the infinity-norm of the vector
        delta = self - other
        for i in xrange(len(delta)):
            delta[i] = abs(delta[i])
        norm = max(delta)
        if norm > self._TINY:
            return True
        else:
            return False
        
#--------------------------------------------------------------------------

class Population:
    """An indexed collection of Member objects"""
    _LARGEVAL = 1.0E20
    _TINY = 1.0E-12
    _defInitRange = (-100.0, 100.0)
    _defGaussParms = (0., 50.)  # (mu, sigma)
#    _cmpMetrics = ('stdDev', 'range', 'RMS', 'fractol', 'chisq')
    _strategies = (('rand',1), ('best',1), ('best',2), ('rand-sort',1),
                   ('rand-trig',1), ('best-trig',1))

    def __init__(self, dim=1, fdim=1, popSize=4, gen=-1, prob=0.1, f=0.5,
                 mt=0.05, initGauss=False):
        self.nDim = dim
        self.fDim = fdim
        self.generation = gen  ## May want to take 'gen' out of constructor arglist
        self.np = popSize
        self.crossProb = prob
        self.F = f
        self.mt = mt       ## only used for rand-trig, best-trig strategies
        self.ibest = ''
        self.iworst = ''
        #self.cmpMetric = Population._cmpMetrics[0]
        self.strategy = Population._strategies[0]
        self.hardConstraints = ()
        self.softConstraints = ()
        self.initRange = [(Population._defInitRange[0],
                           Population._defInitRange[1])]*self.nDim
        self.initGauss = initGauss
        self.initGaussParms = {}
        self.memberColl = {}
        for i in range(self.np):
            tmp = Member(self.nDim, self.fDim, self.generation, i)
            tmpkey = tmp.getKey()
            for j in range(tmp.nDim):
                #tmp[j] = tmp.getRandom(j, lo, hi)
                tmp[j] = 0.0
            self.memberColl[tmpkey] = tmp
            if i == 0:
                self.ibest = tmpkey
                self.iworst = tmpkey


    def __getitem__(self, key):
        # SHOULD RAISE TypeError IF key IS AN INVALID TYPE.
        # SHOULD RAISE KeyError IF key IS NOT IN THIS COLLECTION.
        # No. Let the underlying container class do it.
        return self.memberColl[key]


    def __setitem__(self, key, value):
        # SHOULD *NOT* RAISE KeyError IF key IS NOT IN THIS COLLECTION.
##         if not isinstance(value, Member):
##         if not isinstance(value, self.memberColl[0].__class__):
##             raise TypeError, (
##                   'Objects of type %s cannot be stored by Population instances.' 
##                   % (value.__class__))
        if value.isTrial():
            msg = 'Please set the generation and popIndex of value;'
            msg += (' found (gen, popIdx) = (%s, %s).'
                    % (value.generation, value.popIndex))
            raise ValueError, msg
        self.memberColl[key] = value


    def __contains__(self, key):
        return self.memberColl.has_key(key)


    def __len__(self):
        return len(self.memberColl)


    def keys(self):
        return self.memberColl.keys()


    def setGeneration(self, gen):
        """ Set the Population generation, pass this on to all the Member objects.
        """
        self.generation = gen
        for key in self.keys():
            self[key].generation = gen
        

    def setPopParmBounds(self):
        """ Set the low and high bounds for each parameter of each Member
        object in this Population using the available hard constraints.
        Should be called before Population::initialize().
        """
        # Form a compiled regular expression to capture the parm indices
        # in each constraint (e.g. capture '0' in 'x[0] - 1')
        idx = '\[([0-9]+)\]'
        indexpatt = re.compile(idx)
        #
        pm = '([-+])?'
        oppMult = '(?:\s*\*\s*)?'
        oppAdd = '\s*'+pm+'\s*'
        absnum = '((?:\d+(?:\.\d*)?|\d*\.\d+)(?:[eE][-+]?\d+)?)?'
        bndpatt = re.compile(pm + absnum + oppMult + 'x' + idx + oppAdd + absnum)
        # Loop through hard constraints, look for parameter bounds
        bounds = {}
        for constr in self.hardConstraints:
            j = -1
            strlist = indexpatt.findall(constr)
            if len(strlist) == 1:
                # This constraint only involves one parameter,
                # so it is a boundary candidate.
                j = int(strlist[0])
                if j not in bounds.keys():
                    bounds[j] = [Member._defLoBnd, Member._defHiBnd]
                if j < self.nDim and j >= 0:
                    # Determine the low and high bounds
                    match = bndpatt.search(constr)
                    if (len(match.string[match.start():match.end()]) <
                        len(match.string)):
                        # The length of the matched substring is less than
                        # the full string --> there are characters left
                        # over!
                        continue
                    mgroups = match.groups()
                    if mgroups:
                        # constraint is of the form 's0 * c0 * x[j] + s1 * c1 <= 0'
                        # where the '<= 0' is implicit and c0 > 0, c1 >= 0,
                        # s0, s1 = +- 1.0
                        # --> s0 * x[i] <= - s1 * c1 / c0
                        #
                        # >>> m = bndpatt.search('-55e-2x[0]+3.2')
                        # >>> m.groups()
                        # ('-', '55e-2', '0', '+', '3.2')
                        #
                        c0 = 1.0
                        s1 = 1.0 # sign of c1
                        c1 = 0.0 # c1 >= 0
                        if mgroups[1]:
                            c0 = float(mgroups[1])
                        if mgroups[3] == '-':
                            s1 = -1.0
                        if mgroups[4]:
                            c1 = float(mgroups[4])
                        if mgroups[0] == '+' or mgroups[0] == None:
                            # s0 > 0
                            # we have an upper bound
                            bounds[j][1] = -s1 * c1 / c0
                        elif mgroups[0] == '-':
                            # s0 < 0
                            # we have a lower bound
                            bounds[j][0] = s1 * c1 / c0
        for key in self.keys():
            for j in bounds.keys():
                self[key].setParmBounds(j, bounds[j][0], bounds[j][1])


    def initialize(self):
        """ Initialize the population with random trial vectors;
        should only be done for generation -1!"""
        # Set the parameter upper and lower bounds from the hard constraints
        # given by the user.
        self.setPopParmBounds()
        for key in self.keys():
            # Initialize the Member objects' parameter vectors
            for i in range(self[key].nDim):
                if not self.initGauss:
                    self[key][i] = random.uniform(self.initRange[i][0],
                                                  self.initRange[i][1])
                else:
                    if i in self.initGaussParms.keys():
                        self[key][i] = random.gauss(self.initGaussParms[i][0],
                                                    self.initGaussParms[i][1])
                    else:
                        self[key][i] = random.uniform(self.initRange[i][0],
                                                      self.initRange[i][1])
            # Repair the hard and soft constraints
            self[key].repairHardConstr(self.hardConstraints)
            self[key].repairSoftConstr(self.softConstraints)
            self[key].y = [Member._LARGEVAL]*self.fDim
        self.setGeneration(0)


    def getKey(self):
        """ Get a key for this Population instance for use with the
        shelve module. """
        sep1 = ":"
        sep2 = "::"
        key1 = sep1.join(("gen", `self.generation`))
        key2 = sep1.join(("np", `self.np`))
        key3 = sep1.join(("nDim", `self.nDim`))
        return sep2.join((key1, key2, key3))


    def setStrategy(self, strat):
        if strat not in Population._strategies:
            raise ValueError, ('%s is not an implemented DE strategy.' 
                  % strat)
        self.strategy = strat


    def printAvailStrategies(self):
        print "Available DE strategies:\n%s" % Population._strategies


    def getSortedKeys(self):
        """ Sort the Member objects in this population by increasing cost
        value.
        """
        mbrkeys = self.keys()
        # sort the list of keys in place according to the cost value
        mbrkeys.sort(lambda a,b: Member.__cmp__(self[a], self[b]))
        return mbrkeys


    def getNpTrunc(self, truncFrac=0.15):
        """ Get the truncated number of individuals for use in convergence tests.
        """
        return max(int(math.ceil(truncFrac*self.np)), self.nDim + 1, 3)
    

    def getCovMatrix(self, n=2):
        """ Get the covariance matrix of the parameters for truncated
        population defined by truncFrac. Return value is a tuple:
        (list of mean values, covariance matrix)
        """
        covMat = [[Population._LARGEVAL]*self.nDim]
        for i in range(self.nDim-1):
            covMat.append(copy.deepcopy(covMat[0]))
        if n <= 0:
            return ([Population._LARGEVAL]*self.nDim,
                    covMat)
        sortedKeys = self.getSortedKeys()
        truncKeys = sortedKeys[:n]
        # Calc the mean values
        sum = 0.0
        mu = []
        for iparm in range(self.nDim):
            sum = 0.0
            for key in truncKeys:
                sum += self[key].x[iparm]
            mu.append(sum / float(n))
        # Now calc the covariance matrix
        sum = 0.0
        for iparm in range(self.nDim):
            # Cov matrix should be symmetric, so only calc the
            # upper triangle.
            for jparm in range(iparm, self.nDim):
                sum = 0.0
                for key in truncKeys:
                    sum += ((self[key].x[iparm] - mu[iparm]) *
                            (self[key].x[jparm] - mu[jparm]))
                covMat[iparm][jparm] = sum / float(n)
        # Now set the lower triangle...
        for iparm in range(1, self.nDim):
            for jparm in range(iparm):
                covMat[iparm][jparm] = copy.deepcopy(covMat[jparm][iparm])
        return (mu, covMat)


    def getCovMatRepr(self, covMat):
        """ Print the covariance matrix in a pretty way.
        """
        str = ''
        for i in range(len(covMat)):
            for j in range(len(covMat[i])):
                str += '%.6E   ' % covMat[i][j]
            str += '\n'
        return str
    

    def getStats(self, truncFrac=0.15):
        """ Get the stats for this population.
        """
        # Return a tuple of the form,
        # (bestCost, worstCost, meanCost, stdDev, fracTol, chisq, ndf)
        orderedKeys = self.getSortedKeys()
        # Get the truncated list of keys
        # The number of keys in the truncated list must be at least (nDim+1)
        # so that ndf >= 1 (must have at least one degree of freedom)
        npTrunc = self.getNpTrunc(truncFrac)
        ndf = npTrunc - self.nDim
        truncKeys = []
        truncSuccess = True
        for key in orderedKeys:
            if (self[key].isValid and
                self[key].deltaConstr == 0.0):
                truncKeys.append(key)
        truncKeys = truncKeys[:npTrunc]
        if len(truncKeys) < npTrunc:
            truncSuccess = False
        # Get the best and worst Member objects
        # Ensure that they are viable
        # -- Actually, the viability check for the best Member obejct
        # should not be necessary, since all nonviable objects are
        # automatically assigned a cost value of Population._LARGEVAL
        i = 0
        while (i < len(orderedKeys) and
               (not self[orderedKeys[i]].isValid or
                self[orderedKeys[i]].deltaConstr > 0)):
            i += 1            
        if i < len(orderedKeys):
            self.ibest = orderedKeys[i]
        else:
            # We've traversed the entire list of keys and not found
            # a viable Member object, so choose randomly
            self.ibest = random.choice(orderedKeys)
        i = -1
        while (i >= -len(orderedKeys) and
               (not self[orderedKeys[i]].isValid or
                self[orderedKeys[i]].deltaConstr > 0)):
            i -= 1
        if i >= -len(orderedKeys):
            self.iworst = orderedKeys[i]
        else:
            # we've traversed the entire list of keys and not found
            # a viable Member object
##            self.iworst = self.ibest
            self.iworst = random.choice(orderedKeys)
        #### SADE
        bestCost = self[self.ibest].y[0]
        worstCost = self[self.iworst].y[0]
        if self.ibest == self.iworst:
            # We've got problems -- not enough viable Member objects!
            # Returning Population._LARGEVAL for most stats will ensure that we
            # do not converge early.
            (muParms, covMat) = self.getCovMatrix(npTrunc)
            return (bestCost, worstCost, Population._LARGEVAL,
                    Population._LARGEVAL,
                    Population._LARGEVAL,
                    Population._LARGEVAL,
                    ndf, muParms, covMat)
        # Find the mean cost
        sum = 0.0
        sumsq = 0.0
        sumTrunc = 0.0
        sumsqTrunc = 0.0
        mu = Population._LARGEVAL
        musq = Population._LARGEVAL

        stdDev = Population._LARGEVAL
        varTrunc = Population._TINY
        chisq = Population._LARGEVAL
        fracTol = Population._LARGEVAL
        # Calc the mean and standard deviation together
        n = 0
        for key in self.keys():
            # Only include valid Population Members in the costs;
            # we're only going to use this list of cost values for
            # statistical & convergence purposes, so we don't want to
            # skew the results with non-viable Members.
            if (self[key].isValid and
                self[key].deltaConstr == 0.0):
                sum += self[key].y[0]
                sumsq += self[key].y[0]**2
                n += 1
        if n > 0:
            mu = sum / float(n)        
            musq = sumsq / float(n)
        diff = musq - mu**2
        if diff > 0.0 and n > 1:
            # Calc the standard deviation for the entire population
            stdDev = math.sqrt(n * diff / float(n - 1))
        # Loop through the sorted, truncated list of keys,
        # excluding the best individual
        if truncSuccess:
            #for key in truncKeys[1:]:
            for key in truncKeys:            
                # We've already checked every Member key in truncKeys[] for
                # viability and discarded the unviable ones
                sumTrunc += self[key].y[0]
                sumsqTrunc += self[key].y[0]**2
##            muTrunc = sumTrunc / float(npTrunc - 1)
##            musqTrunc = sumsqTrunc / float(npTrunc - 1)
            muTrunc = sumTrunc / float(npTrunc)
            musqTrunc = sumsqTrunc / float(npTrunc)
            diffTrunc = musqTrunc - muTrunc**2
            if diffTrunc > 0:
##                varTrunc = (npTrunc - 1) * diffTrunc / float(npTrunc - 2)
                varTrunc = npTrunc * diffTrunc / float(npTrunc - 1)
            chisq = 0.0
            for key in truncKeys[1:]:
                chisq += (self[key].y[0] - bestCost)**2
            #chisq /= (varTrunc + Population._TINY)
            #chisq /= (varTrunc)
            if abs(musqTrunc) > Population._TINY:
                chisq /= (musqTrunc)
            else:
                chisq = Population._LARGEVAL
            # 
            ######chisq = (bestCost - muTrunc)**2 / varTrunc
            range = self[truncKeys[-1]].y[0] - bestCost
            fracTol = 2.0 * abs(range) / (abs(self[truncKeys[-1]].y[0]) +
                                          abs(bestCost) + Population._TINY)
        else:
            # Calculate the fractional tolerance
            range = worstCost - bestCost
            fracTol = 2.0 * abs(range) / (abs(worstCost) + abs(bestCost) +
                                          Population._TINY)
        (muParms, covMat) = self.getCovMatrix(npTrunc)
        return (bestCost, worstCost, mu, stdDev, fracTol, chisq, ndf,
                muParms, covMat)
        

    def getRndMembers(self, nMembers, *targets):
        """ Randomly select 'nMembers' from the Population instance.
        They must all be different from each other and from 'targets'.
        Returns a tuple with the selected Member objects' keys. """
        # This is much like the random.sample(population,k) utility,
        # but we need to exclude 'targets' as well.
        if nMembers >= self.np:
            raise ValueError, 'Requested more random members than are in the population!'
        rndMemberKeys = []
        keys = self.keys()
        for target in targets:
            rndMemberKeys.append(target)
        for i in range(nMembers):
            tmp = random.choice(keys)
            while tmp in rndMemberKeys:
                tmp = random.choice(keys)
            rndMemberKeys.append(tmp)
        for target in targets:
            rndMemberKeys.remove(target)
        return tuple(rndMemberKeys)


    def getMutant(self, parentKey):
        """ Generate a mutant Member object according to the current
        Differential Evolution strategy. Accepts a parent key as input
        (not a Member object), and returns a new Member object. """
        if parentKey not in self.keys():
            raise KeyError, ('Key %s is not in this Population!' 
                  % parentKey)
        mutant = self[parentKey].makeTrial()
        rndKeys = ()
        if self.strategy == ('best', 1):
            rndKeys = self.getRndMembers(2, parentKey, self.ibest)
            mutant = self[self.ibest] + self.F * (self[rndKeys[0]] -
                                                  self[rndKeys[1]])
        elif self.strategy == ('best', 2):
            rndKeys = self.getRndMembers(4, parentKey, self.ibest)
            mutant = self[self.ibest] + self.F * (self[rndKeys[0]] +
                                                  self[rndKeys[1]] -
                                                  self[rndKeys[2]] -
                                                  self[rndKeys[3]])
        elif self.strategy == ('best-trig',1):
            #### SADE
            rndKeys = self.getRndMembers(2, parentKey, self.ibest)
            pprime = (abs(self[rndKeys[0]].y[0]) + abs(self[rndKeys[1]].y[0]) +
                      abs(self[self.ibest].y[0]))
            if random.random() <= self.mt and pprime > Population._TINY:
                # 'mutant' is biased toward the region of lower function
                # value in the available parameter space (defined by
                # the 'rndKeys' list)
                p = [abs(self[self.ibest].y[0]) / pprime,
                     abs(self[rndKeys[0]].y[0]) / pprime,
                     abs(self[rndKeys[1]].y[0]) / pprime]
                mutant = ((self[self.ibest] + self[rndKeys[0]] +
                           self[rndKeys[1]]) / 3.0 +
                          (p[1] - p[0]) * (self[self.ibest] - self[rndKeys[0]]) +
                          (p[2] - p[1]) * (self[rndKeys[0]] - self[rndKeys[1]]) +
                          (p[0] - p[2]) * (self[rndKeys[1]] - self[self.ibest]))
            else:                
                mutant = self[self.ibest] + self.F * (self[rndKeys[0]] -
                                                      self[rndKeys[1]])
##        elif self.strategy == ('rand-to-best', 1):
##            rndKeys = self.getRndMembers(3, parentKey, self.ibest)
##            mutant = (self[rndKeys[0]] + self.lmbda * (self[self.ibest] -
##                                                       self[rndKeys[0]]) +
##                      self.F * (self[rndKeys[1]] - self[rndKeys[2]]))
        elif self.strategy == ('rand-sort', 1):
            rndKeys = list(self.getRndMembers(3, parentKey))
            rndKeys.sort(lambda a,b: Member.__cmp__(self[a], self[b]))
            shuffledKeys = rndKeys[1:]
            random.shuffle(shuffledKeys)
            mutant = self[rndKeys[0]] + self.F * (self[shuffledKeys[0]] -
                                                  self[shuffledKeys[1]])
        elif self.strategy == ('rand', 1):
            # assume self.strategy == ('rand', 1):            
            rndKeys = self.getRndMembers(3, parentKey)
            mutant = self[rndKeys[0]] + self.F * (self[rndKeys[1]] -
                                                  self[rndKeys[2]])
        elif self.strategy == ('rand-trig',1):
            rndKeys = self.getRndMembers(3, parentKey)
            pprime = (abs(self[rndKeys[0]].y[0]) + abs(self[rndKeys[1]].y[0]) +
                      abs(self[rndKeys[2]].y[0]))
            if random.random() <= self.mt and pprime > Population._TINY:
                # 'mutant' is biased toward the region of lower function
                # value in the available parameter space (defined by
                # the 'rndKeys' list)
                p = [abs(self[rndKeys[0]].y[0]) / pprime,
                     abs(self[rndKeys[1]].y[0]) / pprime,
                     abs(self[rndKeys[2]].y[0]) / pprime]
                mutant = ((self[rndKeys[0]] + self[rndKeys[1]] +
                           self[rndKeys[2]]) / 3.0 +
                          (p[1] - p[0]) * (self[rndKeys[0]] - self[rndKeys[1]]) +
                          (p[2] - p[1]) * (self[rndKeys[1]] - self[rndKeys[2]]) +
                          (p[0] - p[2]) * (self[rndKeys[2]] - self[rndKeys[0]]))
            else:
                # use standard DE/rand/1 strategy
                rndKeys = self.getRndMembers(3, parentKey)
                mutant = self[rndKeys[0]] + self.F * (self[rndKeys[1]] -
                                                      self[rndKeys[2]])
        else:
            raise ValueError
        return mutant
            

    def __repr__(self):
        stats = self.getStats()
        str = """Object type: class Population
        Generation: %s \t Population Size: %s \t Depth: %s """ % (
            self.generation, self.np, self.nDim)
        str += """\nStats:
        Best cost: %s \t Worst cost: %s
        Mean cost: %s \t Standard Deviation: %s
        Fractional Tolerance: %s \t Chi-Square Tolerance: %s \t NDF: %s
        Mean parameter values: %s
        Covariance Matrix: %s
        """ % stats
        return str


    def isComplete(self):
        """ Test to see whether all Member objects in this population
        have the same generation; returns the number of Member objects
        which have gen index less than the population gen index."""
        # *** should change the name, since 'is' in the name implies
        # *** that the function will return a Boolean!
        if len(self) < self.np:
            print """Ack! Population should have %s Member objects, but
            only found %s!""" % (self.np, len(self))
        count = 0        
        for key in self.keys():
            if self[key].generation < self.generation:
                count += 1
            elif self[key].generation > self.generation:
                # found a Member object with a larger generation index;
                # can either raise an error or take on the new gen index.
                # For now, print a message; if we do see a higher gen index
                # we should call isComplete() recursively, since count needs
                # to be updated.
                print "Ack! Member object with key %s has a higher gen index!" \
                      % key
        if count == 0:
            orderedKeys = self.getSortedKeys()
            self.ibest = orderedKeys[0]
            self.iworst = orderedKeys[-1]
        return count


    def getTrialMember(self, parentKey):
        """ Generate a new trial vector to possibly replace the ith member
        from members of the current population."""
##        if not isinstance(parent, Member):
##            raise TypeError, ('Member.getTrialMember(parent): parent must be a %s instance.' 
##                  % self.__class__)        
        # trial members get a negative generation number
        mutant = self.getMutant(parentKey)
        # Perform the crossover operation -- trialMember will get the
        # popIndex of the parent
        trialMember = self[parentKey].crossOver(self.crossProb, mutant)
        # Enforce Hard constraints
        trialMember.repairHardConstr(self.hardConstraints)
        # Check Soft constraints
        trialMember.repairSoftConstr(self.softConstraints)
        return trialMember


    def saveMe(self, filename='popDB.dbm'):
        """ Append the current Population instance to a DBM file with the
        Python shelve module (filename can refer to a new or existing file). """
        popDB = shelve.open(filename)
        popDB[self.getKey()] = self
        popDB.close()


#--------------------------------------------------------------------------

class SAPopulation(Population):
    """ Subclass of Population class for use with self-adaptive DE.

    x[0] -> F
    x[1] -> Cr

    y[0] -> density objective
    y[1] -> gain objective
    """
    
    _defInitRange = (0.01,0.99)
    
    def __init__(self, dim=1, fdim=1, popSize=4, gen=-1, prob=0.1,
                 f=0.5, mt=0.05, histlen=5, lam=0.2, initGauss=False):
        self.nDim = dim
        self.fDim = fdim
        self.histLen = histlen
        self.lam = lam
        self.fi = 0
        self.generation = gen  ## May want to take 'gen' out of constructor arglist
        self.np = popSize
        self.crossProb = prob
        self.F = f
        self.mt = mt       ## only used for rand-trig, best-trig strategies
        self.initGauss = initGauss
        self.initGaussParms = {}
        self.ibest = ''
        self.iworst = ''
        #self.cmpMetric = Population._cmpMetrics[0]
        self.strategy = ('rand-trig',1)
        self.hardConstraints = ()
        self.softConstraints = ()
        self.initRange = [(SAPopulation._defInitRange[0],
                           SAPopulation._defInitRange[1])]*self.nDim
        self.memberColl = {}        
        for i in range(self.np):
            tmp = SAMember(self.nDim, self.fDim, self.generation, i, self.histLen)
            tmpkey = tmp.getKey()
            for j in range(tmp.nDim):
                #tmp[j] = tmp.getRandom(j, lo, hi)
                tmp[j] = 0.0
            self.memberColl[tmpkey] = tmp
            if i == 0:
                self.ibest = tmpkey
                self.iworst = tmpkey
        

    def getMutant(self, parentKey):
        """ Generate a mutant Member object according to the modified
        DE/rand/1 scheme.  Accepts a parent key as input and returns a new
        Member object.
        """
        if parentKey not in self.keys():
            raise KeyError, ('Key %s is not in this Population!' 
                  % parentKey)
        rndKeys = ()
        mutant = self[parentKey].makeTrial()                
        mu = 0.5
        #sigma = 0.1
        #sigma = 0.3
        sigma = 0.4
##         scale = random.gauss(mu,sigma)
        # Breit-Wigner / Cauchy pdf
        # scale = mu + sigma * math.tan(math.pi * (random.random() - 0.5))
        scale = random.uniform(0.0,1.0)
        if self.strategy == ('rand',1):
            rndKeys = self.getRndMembers(3, parentKey)
            mutant = self[rndKeys[0]] + scale * (self[rndKeys[1]] -
                                               self[rndKeys[2]])
        elif self.strategy == ('rand-trig',1):
            rndKeys = self.getRndMembers(3, parentKey)
            if random.random() <= self.mt:
                # 'mutant' is biased toward the region of lower function
                # value in the available parameter space (defined by
                # the 'rndKeys' list)
                pprime = (abs(self[rndKeys[0]].y[self.fi]) +
                          abs(self[rndKeys[1]].y[self.fi]) +
                          abs(self[rndKeys[2]].y[self.fi]))
                p = [abs(self[rndKeys[0]].y[self.fi]) / pprime,
                     abs(self[rndKeys[1]].y[self.fi]) / pprime,
                     abs(self[rndKeys[2]].y[self.fi]) / pprime]
                mutant = ((self[rndKeys[0]] + self[rndKeys[1]] +
                           self[rndKeys[2]]) / 3.0 +
                          (p[1] - p[0]) * (self[rndKeys[0]] - self[rndKeys[1]]) +
                          (p[2] - p[1]) * (self[rndKeys[1]] - self[rndKeys[2]]) +
                          (p[0] - p[2]) * (self[rndKeys[2]] - self[rndKeys[0]]))
            else:
                # use standard DE/rand/1 strategy
                rndKeys = self.getRndMembers(3, parentKey)
                mutant = self[rndKeys[0]] + scale * (self[rndKeys[1]] -
                                                     self[rndKeys[2]])
        return mutant


    def front(self, SubPopKeys):
        """ Recursive function to identify the non-dominated
        individuals in the population.
        """
        if len(SubPopKeys) == 1:
            return SubPopKeys
        else:
            halfkey = int(math.floor(len(SubPopKeys)/2.))
            # "Top half" of sub-population
            TKeys = self.front(SubPopKeys[:halfkey])
            # "Bottom half" of sub-population
            BKeys = self.front(SubPopKeys[halfkey:])
            # "Merged" sub-population
            MKeys = copy.deepcopy(TKeys)
            # Create the merged set: M := T Union {B_i},
            # where {B_i} is the set of Members in B which
            # are non-dominated w.r.t. all Members in T
            nondomflag = True
            for bkey in BKeys:
                nondomflag = True
                for tkey in TKeys:
                    if self[tkey] < self[bkey]:
                        # tkey dominates bkey, so bkey does not
                        # get added to the merged population.
                        nondomflag = False
                        break
                if nondomflag:
                    MKeys.append(bkey)
            return MKeys


    def nonDomSet(self):
        """ Kung et al.'s Efficient Method for sorting a popultion to
        identify the non-dominated individuals (Ref: Kung et
        al. 1975).
        """
        # First, sort the population by increasing value of the first
        # objective function value
        #
        # There's got to be a better way to do this than making a copy of
        # the underlying dictionary...
        PKeyVals = sorted(self.memberColl.items(),
                          lambda a,b: cmp(a[1].y[0],
                                          b[1].y[0]))
        PKeys = []
        for kv in PKeyVals:
            PKeys.append(kv[0])
        # Now recursively sort the set
        return self.front(PKeys)


##     def update(self):
##         """ Update the population; remove all dominated individuals.
##         """
##         if self.generation % 2 == 0:
##             nonDomKeys = self.nonDomSet()
##             mbrkeys = self.keys()
##             for key in mbrkeys:
##                 if key not in nonDomKeys:
##                     mbrkeys.remove(key)
##             while len(mbrkeys) < 4:
##                 # need to make sure we have enough members to
##                 # perform a mutation for the control parms.
##                 mbrkeys.extend(self.getRndMembers(1,*mbrkeys))
##             for key in self.keys():
##                 if key not in mbrkeys:
##                     del self.memberColl[key]


#     def update(self):
#         """ Update the population; remove all dominated individuals.
#         """
#         if self.generation % 2 == 0:
#             nonDomKeys = self.nonDomSet()
# ##             while len(nonDomKeys) < 4:
#             #minkeys = max(int(math.ceil(0.10*self.np)), 4)
#             minkeys = max(int(math.ceil(0.30*self.np)), 8)
#             maxtrials = len(self.memberColl)
#             ntrials = 0
#             while (len(nonDomKeys) < minkeys and
#                    ntrials < maxtrials):
#                 # need to make sure we have enough members to
#                 # perform a mutation for the control parms.
#                 tmp = self.getRndMembers(1,*nonDomKeys)
#                 uniq = True
#                 for key in nonDomKeys:
#                     if not self.memberColl[tmp[0]].isUnique(self.memberColl[key]):
#                         uniq = False
#                 if uniq:
#                     nonDomKeys.extend(tmp)
#                 ntrials += 1
#             for key in self.keys():
#                 if key not in nonDomKeys:
#                     del self.memberColl[key]


    def update(self):
        """ Update the population; remove all dominated individuals.
        """
        if self.generation % 2 == 0:
            nonDomKeys = self.nonDomSet()
##             while len(nonDomKeys) < 4:
            # minkeys = max(int(math.ceil(0.10*self.np)), 4)
            #minkeys = max(int(math.ceil(0.10*self.np)), 10)
            minkeys = max(int(math.ceil(0.20*self.np)), 8) # min(20,self.np)
            tmpkeys = copy.deepcopy(self.keys())
            random.shuffle(tmpkeys)
            if len(tmpkeys) > minkeys:
                for key in tmpkeys:
                    if key not in nonDomKeys and len(self.keys()) > minkeys:
                        del self.memberColl[key]


    def densityFcn(self, mainPop, child):
        """ The distance objective function for use in self-adapting
        the DE control parms.  The largest hypersphere in which
        'memKey' is the only individual is a measure of the
        crowdedness around individual 'memKey' (NSGA-II's density
        function).  The nearest neighbor will reside on the surface of
        this hypersphere.

        Could use niche count instead (as in NPGA); such kernel-based density
        estimation techniquies have been shown to perform better than
        the 'nearest neighbor' technique implemented here (also in
        NSGA-II).  Reference:
        
        M.Laumanns, E.Zizler, L.Thiele, On the Effects of Archiving,
        Elitism, and Density Based Seelction in Evolutionary
        Multi-objective Optimization.  In Proceedings of EMO 2001,
        pp. 181-196, 2001.
        """
        # NB! This is a minimization objective for SADE
        d = []
        x = [0.0]*mainPop.nDim
        xmin = [0.0]*mainPop.nDim
        xmax = [0.0]*mainPop.nDim
        for i in xrange(mainPop.nDim):
            for key in mainPop.keys():
                x[i] = mainPop[key][i]
            xmin[i] = min(x)
            xmax[i] = max(x)
        for key in mainPop.keys():
            if key != child.getKey():
                d.append(child.dist(mainPop[key], xmin, xmax))
        return -min(d)


    def gainFcn(self, mainPop, parent, child):
        """ The norm'd distance traveled 'downward' in function space
        for use in self-adapting the DE control parms. NB! The
        normalization is w.r.t. the f_min and f_max found in the
        current gen's population.

        INPUT: parent -- parent Member object
                         in the current gen's population
               child -- child Member object which is
                              challenging parent (NB! needs to be
                              a full Member object!)
        """
        # Implementation for arbitrary number of objective functions:
        # NB! Minimization objective! (reflected in sign of 'gain')
        fmin, fmax, gain = 0.0, 0.0, 0.0
        flist = []
        # Collect all the function values for the main population
        # so we can normalize the gain between parent and child.
        for i in xrange(mainPop.fDim):
            flist = []
            flist.append(child.y[i])
            for key in mainPop.keys():
                flist.append(mainPop[key].y[i])
            fmin = min(flist)
            fmax = max(flist)
            delta = fmax - fmin + 1.0E-6
            gain += (child.y[i] - parent.y[i]) / delta
        return gain


    def sortedFront(self, fi=0):
        """ Sort the population according to the given objective
        function index 'fi'.  Returns an ordered list of Member keys.
        """
        sort = sorted(self.memberColl.items(),
                      lambda a,b: cmp(a[1].y[fi],
                                      b[1].y[fi]))
        keys = []
        for item in sort:
            keys.append(item[0])
        return keys


    def getPopFrac(self, zeta=0.1, fi=0):
        """ Utility function; returns a list of keys (from the current
        population) which represent the first 'zeta' fraction of the
        population, sorted according to objective 'fi'.
        """
        sfKeys = self.sortedFront(fi)
        if zeta > 1.0 or zeta < 0.:
            raise ValueError, 'zeta required to be in [0,1].'
        uindex = int(zeta*len(sfKeys))
        if uindex < 1:
            uindex = 1
        return sfKeys[:uindex]
    

    def selCtrlMember(self, zeta=0.1, fi=0):
        """ Select the DE control parameters to use from the Pareto
        front according to parameter 'zeta'.  Returns a SAMember object.
        """
        return copy.deepcopy(self.memberColl[random.choice(self.getPopFrac(zeta,fi))])


    def evalCtrl(self, mainPop, parent, child, ctrlChild):
        """ Evaluate a set of DE control parameters according to the
        parent and child of the main population.

        OUTPUT: y[0] -> density objective
                y[1] -> gain objective
        """
        del ctrlChild.yHist[0:1]
        ctrlChild.yHist.append([self.densityFcn(mainPop, child),
                                self.gainFcn(mainPop, parent, child)])
        ysum = [0.0,0.0]
        n = range(self.histLen)
        w = []
        for i in n:
            w.append(math.e**(-self.lam*i))
            ysum[0] += ctrlChild.yHist[i][0] * w[i]
            ysum[1] += ctrlChild.yHist[i][1] * w[i]
        return ysum


    def __repr__(self):
        sortedKeysA = self.sortedFront(0)
        sortedKeysB = self.sortedFront(1)
        str = 'Most exploratory SAMember:\n'
        str += `self[sortedKeysA[0]]`
        sortedKeys = self.sortedFront(1)
        str += '\nMost exploitative SAMember:\n'
        str += `self[sortedKeysB[0]]`
        str += '\n'
        return str


#--------------------------------------------------------------------------


class DeDriver:
    """ Driver class for the parellelized Differential Evolution algorithm.
    """
    # Dictionary of available farms and their crucial info:
    # entries are like this: (group, available queues)
    # Standard queues and their limits (for now just realtime limit, in seconds)
    _stdQueueData = {}
    _stdQueueData['short'] = 6 * 3600.0
    _stdQueueData['medium'] = 30 * 3600.0
    _stdQueueData['long'] = 72 * 3600.0
    _stdQueues = _stdQueueData.keys()
    # Only dCAFs running Scientific Linux, for now
    # _availFarms dictionary keys are the strings used to specify the
    # farm at submission time; the last element of the value tuple is
    # the farm's 'official' name (really only needed for cnaf!)
    _availFarms = {}
    _availFarms['local'] = ('common', _stdQueues, 'local')
##     _availFarms['ASCAF'] = ('cdf', _stdQueues, 'ascaf')
    _availFarms['BCNCAF'] = ('common', _stdQueues, 'bcncaf')
    _availFarms['RUTCAF'] = ('common', _stdQueues, 'rutcaf')
##     _availFarms['KORCAF'] = ('common', _stdQueues, 'korcaf')
    _availFarms['GroupCAF'] = ('common', _stdQueues, 'caf')
    _availFarms['Fermigrid'] = ('common', _stdQueues, 'fermigrid')
    _availFarms['NAmCAF'] = ('common', _stdQueues, 'namcaf')
##     _availFarms['LCGCAF'] = ('common', _stdQueues, 'lcgcaf')
##     _availFarms['CNAFCAF'] = ('common', _stdQueues, 'cnaf')
##     _availFarms['LyonCAF'] = ('common', ('short', 'medium'), 'lyoncaf')
##     _availFarms['SDSCCAF'] = ('common', _stdQueues, 'sdsccaf')
##     _availFarms['TORCAF'] = ('common', _stdQueues, 'torcaf')
    #### STOP STOP STOP
    #### WHAT ABOUT OTHER GROUPS? IF THE USER SPECIFIES 'MCprod', WE SHOULD
    #### BE ABLE TO SUBMIT JOBS USING 'MCprod' OR 'common'!!
    # Define the available convergence strategies
    _availConvgStrats = ('chisq', 'fractol')
    #_safetyFac = 3600.0
    # Define an additive safety factor, used in determining whether we should
    # worry that the last job submitted has not completed.  Make it rather
    # long to account for busy dCAFs.
    #_safetyFac = 0.75 * _stdQueueData['short']
    _safetyFac = 0.85 * _stdQueueData['short']


    def __init__(self, dim=1, fdim=1, np=4, gen=-1, fname='de_state.dat.gz',
                 email='', costFcn=None, costArgs=None):
        """ Constructor for the DeDriver class.
        """
        # Cost function (python function) for local running:
        self.costFcn = costFcn
        self.costArgs = costArgs
        self.sade = False
        # Create the population(s)
        self.population = Population(dim, fdim, np, gen)
        # Create the control parm population        
##         self.ctrlPop = SAPopulation(2, 2, int(np/2.0), gen)
        self.ctrlPop = SAPopulation(2, 2, np, gen)
        # Set the control parm population's constraints
##         # The radical of D.Zaharie's critical scale factor should be positive:
##         rad = ('2.*x[0]**2*x[1]+(x[1]**2-2.*x[1])/'
##                + `self.ctrlPop.np` + '1.0')
##         # The actual critical variance scale factor:
##         critSF = '-math.sqrt(' + rad + ')+1.0'
        self.ctrlEps = 3E-2
        critSF = '-2.0*%d*x[0]**2+2.0-x[1]+2*%.1E*%d/x[1]' % (self.population.np,
                                                              self.ctrlEps,
                                                              self.population.np)
        self.ctrlPop.hardConstraints = ('-x[0]+0.05','x[0]-0.99',
                                        '-x[1]+0.0', 'x[1]-0.99', critSF)
##         self.ctrlPop.softConstraints = (critSF,)
##         self.ctrlPop.hardConstraints = ('-x[0]+0.05','x[0]-0.99',
##                                         '-x[1]+0.0', 'x[1]-0.5')
        self.ctrlPop.initRange[0] = (0.3,0.99)
        self.ctrlPop.initRange[1] = (0.0,0.3)
        #
        self.stateFName = fname
        # load the population from a file?
        self.initFromFile = False
        self.saveInitPop = False
        self.popFile = ''
        # Differential Evolution Parameters
        self.pareto = False
        self.xi = 0.25 # phase 1 strategy used for first (xi*Gmax) generations
        self.dimThresh = 6
        self.zeta = 0.4
        self.fi = 0
        # Phase 1 (Exploration) parameters
        self.phase1LoDimCr = 0.5
        self.phase1HiDimCr = 0.7
        self.phase1F = 0.8
        self.phase1Strategy = ('rand', 1)             # 'DE/rand/1/bin'
        # Phase 2 (Elitism) parameters
        self.phase2LoDimCr = 0.3
        self.phase2HiDimCr = 0.4
        self.phase2F = 0.6
        self.phase2Strategy = ('best', 2)             # 'DE/best/2/bin'
        # number of previous generations to look at for convergence
        self.m = 5
        # List to hold the statistics of the previous 'm' generations
        # Use integer keys, value should be a tuple:
        # (bestCost, worstCost, meanCost, stdDev, fracTol, chisq, ndf)
        # For multi-objective optimization, each tuple entry should
        # be a tuple whose elements correspond to the various
        # functions.
        self.prevMGenStats = [(Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL)] * self.m
        # Construct a list which will hold the max
        # delta(Stat(i,gen_j), Stat(i,gen_j-1)) for the prev m generations
        self.maxDeltaStats = [Population._LARGEVAL] * len(self.prevMGenStats[0])
        # CONVERGENCE CRITERIA
        self.Gmin = 10    # run 10 generations, minimum
        self.Gmax = 50
        self.convgStrategy = DeDriver._availConvgStrats[1]
        # Use 'tol' to define the tolerance for the convergence strategy,
        # use 'deltaTol' to define the tolerance for the time-stability of the
        # convergence.
        self.tol = 1.0E-3
        self.deltaTol = 1.0E-4
        self.truncFrac = 0.15
        # logistical files
        self.lockFName = 'dedriver.loc'
        self.fcnEvalFName = 'de_fcn.dat'
        # CAF job parameters
        if email != '':
            self.emailAddr = email
        else:
            self.emailAddr = os.environ['USER']+'@fnal.gov'
        self.cafFarm = random.choice(DeDriver._availFarms.keys())
        # Get the absolute path of the working directory;
        # NOTE: DO NOT USE THIS ON THE dCAFs!! THE GLIDE-CAFs ESPECIALLY USE
        # SUPER-LONG PATHNAMES WHICH ARE LIABLE TO SCREW STUFF UP (e.g. Fortran)
        self.workDir = os.getcwd()
        self.cafBuildDir = os.path.join(self.workDir, 'build_de')
        self.scriptsDir = os.path.join(self.workDir, 'scripts')
        self.cafQueue = 'short'
        self.cafGroup = 'common'
        self.cafOutLocale = '' # e.g. 'galyardt@fcdflnx9.fnal.gov:~/'
        self.cafDataSrc = 'None' # should only be changed for SAM access
        self.cafNEvalsPerSeg = 1
        self.cafCmd = ''
        # This is only needed for monolithic tarball submission:
        self.cafFcnTarball = ''
        # Split tarball submission:
        self.cafFcnTarballUrl = ''
        # This is the submission tarball in the split framework; should contain
        # the master script and a nested 'user' tarball (self.cafNestedTarball)
        self.cafSubTarball = 'deDriveCaf.tgz'
        # This is the 'user' tarball in the split tarball framework; should contian
        # a pickled DeDriver instance (the state file), plus
        # any other scripts required by DeDriver
        self.cafNestedTarball = 'deFiles.tgz'
        self.cafSubmitScript = os.path.join(self.scriptsDir, 'submit_DE')
        self.cafBuildTarScript = os.path.join(self.scriptsDir, 'build_user_tarball')
        self.cafMasterScript = 'deCaf.sh'
        self.cafSegmentMap = {}
        self.cafSegment = 0
        self.cafOutFName = ''
        self.cafSubmitTime = 0.0 # in seconds since the epoch (a float)
        self.cafJobCompltFrac = 0.2
        self.local = False
        self.runLocalDir = 'de_local'
        # Web Monitoring Parameters
##        self.monitorNode = ''
##        self.monitorDir = ''
        self.monitorLoc = '' # e.g. 'user@node.fnal.gov:~/'
        self.monitorUrl = ''
        self.verbose = False
        self.debug = False
        self.logFName = 'de.log'
        

    def __repr__(self):
        """ Overload __repr__() so that `deDriver` will print the algorithm
        state to STDOUT.
        """
        stats = self.population.getStats(self.truncFrac)
        outStr = ("""DeDriver state:
        Generation: %s
        best trial solution's cost value: %s
        best trial solution: %s\n
        worst trial solution's cost value: %s
        worst trial solution: %s\n
        mean cost value: %s
        standard deviation of cost values: %s        
        fractional tolerance: %s
        chi-square: %s
        ndf: %s
        np_trunc: %s\n
        Parameter mean values:\n%s\n
        Covariance Matrix for np_trunc members:\n%s
        Name of dCAF last used: %s
        Time of last CAF job submission: %s""" %
                  (self.population.generation,
                   stats[0],
                   self.population[self.population.ibest].x,
                   stats[1],
                   self.population[self.population.iworst].x,
                   stats[2],
                   stats[3],
                   stats[4],
                   stats[5],
                   stats[6],
                   self.population.getNpTrunc(self.truncFrac),
                   stats[7],
                   self.population.getCovMatRepr(stats[8]),
                   self.cafFarm,
                   time.ctime(self.cafSubmitTime)))
        return outStr


    def printFull(self):
        """ Return a string with the full instance state.
        """
        ph1Cr = 0.0
        ph2Cr = 0.0
        if self.population.nDim < self.dimThresh:
            ph1Cr = self.phase1LoDimCr
            ph2Cr = self.phase2LoDimCr
        else:
            ph1Cr = self.phase1HiDimCr
            ph2Cr = self.phase2HiDimCr
        outStr = """DeDriver Algorithmic Parameters:
        --------------------------------
        Pareto: %s
        Generation: %s
        Population size (Np): %s
        Population Depth (D): %s
        Hard Constraints: %s \n
        Soft Constraints: %s \n
        Initial Ranges: %s \n
        Save Initial Population: %s
        Load Initial Population: %s
        Phase 1 Cr: %s
        Phase 1 F:  %s
        Phase 1 Strategy: %s
        Phase 2 Cr: %s
        Phase 2 F:  %s
        Phase 2 Strategy: %s
        Xi: %s
        trig-prob(m_t): %s

        DeDriver Convergence Parameters:
        --------------------------------
        Gmin: %s
        Gmax: %s
        Convergence strategy: %s
        m: %s
        Tolerance: %s
        Max tolerance for last 'm' generations: %s
        ChiSquare truncation fraction: %s

        DeDriver CAF Parameters:
        ------------------------
        Objective Function Command: %s
        Objective Function Tarball URL: %s
        Queue: %s
        Group: %s
        Output location: %s
        Number of function evaluations per CAF segment: %s
        Minimum fraction of completed CAF segments: %s

        DeDriver Misc Logistical Parameters:
        ------------------------------------
        Email: %s
        Work dir: %s
        State file: %s
        Initial Population file: %s
        Lock file: %s
        Objective function output file: %s\n
        """ % (self.pareto, self.population.generation, self.population.np,
               self.population.nDim, ', '.join(self.population.hardConstraints),
               ', '.join(self.population.softConstraints),
               `self.population.initRange`, self.saveInitPop, self.initFromFile,
               ph1Cr, self.phase1F, self.phase1Strategy,
               ph2Cr, self.phase2F, self.phase2Strategy, self.xi, self.population.mt,
               self.Gmin, self.Gmax, self.convgStrategy, self.m, self.tol,
               self.deltaTol, self.truncFrac,
               self.cafCmd, self.cafFcnTarballUrl, self.cafQueue, self.cafGroup,
               self.cafOutLocale, self.cafNEvalsPerSeg, self.cafJobCompltFrac,
               self.emailAddr, self.workDir, self.stateFName, self.popFile,
               self.lockFName, self.fcnEvalFName)
        return outStr


    def getKey(self):
        """ Return a string which can be used as a dictionary key for the
        current DeDriver instance.  Form: 'de_state:5'
        """
        return ':'.join((self.stateFName.split('.')[0],
                         `self.population.generation`))


    def saveState(self):
        """ Save the algorithm state to a gzipped file.
        """
        try:
            gfile = gzip.open(self.stateFName, 'wb')
            cPickle.dump(self, gfile)
            gfile.close()
        except(IOError, gzip.zlib.error, cPickle.PicklingError):
            return False
        else:
            return True


    def shelveState(self):
        """ Save the algorithm state to a gzipped shelve file;
        this state can be retrieved later by keyed DB lookup (see
        Python's shelve module).  Currently, the DBM file is not
        compressed.
        """
        fname = self.stateFName.split('.')[0] + '_history.dbm'
        try:
            stateDB = shelve.open(fname)
            stateDB[self.getKey()] = self
            stateDB.close()
        except:
            return False
        else:
            return True


    def pause(self):
        """ Create the lock file, pausing the algorithm.
        """
        Pause(self.lockFName)


    def sendMail(self, subject, message):
        """ Send mail to the user [wrapper for DeDriver.SendMail()]
        """
        msg = ("DeDriver State File: %s\nDeDriver Function Command: %s\n\n"
               % (self.stateFName, self.cafCmd))
        msg += '-' * 30 + '\n\n' + message + '\n\n' + '-' * 30 + '\n\n' + `self`
        return SendMail(self.emailAddr, subject, msg)
    

    def writeLog(self, message):
        """ Write a message to the log file.
        """
        #msg = '\n\n' + '-' * 80 + '\n\n'
        msg = '\n\n'
        msg += '<<' + time.ctime(time.time()) + '>>' + '\n'
        #msg += message + `self`
        msg += message
        WriteLog(self.logFName, msg)        

        
    def getCafOutFileName(self, generation, segment):
        """ Return a string var containing the CAF output filename for the
        given segment; assume that this file is gzipped.
        """
        return ('gen' + `generation` +
                '_seg' + `segment` + '.dat.gz')


    def getCafStateFileName(self):
        """ Return a string defining the state file name to be used
        for CAF jobs.  Hopefully, this will save the user from
        inadvertently overwriting the current state file if they
        should happen to expand a CAF output tarball in the DeDriver
        working directory.
        """
        return ('gen' + `self.population.generation` + '_' + self.stateFName)


    def setCafOutFileName(self):
        """ Set the CAF output filename for the current CAF segment.
        """
        self.cafOutFName = self.getCafOutFileName(self.population.generation+1,
                                                  self.cafSegment)


    def setDEStrategy(self, xi=0.7, dimThresh=6,
                      ph1LDCr=0.5, ph1HDCr=0.7, ph1F=0.8,
                      ph2LDCr=0.1, ph2HDCr=0.2, ph2F=0.6, mt=0.05,
                      ph1StratStr='rand-trig1', ph2StratStr='best1',
                      sade=False, zeta=0.5):
##                      ph1Strat=('rand', 1), ph2Strat=('best', 2)):
        """ Set the parameters for the DE strategy.
        """
        # Ugly implementation of DE strategy
        ph1Strat = Population._strategies[-2]
        ph2Strat = Population._strategies[1]
        if ph1StratStr == 'rand-trig1':
            ph1Strat = ('rand-trig',1)
        elif ph1StratStr == 'rand1':
            ph1Strat = ('rand',1)
        elif ph1StratStr == 'best1':
            ph1Strat = ('best',1)
        elif ph1StratStr == 'best2':
            ph1Strat = ('best',2)
        elif ph1StratStr == 'rand-sort1':
            ph1Strat = ('rand-sort',1)
        elif ph1StratStr == 'best-trig1':
            ph1Strat = ('best-trig',1)
        #
        if ph2StratStr == 'rand-trig1':
            ph2Strat = ('rand-trig',1)
        elif ph2StratStr == 'rand1':
            ph2Strat = ('rand',1)
        elif ph2StratStr == 'best1':
            ph2Strat = ('best',1)
        elif ph2StratStr == 'best2':
            ph2Strat = ('best',2)
        elif ph2StratStr == 'rand-sort1':
            ph2Strat = ('rand-sort',1)
        elif ph2StratStr == 'best-trig1':
            ph2Strat = ('best-trig',1)
        #
        self.xi = xi
        self.dimThresh = dimThresh
        # Set Phase 1 (Exploration) parameters
        if ph1Strat in Population._strategies:
            self.phase1Strategy = ph1Strat
        else:
            self.phase1Strategy = ('rand', 1) 
        self.phase1LoDimCr = ph1LDCr
        self.phase1HiDimCr = ph1HDCr
        self.phase1F = ph1F
        if self.population.nDim < dimThresh:
            self.population.crossProb = self.phase1LoDimCr
        else:
            self.population.crossProb = self.phase1HiDimCr
        self.population.strategy = self.phase1Strategy
        self.population.F = self.phase1F
        self.population.mt = mt
        # Set Phase 2 (convergence) parameters
        if ph2Strat in Population._strategies:
            self.phase2Strategy = ph2Strat
        else:
            self.phase2Strategy = ('best', 2)
        self.phase2LoDimCr = ph2LDCr
        self.phase2HiDimCr = ph2HDCr
        self.phase2F = ph2F
        self.sade = sade
        self.zeta = zeta
        

    def setConstraints(self, softConstr=[], hardConstr=[]):
        """ Set the constraints.  Assume they are well-formed.
        """
        self.population.hardConstraints = tuple(hardConstr)
        self.population.softConstraints = tuple(softConstr)
        

    def setConvergenceStrategy(self, strat, gmin=10, gmax=50,
                               tol=1.0E-3, dtol=1.0E-4, m=5, truncFrac=0.15):
        """ Set the convergence strategy for the DE algorithm.
        """
        if gmax < gmin and gmin > 0:
            self.Gmax = gmin
        elif gmax > 0:
            self.Gmax = gmax
        self.Gmin = gmin
        if strat not in DeDriver._availConvgStrats:
            self.convgStrategy = DeDriver._availConvgStrats[1]
        else:
            self.convgStrategy = strat
        if tol > 0:
            self.tol = tol
        if dtol > 0:
            self.deltaTol = dtol
        if m < self.Gmax and m > 0:
            self.m = m
        else:
            # setting self.m = self.GMax will effectively
            # ensure that the algorithm will not converge before
            # reaching Gmax generations.
            self.m = self.Gmax
        if truncFrac > 0.0 and truncFrac <= 1.0:
            self.truncFrac = truncFrac
        self.prevMGenStats = [(Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL)] * self.m
        self.maxDeltaStats = [Population._LARGEVAL] * len(self.prevMGenStats[0])


    def setCafParms(self, cmd, neval=1, url='', tarfile='', queue='short',
                    group='common', outloc='', dataSrc='None', frac=0.2,
                    local=False, runLocalDir=''):
        """ Set CAF parameters.
        """
##         ms='./scripts/deCaf.sh',
##                     ss='./scripts/submit_DE',
##                     bs='./scripts/build_user_tarball',
        # Error handling!
        if not cmd:
            raise ValueError, 'CAF Command not set!'
        if not outloc:
            # Construct the output location from the current nodename, username,
            # and working directory
            # Should probably check that fcpd is running on the local host...
            self.cafOutLocale = (os.environ['USER'] + '@' + os.environ['HOST'] +
                                 ':' + self.workDir)
        else: 
            self.cafOutLocale = outloc
        self.cafFcnTarball = tarfile
        self.cafFcnTarballUrl = url
        self.cafCmd = cmd
        self.cafNEvalsPerSeg = neval
        self.cafQueue = queue
        self.cafGroup = group
        self.cafDataSrc = dataSrc
##         self.cafMasterScript = ms
##         self.cafSubmitScript = ss
##         self.cafBuildTarScript = bs
        self.cafJobCompltFrac = frac
        self.local = local
        if self.local:
            self.cafFarm = 'local'
            self.runLocalDir = runLocalDir
        # now build self.cafSegmentMap, which need only be done once
        # per minimization run.
        nSeg = int(math.ceil(float(self.population.np) / float(self.cafNEvalsPerSeg)))
        mPtr = 0
        upper = 0
        for i in range(1,nSeg+1):
            if mPtr+self.cafNEvalsPerSeg < len(self.population):
                upper = mPtr+self.cafNEvalsPerSeg
            else:
                upper = len(self.population)
            self.cafSegmentMap[i] = tuple(range(mPtr, upper))
            mPtr += self.cafNEvalsPerSeg
        if self.debug:
            print '\nCaf Segment Map:  %s\n' % `self.cafSegmentMap`


    def initialize(self, ir=[], gaussParms=[]):
        """ Initialize the Population.
        """
        # Expect initRange to be a list of tuples, each tuple
        # of the form (index, lobnd, hibnd); assume that
        # the tuple values are valid.
        # NOTE: if len(ir) < self.population.nDim, the remainder
        # of the parameters will take their initial ranges from
        # the default range, defined in Population.__init__()
        #if len(self.popFile) == 0:
        #
        # SADE --> INTIALIZE CONTROL PARM POPULATION!!
        #
        # NOTE: Gaussian parms will be used only for those
        #       decision vector elements for which they are
        #       defined; other elements will use a random uniform
        #       distribution using the bounds set in 'ir'.
        #
        if len(gaussParms) > 0:
            self.population.initGauss = True
        if not self.initFromFile:
            for item in ir:
                self.population.initRange[item[0]] = tuple(item[1:])
            for item in gaussParms:
                self.population.initGaussParms[item[0]] = tuple(item[1:])
            self.population.initialize()
            self.ctrlPop.initialize()
            if self.saveInitPop:
                try:
                    pfile = gzip.open(self.popFile, 'wb')
                    cPickle.dump(self.population, pfile)
                    cPickle.dump(self.ctrlPop, pfile)
                    pfile.close()
                except:
                    print """%s: Warning: Failed to save initial population to file '%s'.
                    Continuing... """
        else:
            # load the population from a file
            # population should be in a gzipped pickle file, not a DBM (shelve)
            try:
                pfile = gzip.open(self.popFile, 'rb')
                # Weak design: order of objects in file should not be
                # important!
                tmp = cPickle.load(pfile)
                tmpCtrlPop = cPickle.load(pfile)
                pfile.close()
                if not isinstance(tmp, Population):
                    raise (TypeError, 'Expected a pickled instance of class Population,'
                           +' found a pickled object of type %s' % tmp.__class__)
                if not isinstance(tmpCtrlPop, SAPopulation):
                    raise (TypeError, 'Expected a pickled instance of class '
                           +'SAPopulation, found an object of type %s'
                           % tmpCtrlPop.__class__)
                if tmp.nDim != self.population.nDim:
                    raise ValueError, 'Expected D=%s, found D=%s' % (
                        self.population.nDim, tmp.nDim)
                if tmp.np != self.population.np:
                    raise ValueError, 'Expected Np=%s, found Np=%s' % (
                        self.population.np, tmp.np)
                if tmpCtrlPop.np < 1:
                    raise (ValueError, 'Need at least ONE member in the'
                           +'control parm population.')
                # If we've gotten this far, it should be safe to use the
                # loaded Population object.
                self.population = tmp
                self.ctrlPop = tmpCtrlPop
            except ((IOError, gzip.zlib.error,
                     cPickle.PicklingError, TypeError, ValueError), sys.exc_info()[1]):
                print """%s: ERROR: %s.\nFailed to load initial population from file %s.
                Exiting...
                """ % (sys.argv[0], data, self.popFile)
                sys.exit(1)


    def updateMonitor(self):
        """ Update the web-based monitor's data.
        """
        msg = '\n' + `self` + '\n'
        self.writeLog(msg)
        # print "  Monitor URL for this job: %s" % deDriver.monitorUrl


    def getBestSubmit(self):
        """ Get the best-submit dCAF for the current username.
        """
        if self.local:
            self.cafFarm = 'local'
            return
        global farmdict
        try:
            # FAIL!
            raise ValueError
            farmdict = {}
            bestfarm = {'name':'', 'group':'', 'score':-500.0, 'free_vms':0}
            # create the XML parser object
            xp = xml.parsers.expat.ParserCreate()
            xp.returns_unicode = False
            xp.StartElementHandler = XmlStartElement
            # Open the URL for the XML-based best-submit CGI script.
##             urlstream = urllib.urlopen(
##                 'http://dcafmon.fnal.gov/cgi-bin/dcafmon/xml/best_submit.py?user=%s'
##                 % os.environ['USER'])
            urlstream = urllib.urlopen(
                'http://cdfcaf.fnal.gov/cgi-bin/dcafmon/xml/best_submit.py?user=%s'
                % os.environ['USER'])
            # Get the actual data.
            bsxml = urlstream.read()
            xp.Parse(bsxml, True)
            #
            groups = ['common', 'group_CDFgeneric']
            if self.cafGroup == 'MCprod':
                groups.append('group_MCprod')
            # Loop through the farms
            for farm in farmdict.keys():
                # Consider only those farms which are 'available' (ie. they
                # are dCAFs and the desired queue is available on each farm)
                if (farm in DeDriver._availFarms.keys() and
                    self.cafQueue in DeDriver._availFarms[farm][1]):
                    for grp in farmdict[farm]:                        
                        # farmdict[farm] is a list of dictionaries.
                        # grp will be a dictionary of attributes for the
                        # current account group for farm 'farm'
                        # --> grp['name'] is the account group name
                        if (grp['name'] in groups and
                            float(grp['score']) > bestfarm['score'] and
                            int(grp['free_vms']) >= bestfarm['free_vms']):
                            bestfarm['name'] = farm
                            bestfarm['group'] = grp['name']
                            bestfarm['score'] = float(grp['score'])
                            bestfarm['free_vms'] = int(grp['free_vms'])
            self.cafFarm = bestfarm['name']
        except:
            # Problems with the best submit script, so choose
            # the farm randomly
            farm = ''
            farmKeys = DeDriver._availFarms.keys()
            if len(farmKeys) == 2:
                # only one choice; choose the one that's not 'local'
                for f in farmKeys:
                    if f != 'local':
                        farm = f
                        break
            elif len(farmKeys) > 2:
                # we have choices!
                while True:
                    farm = random.choice(farmKeys)
                    if (farm != self.cafFarm and farm != 'local' and
                        self.cafQueue in DeDriver._availFarms[farm][1]):
                        # Either we chose the same farm as used in the last
                        # job, or the requested queue is not available for
                        # this farm, so generate a new one.
                        break
##             else:
##                 raise ValueError, 'No valid farms to choose from!'
            if farm == '' or farm == 'local':
                raise ValueError, 'No valid farms to choose from!'
            self.cafFarm = farm
        # OTHER ERROR HANDLING? .k5login checking?

        
    def setupCafJob(self):
        """ Setup the CAF job.
        """
        self.getBestSubmit()
        self.cafSubTarball = 'deDriveCaf_gen%s.tgz' % self.population.generation
        self.cafNestedTarball = 'deFiles_gen%s.tgz' % self.population.generation
        # Build the job submission tarball.
        # MONOLITHIC Tarball submission:
        #   The submission tarball should contain both the
        #   user's tarball (which contains everything necessary
        #   for the objective function eval) and the DeDriver files
        #   (the pickled DeDriver instance + necessary scripts)
        # SPLIT TARBALL FRAMEWORK:
        #   The submission tarball should *only* contain the
        #   DeDriver files (state file + scripts) + the master
        #   script which handles downloading and expanding the
        #   user's tarball (analogous to mcProduction/scripts/cafmc.sh).
        #   ** The tarball which contains all the files for function evaluation
        #   ** should be housed on a public web server (URL provided by user in
        #   ** DeDriver setup mode).
        #--
        if self.cafFcnTarballUrl != '':
            # Assume the user want's to use the split tarball framework.
            # Copy the files needed for the CAF job to the build directory
            if os.path.isdir(self.cafBuildDir):
                # build directory exists; remove it, ignoring errors
                try:
                    shutil.rmtree(self.cafBuildDir, True)
                except:
                    pass
            os.mkdir(self.cafBuildDir)
            # Copy the Python distribution to the build dir            
            # REALLY need a check on the Python version available
            # from the system here!!
            try:
                shutil.copytree(os.path.join(self.workDir, 'bin'),
                                os.path.join(self.cafBuildDir, 'bin'))
                shutil.copytree(os.path.join(self.workDir, 'lib'),
                                os.path.join(self.cafBuildDir, 'lib'))
                shutil.copytree(os.path.join(self.workDir, 'include'),
                                os.path.join(self.cafBuildDir, 'include'))
            except:
                # STOP STOP STOP
                msg = """DeDriver::setupCafJob() Failed to find the proper
                Python distribution; job submission is impossible.  The
                DeDriver algorithm will be paused (lock file: %s.)\n""" % (
                    self.lockFName)
                if self.debug:
                    print msg
                else:
                    self.writeLog('ERROR -- DeDriver::setupCafJob()\n'+msg)
                    self.sendMail('DeDriver::setupCafJob() ERROR', msg)
                self.pause()
                return False
                return True
            # May want to be more flexible about the path to the scripts...
            deScript = os.path.basename(sys.argv[0])            
            try:
                shutil.copy(os.path.join(self.workDir, deScript),
                            self.cafBuildDir)
            except:
                msg = """DeDriver.setupCafJob() failed to find the
                necessary shell scripts; job submission is impossible.
                The DeDriver algorithm will be paused."""
                if self.debug:
                    print msg
                else:
                    self.writeLog('ERROR -- DeDriver.setupCafJob()\n'+msg)
                    self.sendMail('DeDriver::setupCafJob() ERROR', msg)
                self.pause()
                return False
            try:
                shutil.copy(os.path.join(self.workDir, self.stateFName),
                            os.path.join(self.cafBuildDir,
                                         self.getCafStateFileName()))
            except:
                msg = """DeDriver.setupCafJob() failed to find the DE
                state file %s; job submission is impossible.
                The DeDriver algorithm will be paused (lock file: %s)""" % (
                    self.stateFName, self.lockFName)
                if self.debug:
                    print msg
                else:
                    self.writeLog('ERROR -- DeDriver.setupCafJob()\n' + msg)
                    self.sendMail('DeDriver::setupCafJob() ERROR', msg)
                self.pause()
                return False
            # Copy the fcp stuff to the build directory...
            try:
                shutil.copy(os.path.join(self.scriptsDir, 'setup_fcp.sh'),
                            self.cafBuildDir)
                shutil.copytree(os.path.join(self.workDir, 'fcp'),
                                os.path.join(self.cafBuildDir, 'fcp'))
                shutil.copytree(os.path.join(self.workDir, 'fcslib'),
                                os.path.join(self.cafBuildDir, 'fcslib'))
            except:
                msg = """DeDriver.setupCafJob() failed to copy the fcp-
                related files to the build directory."""
                if self.debug:
                    print msg
                else:
                    self.writeLog('ERROR --DeDriver.setupCafJob()\n' + msg)
                    self.sendMail('DeDriver.setupCafJob() ERROR', msg)
                self.pause()
                return False
            # any other scripts to copy to the build dir?
            # (The master logistical script run on the CAF, deCaf.sh, will be copied
            # automagically by build_user_tarball).
            buildCmd = ('%s -b %s -s %s -t %s -u %s' %
                        (self.cafBuildTarScript,
                         self.cafBuildDir,
                         os.path.join(self.scriptsDir, self.cafMasterScript),
                         self.cafSubTarball,
                         self.cafNestedTarball))
            if self.debug:
                print 'DeDriver::setupCafJob(): buildCmd = %s' % buildCmd
            child = popen2.Popen4(buildCmd)            
            retVal = child.wait()
            if retVal != 0:
                # send mail to user with the error message
                subject = 'DeDriver Job Submission Error'
                message = """DeDriver Job Submission Error:
                Problems with the command\n\n%s\n\nThe error:\n\n%s""" % (
                    buildCmd, child.fromchild.read())
                child.tochild.close()
                child.fromchild.close()
                self.writeLog(subjet+'\n\n'+message)
                self.sendMail(subject, message)
                # Write the lock file
                self.pause()
                return False
            else:
                # Remove the nested tarball build directory, ignoring errors
                child.tochild.close()
                child.fromchild.close()
                shutil.rmtree(self.cafBuildDir, True)
                return True
        else:
            # use the Monolithic tarball framework
            # Disable, for now
            return False


    def submitCafJob(self):
        """ Submit the CAF job.
        """
        # Check that the proper environment is setup
        # Build the tarball, etc.
        if not self.setupCafJob():
            return False
        else:
            # Set the CAF group to use (could be 'common' or 'MCProd')
            group = self.cafGroup
            # Is this bit necessary? or will 'common' be translated
            # correctly by the dCAFs which have a different name for
            # that group?
            if self.cafGroup == 'common':
                group = DeDriver._availFarms[self.cafFarm][0]
            #
            # Construct the DeDriver command for caf mode; note that
            # we have to use a raw string with the escaped '$' in order
            # to make sure the CAF segment number gets passed.
            deScript = os.path.basename(sys.argv[0])
            #deCmd = r'./%s --state-file %s -j $ caf' % (deScript, self.stateFName)
            deCmd = r'./%s --state-file %s -j $ caf' % (deScript,
                                                        self.getCafStateFileName())
            segments = '1:%s' % len(self.cafSegmentMap)
            # Construct the submission command
            submitCmd = '%s' % self.cafSubmitScript
            if self.debug:
                submitCmd += ' -v debug_only'
            if self.verbose:
                submitCmd += ' -v yes'
            submitCmd += r' -n -c %s -f %s -g %s -q %s -S %s -s %s -t %s -w %s -o %s' % (
                self.cafNestedTarball, DeDriver._availFarms[self.cafFarm][2],
                group, self.cafQueue,
                self.cafMasterScript, segments, self.cafSubTarball,
                self.cafFcnTarballUrl, self.cafOutLocale)
            if self.local:
                submitCmd += r' -R %s' % (self.runLocalDir)
            submitCmd += r' %s' % (deCmd)
            if self.debug or self.verbose:
                self.writeLog('CAF Job Submission command:\n\n%s\n\n' % submitCmd)
                #self.cafSubmitTime = time.time()
                # NOTE: submission tarball is *not* removed.
                #return True
            # Submit the job
            child = popen2.Popen4(submitCmd)
            # Do we need / want to capture the JOB ID?
            # Close the child's standard input
            child.tochild.close()
            retVal = child.wait()
            # Trap the output of CafSubmit, send it to the user
            # if CafSubmit terminates with an error condition.
            if retVal != 0:
                msg = child.fromchild.read()
                subj = 'DeDriver Job Submission Error'
                self.writeLog(subj+'\n\n'+msg)
                self.sendMail(subj, msg)
                # Write the lock file
                self.pause()
                return False
            if self.debug:
                self.writeLog('\nOutput of Job submission cmd:\n%s'
                              % child.fromchild.read())
            child.fromchild.close()
            # Remove the submission tarball...
            if not self.debug:
                try:                
                    os.remove(self.cafSubTarball)
                except:
                    pass
            # Lastly, record the submission time in seconds since the epoch
            self.cafSubmitTime = time.time()
            return True
        

    def saveStatistics(self):
        """ Save the stats of the current population.
        """
        # Find the best and worst population members
        # Best and worst members are found in population.getStats(),
        # which is called below -- should happen for every generation.
        # Shift the current values back one index
        #for i in range(self.m-1):
        del self.prevMGenStats[0:1]
        self.prevMGenStats.append(self.population.getStats(self.truncFrac)[:7])
        # Loop over the statistics defined in Population::getStats()
        for j in range(len(self.prevMGenStats[0])):
            maxDelta = -Population._LARGEVAL
            delta = -Population._LARGEVAL
            #for i in range(self.m-1):
            for i in range(len(self.prevMGenStats)-1):
                # Take the abs of the difference in statistic j
                # for generation i and generation i+1
                delta = abs(self.prevMGenStats[i][j] -
                            self.prevMGenStats[i+1][j])
                if delta > maxDelta:
                    maxDelta = delta
            self.maxDeltaStats[j] = maxDelta


    def converged(self):
        """ Check for convergence of the DE algorithm.
        """
        viCtr = 0
        npTrunc = self.population.getNpTrunc(self.truncFrac)
        chi2ndf = self.prevMGenStats[-1][5] / float(self.prevMGenStats[-1][6])
        maxChi2ndf = self.maxDeltaStats[5] / float(self.prevMGenStats[-1][6])
        # Count the viable Member objects
        for key in self.population.keys():
            if (self.population[key].isValid == True and
                self.population[key].deltaConstr == 0.0):
                viCtr += 1
        # Check the fractional tolerance of the current generation,
        # as well as the fractional tolerance of the last m generations
        if self.population.generation < self.Gmin:
            return False
        elif self.population.generation >= self.Gmax:
            return True
        elif (self.convgStrategy == 'fractol' and
              viCtr >= 2 and
              self.prevMGenStats[-1][4] < self.tol and
              self.maxDeltaStats[4] < self.deltaTol):
            return True
        elif (self.convgStrategy == 'chisq' and
              viCtr >= npTrunc and
              chi2ndf < (1.0 + self.tol) and
              maxChi2ndf < self.deltaTol):
            return True
        else:
            return False


    def loadMembers(self, filename, gen):
        """ Load Member objects from a file without a priori
        knowledge of how many Member objects will be in the file. Expect the
        Member objects to be from generation 'gen'.
        """
        mainSubPop = {}
        ctrlSubPop = {}
        uniq = True
        try:
            gfile = gzip.open(filename, 'rb')
        except IOError:
##             return False
            return (mainSubPop, ctrlSubPop)
##         mbrCtr = 0
        while True:
            uniq = True
            try:
                tmpObj = cPickle.load(gfile)
                # Check the class name of the object just retrieved
                if isinstance(tmpObj, SAMember):
                    if (tmpObj.generation != gen or
                        tmpObj.nDim != self.ctrlPop.nDim):
                        raise ValueError
                    for key in ctrlSubPop.keys():
                        if not tmpObj.isUnique(ctrlSubPop[key]):
                            # we've loaded a non-unique individual (w.r.t. decision
                            # vector), but this individual may have performed better
                            # than it's twin(s).
                            uniq = False
                            if tmpObj < ctrlSubPop[key]:                                
                                # the newer individual dominates the old (strictly), so
                                # replace the old with the new
                                if tmpObj.getKey() != key:
                                    del ctrlSubPop[key]
                                uniq = True
                            break
                    if uniq:
                        ctrlSubPop[tmpObj.getKey()] = tmpObj
                elif isinstance(tmpObj, Member):
                    if (tmpObj.generation != gen or
                        tmpObj.nDim != self.population.nDim):
                        raise ValueError
##                     mbrCtr += 1
                    mainSubPop[tmpObj.getKey()] = tmpObj
                else:
                    raise TypeError, 'Attempted to load an object of type other than Member or SAMember'
            except EOFError:
                break
            except ValueError:
                continue
            except TypeError:
                # Found an object which is not of type Member;
                # can either close the file, or look for more
                # Member objects...
                continue
            except cPickle.UnpicklingError:
                continue
        gfile.close()
        return (mainSubPop, ctrlSubPop)


    def loadNextGen(self, nseg):
        """ Load the members of the next generation from the output files of
        the last CAF job. Expect at least nseg segments to have returned.
        """
        nextGen = self.population.generation + 1
        nextGenPop = copy.deepcopy(self.population)
        nextGenPop.setGeneration(nextGen)
        tmpCtrlPop = {}
        # Get the Member objects for this CAF segment;
        # the file will be gzipped
        success = False
        mbrCtr = 0
        segCtr = 0
        loadedFiles = []
        lastGenFiles = []
        for seg in self.cafSegmentMap.keys():            
            fname = self.getCafOutFileName(nextGen, seg)
            if nextGenPop.generation > 1:
                lastGenFiles.append(self.getCafOutFileName(nextGen-1, seg))
            mainSubPop, ctrlSubPop = self.loadMembers(fname, nextGen)
            n = len(mainSubPop)
            if n > 0:
                loadedFiles.append(fname)
                segCtr += 1
                mbrCtr += n
                nextGenPop.memberColl.update(mainSubPop)
                if self.sade:
                    for mkey in ctrlSubPop.keys():
                        if mkey in tmpCtrlPop.keys():
                            # This key already loaded! (can happen because
                            # we're sampling from the Pareto front in a
                            # stochastic fasion)
                            uniq = ctrlSubPop[mkey].isUnique(tmpCtrlPop[mkey])
                            if ctrlSubPop[mkey] < tmpCtrlPop[mkey]:
                                # The new SAMember object dominates the old
                                # object with the same key, so replace
                                # the old with the new (regardless of whether
                                # the new is actually unique w.r.t. the old).
                                tmpCtrlPop[mkey] = ctrlSubPop[mkey]
                            elif uniq and ctrlSubPop[mkey] == tmpCtrlPop[mkey]:
                                # we have two non-dominated individuals,
                                # create a new dict key for the SAMember
                                # just loaded if the new guy is unique.
                                keylist = copy.deepcopy(self.ctrlPop.keys())
                                keylist.extend(tmpCtrlPop.keys())
                                maxkey = max(keylist)
                                ctrlSubPop[mkey].popIndex = maxkey+1
                                tmpCtrlPop[maxkey+1] = ctrlSubPop[mkey]
                        else:
                            # check that all members of ctrlSubPop are unique
                            # w.r.t. tmpCtrlPop
                            uniq = True
                            for tkey in tmpCtrlPop.keys():
                                if not ctrlSubPop[mkey].isUnique(tmpCtrlPop[tkey]):
                                    uniq = False
                                    if ctrlSubPop[mkey] < tmpCtrlPop[tkey]:
                                        # the newer individual dominates the old (strictly), so
                                        # replace the old with the new
                                        if mkey != tkey:
                                            del tmpCtrlPop[tkey]
                                        uniq = True
                                    break
                            if uniq:
                                tmpCtrlPop[mkey] = ctrlSubPop[mkey]
        if segCtr >= nseg:
            self.population = nextGenPop
            if self.sade:
                self.ctrlPop.memberColl.update(tmpCtrlPop)
                self.ctrlPop.setGeneration(nextGen)
        # remove the files from the last generation
        if not self.debug:
            for fname in lastGenFiles:
                try:
                    os.remove(fname)
                except:
                    pass
        return mbrCtr


    def cafJobComplete(self):
        """ Check whether the last CAF Job submitted has completed
        """
        # deDriver.cafJobComplete() should return True if we have all
        # output files from the CAF job or we have exceeded the real-time
        # limit for the queue used (in this case, send an error report to
        # the user; proceed with the next generation's CAF job only if at
        # least one CAF job segment came back okay).
        nextGen = self.population.generation + 1
        elapsedTime = time.time() - self.cafSubmitTime
        # The safety factor is additive because we're trying to account
        # for variations in file transfer times; we want the same safety
        # factor no matter the length of the queue real time limit.
        queueTime = DeDriver._safetyFac + DeDriver._stdQueueData[self.cafQueue]
        fileCtr = 0
        nseg = len(self.cafSegmentMap)
        for seg in self.cafSegmentMap.keys():
            fname = self.getCafOutFileName(nextGen, seg)
            if os.path.exists(fname):
                fileCtr += 1
        if fileCtr == nseg:
            return True
        elif (fileCtr >= self.cafJobCompltFrac * nseg and
              elapsedTime > queueTime):
            # Some fraction of our total number of segments have returned,
            # but there should have been plenty of time for the remaining
            # jobs to transfer their output.
            # Assume a fraction of the jobs crashed and will not give
            # any output -- continue with the minization alg.
            return True
        elif (fileCtr < self.cafJobCompltFrac * nseg and
              elapsedTime > queueTime):
            # No files have returned yet and
            # too much time has elapsed since job submission,
            # create the lockfile and send email to user.
            self.pause()
##            bakfname = self.stateFName + '.bak'
##            msg = """DeDriver Error: \t Low output file count for elapsed time
##            greater than the queue real-time limit plus a safety factor
##            (%s seconds).\n
##            The Current algorithmic state will be saved to %s and what files
##            did return from the CAF job will be loaded to %s and used when the
##            algorithm is resumed by the removal of the lock file (%s).""" % (
##                DeDriver._safetyFac, bakfname, self.stateFName, self.lockFName)
            msg = """DeDriver Error: \t Low output file count for elapsed time
            greater than the queue real-time limit plus a safety factor
            (%s seconds).  The algorithm may be resumed by removing the lock
            file (%s)""" % (DeDriver._safetyFac, self.lockFName)
            self.writeLog('ERROR -- DeDriver.py CAF Job Incomplete\n'+msg)
            self.sendMail('ERROR -- DeDriver.py CAF Job Incomplete', msg)
            # Backup the last generation's state
##            try:
##                shutil.copy(os.path.join(self.workDir, self.stateFName),
##                            os.path.join(self.workDir, bakfname))
##            except:
##                pass
##            # Load what new Member objects we can find
##            nloaded = self.loadNextGen(1)
##            self.writeLog('\nLoaded a total of %s CAF output files of generation %s.\n' % (
##                nloaded, self.population.generation))
##            # Save the state of the algorithm
##            if not self.saveState():
##                msg = """DeDriver ERROR: failed to save state (state file: %s).
##                Algorithm is paused (lock file: %s)""" % (self.stateFName,
##                                                          self.lockFName)
##                self.writeLog('ERROR -- DeDriver::cafJobComplete()\n'+msg)
##                self.sendMail("DeDriver ERROR", msg)
##                self.pause()
            return False
        else:
            # Not enough time has passed since the submission of the job,
            # so we should not expect to have a complete set of output files.
            return False
    

    def updateStrategy(self):
        """ Update the current DE strategy based upon the fraction of
        the maximum number of generations completed so far.
        """        
        if self.population.generation < self.xi * self.Gmax:
            self.fi = 0
            self.ctrlPop.fi = 0
            self.population.strategy = self.phase1Strategy
            self.population.F = self.phase1F
            if self.population.nDim < self.dimThresh:
                self.population.crossProb = self.phase1LoDimCr
            else:
                self.population.crossProb = self.phase1HiDimCr
        else:
            self.fi = 1
            self.ctrlPop.fi = 1
            self.population.strategy = self.phase2Strategy
            self.population.F = self.phase2F
            if self.population.nDim < self.dimThresh:
                self.population.crossProb = self.phase2LoDimCr
            else:
                self.population.crossProb = self.phase2HiDimCr
            

    def getParentList(self):
        """ Return the list of population member indices which are up
        for replacement in the current CAF segment.
        """
        return self.cafSegmentMap[self.cafSegment]


    def getCtrlParentKeys(self, zeta):
        """ Returns a list of keys for the ctrlPop member indices
        which are up for replacement in the current CAF segment.
        """
        allKeys = self.ctrlPop.getPopFrac(zeta, self.fi)
        segKeyMap = {}
        mptr = 0
        nseg = len(self.cafSegmentMap)
        nCtrlParentsPerSeg = int(math.ceil(float(len(allKeys)) / float(nseg)))
        for iseg in range(1,nseg+1):
            if mptr+nCtrlParentsPerSeg < len(allKeys):
                upper = mptr + nCtrlParentsPerSeg
            else:
                upper = len(allKeys)
            segKeyMap[iseg] = allKeys[mptr:upper]
            mptr += nCtrlParentsPerSeg
            if mptr >= len(allKeys):
                mptr = 0
        return segKeyMap[self.cafSegment]


    def evalFcn(self, trial):
        """ Evaluate the objective function for Member object 'trial'.
        """
        result = [Population._LARGEVAL] * self.population.fDim
        if self.local:
            # DO WE NEED TO HAVE BOTH PYTHON CALLABLE FUNCTIONS
            # AND COMPLETELY EXTERNAL FUNCTIONS?  IF SO, NEED A
            # WAY TO DISTINGUISH BETWEEN THESE TWO MODES.
            #
            # NEED A WAY TO PASS EXTRA ARGS TO self.costFcn!!!
            try:
                result = self.costFcn(trial.x, self.costArgs)
            except NameError, TypeError:
                print 'Please set the cost function.'
                raise TypeError
            except:
                pass
        else:
            # Create the actual command to execute
            # First, convert the parameters to strings
            args = ''
            for arg in trial.x:
                args += ' %s' % arg
            # Now append the argument string to the command supplied by the user
            cmd = self.cafCmd + args
            # Spawn a child process and block while it executes        
            child = popen2.Popen4(cmd)
            # close child's STDIN
            child.tochild.close()
            # close child's STDOUT + STDERR -- don't care about log info
            child.fromchild.close()
            # Block while the child process is running, capture the return value
            retval = child.wait()
            if retval == 0:
                # Read the fucntion values from the expected file            
                # What about child processes which return very
                # quickly with a non-zero value?  Do we generate a different
                # parameter vector and try again? (up to some max # of tries)
                try:
                    file = open(self.fcnEvalFName, 'r')
                except IOError:
                    # Problems opening the file with the function values,
                    # so return Population._LARGEVAL for each function
                    return result
##                 if self.pareto:
##                     # Get a list object with each element corresponding to
##                     # a line from the file (including it's newline char)
##                     tmp = file.readlines()        
##                     for i in range(len(tmp)):
##                         # Trim off the trailing newline character
##                         result[i] = float(tmp[i].rstrip('\n'))
##                 else:
##                     tmp = file.readline()
##                     result = float(tmp.rstrip('\n'))
                # Get a list object with each element corresponding to
                # a line from the file (including it's newline char)
                tmp = file.readlines()        
                for i in range(len(tmp)):
                    # Trim off the trailing newline character
                    result[i] = float(tmp[i].rstrip('\n'))
                file.close()
                # Remove the function eval results file to avoid
                # collisions.
                try:
                    os.remove(self.fcnEvalFName)
                except:
                    pass
        if len(result) < self.population.fDim:
            msg = 'Cost function output length < %s.' % self.population.fDim
            raise TypeError, msg
        return result


    def driveCron(self):
        """ Drive the algorithm while in cron mode.
        """
        if self.cafJobComplete():
            logsep = '='*50
            genmsg = 'DeDriver (cron mode) Generation %s' % (
                self.population.generation+1)
            msg = logsep.center(80) + '\n\n' + genmsg.center(80) + '\n'
            self.writeLog(msg)
            # Load the Member objects of the next generation
            nloaded = self.loadNextGen(int(round(self.cafJobCompltFrac *
                                                 len(self.cafSegmentMap))))
            #########
            if self.population.generation == 1:
                # Just loaded the initialized population, so check the
                # cost functions
                ctr = 0
                for key in self.population.keys():
                    if self.population[key].y == ([Population._LARGEVAL]*
                                                  self.population.fDim):
                        ctr += 1
                if ctr == nloaded:
                    self.pause()
                    msg = 'ERROR: The population initialization failed!'
                    msg += ('       All loaded members had cost values of %s!'
                            % Population._LARGEVAL)
                    self.writeLog(msg)
                    self.sendMail('ERROR: Population Init Failed!', msg)
                    sys.exit(1)
            #########
            self.writeLog('\nLoaded a total of %s Members of generation %s.\n' % (
                nloaded, self.population.generation))
            # Update the DE strategy based upon the current fraction of Gmax
            # generations we have already completed.
            self.updateStrategy()
            if self.sade:
                if self.ctrlPop.generation > 2:
                    # Update the control parm population
                    # (remove dominated members)
                    self.ctrlPop.update()
            # Save the statistics for the current generation
            self.saveStatistics()
            # Backup the current state file -- guaranteed to be the state 
            # file of the previous generation 
            bakfname = self.stateFName + '.bak'
            try:
                shutil.copy(self.stateFName, bakfname)
            except:
                pass
            if self.converged():
                # SEND MAIL stating that we have converged!
                msg = """DeDriver has converged!
                \nMax Delta(stat) for the last %s generations:
                Max Delta(best cost value):             %s
                Max Delta(worst cost value):            %s
                Max Delta(mean cost value):             %s
                Max Delta(standard deviation of cost):  %s
                Max Delta(fractional tolerance):        %s
                Max Delta(chi square):                  %s\n""" % (
                    self.m,
                    self.maxDeltaStats[0],
                    self.maxDeltaStats[1],
                    self.maxDeltaStats[2],
                    self.maxDeltaStats[3],
                    self.maxDeltaStats[4],
                    self.maxDeltaStats[5])
                self.writeLog('DeDriver Converged!\n'+msg)
                self.sendMail('DeDriver Converged!', msg)
                # Set the lock file so that no more CAF jobs
                # are submitted, in case the user forgets to
                # remove the cronjob; also dump the algorithm state
                self.pause()
            else:                
                # save the state prior to job submission
                if not self.saveState():
                    msg = """DeDriver ERROR: failed to save state (state file: %s).
                    Algorithm is paused (lock file: %s)""" % (self.stateFName,
                                                              self.lockFName)
                    self.writeLog('ERROR -- DeDriver.driveCron()\n'+msg)
                    self.sendMail("DeDriver ERROR", msg)
                    self.pause()
                    sys.exit(1)
                # Submit a CAF job
                if not self.submitCafJob():
                    # mail is sent and alg is paused within submitCafJob()
                    sys.exit(1)
            # Save the state of the algorithm after the CAF job has been submitted
            # or we have converged.
            if not self.saveState():
                msg = """DeDriver ERROR: failed to save state (state file: %s).
                Algorithm is paused (lock file: %s)""" % (self.stateFName,
                                                          self.lockFName)
                self.writeLog("ERROR -- DeDriver.driveCron()\n"+msg)
                self.sendMail("DeDriver ERROR", msg)
                self.pause()
                sys.exit(1)
            if not self.shelveState():
                msg = """Failed to shelve state to DBM file.  Execution continues...
                """
                self.writeLog("WARNING -- DeDriver.driveCron()\n"+msg)
            # Update the web-monitor files
            self.updateMonitor()


    def getCtrlMember(self, parentKey, ctrlParentKey):
        """ Get the DE control parameters to use for a given pop
        index.
        """
        # SADE
##         if parentKey % 2 == 0:
        if self.population.generation % 2 == 0:
            # Generate a new control parameter individual
            ctrlChild = self.ctrlPop.getTrialMember(ctrlParentKey)
            while not ctrlChild.isValid or ctrlChild.deltaConstr != 0:
                ctrlChild = self.ctrlPop.getTrialMember(ctrlParentKey)
        else:
            # Select a control parameter individual from the
            # corresponding Pareto front.
            ctrlChild = self.ctrlPop.selCtrlMember(self.zeta, self.fi)
        return ctrlChild


    def select(self, parent, child=None, ctrlParent=None,
               ctrlChild=None, toFile=True):
        """ Between two Members, select which survives to the next population.
        """
        nextGen = self.population.generation + 1
        selMbr = None
        selCtrlMbr = None
        if (isinstance(ctrlParent,SAMember)):
            sade = True
            selCtrlMbr = ctrlParent
        else:
            sade = False
        if not isinstance(child,Member):
            # This must be the initialization generation,
            # so just accept the parent.
            try:
                parent.y = self.evalFcn(parent)
            except:
                pass
            parent.generation = nextGen
            if toFile:
                try:
                    gfile = gzip.open(self.cafOutFName, 'ab')
                    cPickle.dump(parent, gfile)
                    gfile.close()
                except:
                    pass
            else:
                self.population[parent.getKey()] = parent
            return
        if child.isValid:
            if child.deltaConstr == 0:
                # no soft constraint violations, so evaluate the cost
                # function
                try:
                    child.y = self.evalFcn(child)
                    # Need to wait until *after* the function eval is done
                    # to open the output file!
                    if child <= parent:
                        # Save the child ; it will survive to the next
                        # generation of the population.                  
                        # Use '<=' rather than '<' for the selection
                        # test to avoid stagnation of the algorithm.
                        selMbr = child
                    else:
                        # Parent is better, so parent survives to the
                        # next generation.                        
                        selMbr = parent
                except:
                    selMbr = parent
                if sade:
                    if isinstance(ctrlChild, SAMember):
                        try:
                            # Evaluate the DE control parameters used
                            ctrlChild.y = self.ctrlPop.evalCtrl(self.population,
                                                                parent, child, ctrlChild)
                        except:
                            ctrlChild.isValid = False
                        if ctrlChild.isValid:
                            if ctrlChild < ctrlParent:
                                # Place the child into the current ctrl parm
                                # population, replacing the parent
                                selCtrlMbr = ctrlChild
                            elif ctrlChild == ctrlParent:
                                # child and parent are non-dominated w.r.t.
                                # each other
                                if ctrlChild.isUnique(ctrlParent):
                                    # Child is unique w.r.t. parent, so accept both.
                                    ctrlChild.popIndex = max(self.ctrlPop.keys()) + 1
                                    selCtrlMbr = ctrlChild
                                else:
                                    # child is not unique, so select the member
                                    # which has the best perfomance on the current
                                    # objective. NOTE: child should have parent's
                                    # population index at this point.
                                    if ctrlChild.y[self.fi] <= ctrlParent.y[self.fi]:
                                        selCtrlMbr = ctrlChild
                    else:
                        try:
                            ctrlParent.y = self.ctrlPop.evalCtrl(self.population,
                                                                 parent, child, ctrlParent)
                        except:
                            # We have no ctrl child, so we accept the ctrl parent,
                            # even though it's not valid.
                            ctrlParent.isValid = False
            else:
                # 'child' violates one or more soft constraints
                # 'child' will have a cost value of Member._LARGEVAL
                if child.deltaConstr < parent.deltaConstr:
                    # child violates the soft constraints less
                    # than the parent, so child survives.
                    # child.generation = parent.generation + 1
                    selMbr = child
                else:
                    # parent is either viable or more viable (violates
                    # constraints less) than the child, so parent survives
                    # parent.generation += 1
                    selMbr = parent
##                 if sade:
##                     # Since the child was not evaluated, ctrlChild
##                     # cannot be evaluated; ctrlParent survives.
##                     selCtrlMbr = ctrlParent
        else:
            # child is not viable w.r.t the hard constraints, so parent
            # survives.
            selMbr = parent
##             if sade:
##                 # Since the child was not evaluated, ctrlChild
##                 # cannot be evaluated; ctrlParent survives.
##                 selCtrlMbr = ctrlParent
        selMbr.generation = nextGen
        if sade:
            selCtrlMbr.generation = nextGen
        if toFile:
            try:
                gfile = gzip.open(self.cafOutFName, 'ab')
                cPickle.dump(selMbr, gfile)
                if sade and selCtrlMbr is not None:
                    cPickle.dump(selCtrlMbr, gfile)
                gfile.close()
            except:
                pass
        else:
            self.population[selMbr.getKey()] = selMbr
            if sade and selCtrlMbr is not None:
                self.ctrlPop[selCtrlMbr.getKey()] = selCtrlMbr

        
    def driveCaf(self):
        """ Drive the algorithm while in CAF mode.
        """
        ctrlParent = None
        ctrlChild = None
        self.setCafOutFileName()
        nextGen = self.population.generation + 1
        if self.sade:
            zeta = self.zeta
            if self.population.generation < 2:
                # Eval all the ctrlPop members
                zeta = 1.0
##             ctrlParentKeys = self.ctrlPop.getPopFrac(zeta, self.fi)
            ctrlParentKeys = self.getCtrlParentKeys(zeta)
            ctrlParentIndex = -1
        # Loop over the parent Member objects associated
        # with the current CAF segment.
        for parentIndex in self.getParentList():
            parent = copy.deepcopy(self.population[parentIndex])
            if self.sade:
                ctrlParentIndex += 1
                if ctrlParentIndex >= len(ctrlParentKeys):
                    ctrlParentIndex = 0
            if nextGen == 1:
                child = None
            else:
                if self.sade:
                    # Get the DE control parameters
                    ctrlParent = self.ctrlPop[ctrlParentKeys[ctrlParentIndex]]
                    if nextGen >= 2:
                        # The following call is correct in using parent.getKey();
                        # ctrlChild is either generated or sampled, depending
                        # on parent.getKey().
                        ctrlChild = self.getCtrlMember(parent.getKey(),
                                                       ctrlParent.getKey())
                        self.population.F = ctrlChild[0]
                        self.population.crossProb = ctrlChild[1]
                    else:
                        self.population.F = ctrlParent[0]
                        self.population.crossProb = ctrlParent[1]
                child = self.population.getTrialMember(parent.getKey())
            self.select(parent, child, ctrlParent, ctrlChild, True)


    def minimize(self):
        """ A complete minimization driver; for use in local running only!
        Assumes all parameters (including self.costFcn) have been set.
        """
        ctrlParent = None
        ctrlChild = None
        pkeys = self.population.keys()
        pkeys.sort()
        while not self.converged():
            self.updateStrategy()
            nextGen = self.population.generation + 1
            if self.sade:
                zeta = self.zeta
                if self.population.generation < 2:
                    # Eval all the ctrlPop members
                    zeta = 1.0
                ctrlParentKeys = self.ctrlPop.getPopFrac(zeta, self.fi)
                ctrlParentIndex = -1
            for parentIndex in pkeys:
                parent = copy.deepcopy(self.population[parentIndex])
                if self.sade:
                    ctrlParentIndex += 1
                    if ctrlParentIndex >= len(ctrlParentKeys):
                        ctrlParentIndex = 0
                if nextGen == 1:
                    child = None
                else:
                    if self.sade:
                        # Get the DE control parameters
                        ctrlParent = self.ctrlPop[ctrlParentKeys[ctrlParentIndex]]
                        if nextGen >= 2:
                            # The following call is correct in using parent.getKey();
                            # ctrlChild is either generated or sampled, depending
                            # on parent.getKey().
                            ctrlChild = self.getCtrlMember(parent.getKey(),
                                                           ctrlParent.getKey())
                            self.population.F = ctrlChild[0]
                            self.population.crossProb = ctrlChild[1]
                        else:
                            self.population.F = ctrlParent[0]
                            self.population.crossProb = ctrlParent[1]
                    child = self.population.getTrialMember(parent.getKey())
                self.select(parent, child, ctrlParent, ctrlChild, False)
            self.population.setGeneration(nextGen)
            self.saveStatistics()
            if self.sade:
                self.ctrlPop.setGeneration(nextGen)
                if nextGen > 1:
                    # Update the control parm population
                    # (remove dominated members)
                    print 'About to update ctrlPop... np = %s' % len(self.ctrlPop.keys())
                    self.ctrlPop.update()
                    print 'CtrlPop updated.  np = %s' % len(self.ctrlPop.keys())
##             if nextGen % 10 == 0:
            if nextGen > 0:
                logsep = '='*50
                genmsg = 'DeDriver (cron mode) Generation %s' % (
                    self.population.generation)
                msg = logsep.center(80) + '\n\n' + genmsg.center(80) + '\n'
                msg += `self`
                print msg
                print '\nControl Parm Population Summary:'
                print `self.ctrlPop`
                sys.stdout.flush()
##                 self.writeLog(msg)
        # We've converged, so return the best cost value
        return self.population[self.population.ibest].y
    
#--------------------------------------------------------------------------


def monitorMode(parser, options):
    """ Monitor the progress of a DeDriver run.
    """
    try:
        stateFile = gzip.open(options.stateFile, 'rb')
        deDriver = cPickle.load(stateFile)
    except (IOError, gzip.zlib.error, EOFError, cPickle.UnpicklingError):
        # handle the IO error
        # Exit with an error condition
        print 'DeDriver Error: Failed to open DeDriver state file %s' % options.stateFile
        sys.exit(1)
    stateFile.close()
    print '\n'
    print `deDriver`


def cronMode(parser, options):
    """ Operate the DeDriver class in cronjob mode.
    """
    # Check for the presence of a lock file (created by the user);
    # if it exists, do nothing (i.e. do *not* execute any part of the
    # optimization algorithm -- no CAF jobs will be submitted!).
    ########## WHAT ABOUT RESTARTING THE ALG FROM WHERE WE LEFT OFF? ########
    ### Will need to set the caf submit time so that DeDriver::cafJobComplete()
    ### does not set the lock file again
    ### self.cafSubmitTime = time.time()
    lockFName = './dedriver.loc'  # should this be a user option?
    if not os.path.exists(lockFName):
        try:
            # This should probably be split into two separate try statements...
            stateFile = gzip.open(options.stateFile, 'rb')
            deDriver = cPickle.load(stateFile)
        except (IOError, gzip.zlib.error, EOFError, cPickle.UnpicklingError):
            # handle the IO error -- What about error message from Python?
            # set the lock file
            Pause(lockFName)
            # SEND MAIL
            addr = '%s@fnal.gov' % os.environ['USER']
            subj = 'DeDriver StateFile error'
            msg = 'Failed to open the DeDriver state file, %s' % options.stateFile
            msg += """\n\nThe DE algorithm will be paused until the user removes
            the lock file, %s""" % lockFName
            WriteLog(options.logFile, 'ERROR -- DeDriver State file error\n'+msg)
            SendMail(addr, subj, msg)
            sys.exit(1)
        stateFile.close()
        deDriver.driveCron()


def cafMode(parser, options):
    """ Operate the DeDriver class in CAF mode.
    """
    if options.segment < 0:
        parser.error("Invalid CAF segment number.")
    try:
        stateFile = gzip.open(options.stateFile, 'rb')
        deDriver = cPickle.load(stateFile)
    except (IOError, gzip.zlib.error, EOFError, cPickle.UnpicklingError):
        # handle the IO error
        print 'DeDriver (CAF mode) ERROR: Failed to open the state file, %s' % options.stateFile
        sys.exit(1)
    stateFile.close()
    deDriver.cafSegment = options.segment
    deDriver.driveCaf()
            

def setupMode(parser, options):
    """ Operate the DeDriver class in setup mode.
    """
    # Regular expression used for checking for the absence of
    # a required argument to an option (e.g. user inputs
    # '-f -F blah.dat', in which case the argument of
    # the '-f' option is '-F', rather than a valid filename)
    # This 'manual' check is only necessary for options
    # expecting string type arguments; other types will
    # be checked by the optparse module.
    if options.debug:
        print '%s: Entering setup mode...' % sys.argv[0]
    spatt = re.compile(r"^-")
    # CAF Parameters
    if options.debug:
        print '%s: Testing options.cmd ...' % sys.argv[0]
    if not options.cmd or spatt.search(options.cmd):
        parser.error("Please specify a command for function evaluation")
    if options.nDim < 1:
        parser.error("Parameter vector length must be > 0")
    if options.dtol < 0.0:
        parser.error("delta-tol must be >= 0")
    # Check for mal-formed email address
    if options.debug:
        print '%s: Testing options.email ...' % sys.argv[0]
    epatt = re.compile(r"^\S+@(?:\S+|\.)+")
    if not epatt.search(options.email):
        parser.error("Malformed email address: %s" % options.email)
    if options.debug:
        print '%s: Testing options.stateFile ...' % sys.argv[0]
    gzpatt = re.compile(r"\.gz$")
    if not gzpatt.search(options.stateFile):
        print """DeDriver Warning: State file will be gzipped; appending
        '.gz' extension to %s""" % options.stateFile
        options.stateFile += '.gz'
    if options.debug:
        print '%s: Testing options.fcnEvalFile ...' % sys.argv[0]
    if spatt.search(options.fcnEvalFile):
        parser.error("""Please specify a file to use for function value
        passing""")
    if options.debug:
        print '%s: Testing options.group ...' % sys.argv[0]
    if spatt.search(options.group):
        parser.error("Malformed CAF group name %s" % options.group)
    if options.Gmax < 0:
        parser.error("Gmax must be > 0")
    ### Test the hard and soft constraints together
    if options.debug:
        print '%s: Testing options.hardConstr ...' % sys.argv[0]
    relOps = re.compile("[<>=]")
    parmPatt = re.compile("x\[\d\]")
    for constr in (options.hardConstr+options.softConstr):
        if relOps.search(constr):
            parser.error("Constraint %s contains relational operators" % constr)
        if not parmPatt.search(constr):
            parser.error("Constraint %s contains no parameter variables like x[0]"
                         % constr)
    if spatt.search(options.logFile):
        parser.error("Malformed log file %s" % options.logFile)
    if options.m < 1 or options.m > 10:
        parser.error("""Number of previous generations to consider for
        convergence should be in the domain [1, 10]""")
    if options.nEvals < 1:
        parser.error("Must have at least one function eval per job segment")
    # Check the basic format of the output location
    if options.debug:
        print '%s: Testing options.outLoc ...' % sys.argv[0]
        print '\toptions.outLoc = %s' % options.outLoc
    olpatt = re.compile(r"^\S+@(?:\S+|\.)+:(?:/|~)")    
##    if not olpatt.search(options.outLoc):
##        parser.error("Please specify a directory for job output in the form user@host.edu:/path")
    if options.debug:
        print '%s: Testing options.queue ...' % sys.argv[0]
    if options.queue not in DeDriver._stdQueueData.keys():
        parser.error("Invalid CAF queue type %s" % options.queue)
    # Check the format of the initial range list
    # should be of the form 'int:float:float'
    if options.debug:
        print '%s: Testing options.initRange ...' % sys.argv[0]
    irdata = []
    irctr = 0
    for item in options.initRange:
        irlist = item.split(':', 2)
        index = -1
        lo = Member._defLoBnd
        hi = Member._defHiBnd
        if len(irlist) != 3:
            parser.error("Malformed initial range string %s" % item)
        try:
            index = int(irlist[0])
            lo = float(irlist[1])
            hi = float(irlist[2])
        except (ValueError, TypeError):
            parser.error("Malformed initial range string %s" % item)
        else:
            if irctr <= options.nDim:
                irctr += 1
                irdata.append((index, lo, hi))    
    irctr = 0
    gaussParms = []
    for item in options.gaussParms:
        gaussList = item.split(':',2)
        index = -1
        mu = Population._defGaussParms[0]
        sigma = Population._defGaussParms[1]
        if len(gaussList) != 3:
            parser.error("Malformed gaussian parameter string %s" % item)
        try:
            index = int(gaussList[0])
            mu = float(gaussList[0])
            sigma = float(gaussList[0])
        except:
            parser.error("Malformed gaussian parameter string %s" % item)
        else:
            if irctr <= options.nDim:
                irctr += 1
                gaussParms.append((index, mu, sigma))
    if options.popSize < 5:
        parser.error("Population size must be >= 4")
    if options.popSize < options.nDim:
        parser.error("Population size must be >= dimension of parameter vector")
    if options.debug:
        print '%s: Testing options. ...' % sys.argv[0]
    if options.tarfile:
        parser.error("""%s: Monolithic tarball submission is not currently
        supported.  Exiting...""" % sys.argv[0])
##         if options.debug:
##             print '%s: Testing options.tarfile ...' % sys.argv[0]
##         if spatt.search(options.tarfile):
##             parser.error("Please specify a function evaluation tarball.")
##         if not os.path.exists(options.tarfile):
##             parser.error("Tarball for function evaluation %s does not exist!"
##                          % options.tarfile)
    if options.tol < 0.0:
        parser.error("tol must be >= 0")
    # Check the basic format of the function tarball URL
    if options.debug:
        print '%s: Testing options.url ...' % sys.argv[0]
    urlpatt = re.compile(r"(?:^http://)|(?:^ftp://)")
    if not urlpatt.search(options.url):
        parser.error("Malformed objective function tarball URL.")
    if options.url and options.tarfile:
        print """DeDriver Warning: Split tarball framework and monolithic
        tarball framework are mutually exclusive.  Will use split framework."""
        options.tarfile = None
    elif not options.url and not options.tarfile:
        parser.error("""Please specify either a URL for the function eval
        tarball (for split tarball framework) or a filename for monolithic
        tarball framework.""")
    if options.xi < 0 or options.xi > 1.0:
        parser.error("Exploratory strategy fraction xi must be in [0.0, 1.0]")
    if options.cmpltSegFrac < 0.0 or options.cmpltSegFrac > 1.0:
        parser.error("""Fraction of complete segments returned for each job
        should be in [0.0, 1.0].""")
    if spatt.search(options.ph1Strat):
        parser.error("Please specify a valid phase 1 DE strategy.")
    if options.ph1LDCr < 0.0 or options.ph1LDCr > 1.0:
        parser.error("""Phase 1 low-dim crossover probabilty should
        be in [0.0, 1.0].""")
    if options.ph1HDCr < 0.0 or options.ph1HDCr > 1.0:
        parser.error("""Phase 1 high-dim crossover probabilty should
        be in [0.0, 1.0].""")
    if options.ph1F < 0.0 or options.ph1F > 2.0:
        parser.error("""Phase 1 DE scale factor F should be in
        [0.0, 2.0].""")        
    if spatt.search(options.ph2Strat):
        parser.error("Please specify a valid phase 2 DE strategy.")
    if options.ph2LDCr < 0.0 or options.ph2LDCr > 1.0:
        parser.error("""Phase 2 low-dim crossover probabilty should
        be in [0.0, 1.0].""")
    if options.ph2HDCr < 0.0 or options.ph2HDCr > 1.0:
        parser.error("""Phase 2 high-dim crossover probabilty should
        be in [0.0, 1.0].""")
    if options.ph2F < 0.0 or options.ph2F > 2.0:
        parser.error("""Phase 2 DE scale factor F should be in
        [0.0, 2.0].""")        
    if options.debug:
        print '%s: Testing options.dataSrc ...' % sys.argv[0]
    if spatt.search(options.dataSrc):
        parser.error("Please specify a valid CAF data source.")
    if options.debug:
        print '%s: Testing options.convgStrat ...' % sys.argv[0]
    if options.convgStrat not in DeDriver._availConvgStrats:
        parser.error("Invalid convergence strategy %s" % options.convgStrat)
    if options.Gmin < 0:
        parser.error("Gmin must be > 0")
    if options.Gmax < options.Gmin:
        #parser.error("Gmax must be greater than Gmin")
        print "DeDriver Warning: Gmax < Gmin --> Gmax = Gmin"
    if options.dimThresh < 0 or options.dimThresh > 10:
        parser.error("Invalid dimensional threshold for DE crossover parameter.")
    if options.truncFrac < 0 or options.truncFrac > 1.0:
        parser.error("Truncation fraction for chi-square must be in [0,1].")
##    if spatt.search(options.lockFile):
##        parser.error("Invalid lock file %s" % options.lockFile)
    if options.mt < 0 or options.mt > 1.0:
        parser.error("Trigonometric mutation probability must be in [0,1].")
    if spatt.search(options.buildDir):
        parser.error("Please specify a build directory for CAF job submission.")    
    if spatt.search(options.runLocalDir):
        parser.error("Please specify a working directory for local function eval.")
##    if options.monitorLoc == '':
##        print "Warning:  No monitor location specified, so web-based monitoring will not be available."
        #parser.error("Please specify a monitor location which is served to the web")

    ### Create a DeDriver instance
    if options.debug:
        print '%s: Creating the DeDriver instance...' % sys.argv[0]
    fdim = 1
    deDriver = DeDriver(options.nDim, fdim, options.popSize, -1,
                        options.stateFile, options.email)
    deDriver.verbose = options.verbose
    deDriver.debug = options.debug
    deDriver.logFName = options.logFile
    if options.saveInitPopFile != '':
        deDriver.saveInitPop = True
        deDriver.popFile = options.saveInitPopFile
    if options.loadInitPopFile != '':
        deDriver.initFromFile = True
        deDriver.popFile = options.loadInitPopFile
    # Set the hard and soft constraints
    if options.debug:
        print '%s: Setting DeDriver constraints...' % sys.argv[0]
    deDriver.setConstraints(options.softConstr, options.hardConstr)
    # Set phase1 / phase2 DE parms.
    # NOTE: currently, we use the default phase1/2 strategies  
    if options.debug:
        print '%s: Setting DeDriver DE strategy...' % sys.argv[0]
    
    deDriver.setDEStrategy(options.xi, options.dimThresh, options.ph1LDCr,
                           options.ph1HDCr, options.ph1F, options.ph2LDCr,
                           options.ph2HDCr, options.ph2F, options.mt,
                           options.ph1Strat, options.ph2Strat,
                           options.sade, options.zeta)
##                           ('rand-trig',1), ('best',1))
##                           ('rand-trig',1), ('best-trig',1))
##                           ('rand-sort',1), ('best',1))
##                           ('rand-sort',1), ('best',2))
    # Set the convergence strategy parameters
    if options.debug:
        print '%s: Setting DeDriver convergence strategy...' % sys.argv[0]
    deDriver.setConvergenceStrategy(options.convgStrat, options.Gmin,
                                    options.Gmax, options.tol, options.dtol,
                                    options.m, options.truncFrac)
    # Set the CAF parameters
    if options.debug:
        print '%s: Setting DeDriver CAF parameters...' % sys.argv[0]
    local = False
    if options.runLocalDir != '':
        local = True        
    deDriver.setCafParms(options.cmd, options.nEvals, options.url, options.tarfile,
                         options.queue, options.group, options.outLoc,
                         options.dataSrc, options.cmpltSegFrac,
                         local, options.runLocalDir)
    if options.debug:
        print '%s: Initializing DeDriver instance...' % sys.argv[0]
    deDriver.initialize(irdata, gaussParms)
    msg = (('='*50).center(80) + '\n\n' + ('DeDriver Initialized').center(80)
           + '\n\n' + deDriver.printFull() + '\n')
    deDriver.writeLog(msg)
    print 'DeDriver initialized!'
    print '** Please setup a cronjob to run %s in cron mode.\n' % sys.argv[0]
    print `deDriver`
    if options.debug:
        print '%s: Saving the DeDriver state...' % sys.argv[0]
    if not deDriver.saveState():
        print """DeDriver ERROR: failed to save state.  Algorithm is paused
        (lock file is set)."""
        deDriver.pause()
        sys.exit(1)
    else:
        if deDriver.debug:
            print '%s: state saved to %s' % (sys.argv[0], deDriver.stateFName)
            print '%s: Submitting a CAF job...' % sys.argv[0]
        if deDriver.submitCafJob():
            if not deDriver.saveState():
                print """%s: ERROR: failed to save state after CAF job
                submission.  Algorithm is paused (lock file is set).""" % (
                    sys.argv[0])
                deDriver.pause()
                sys.exit(1)
            if not deDriver.shelveState():
                print """%s: WARNING: failed to shelve state to DBM file.
                Execution continues... """ % (sys.argv[0])
        else:
            print """%s: ERROR Submitting CAF job (check the log file).
            Exiting...""" % sys.argv[0]
            sys.exit(1)
        if options.debug:
            print '%s: Updating the monitor info...' % sys.argv[0]
        deDriver.updateMonitor()
        


if __name__ == "__main__":
    # Set the random seed from the system clock or /dev/urandom
    random.seed()
    #farmdict = []
    # determine the operational mode from the command-line args,
    # do the appropriate setup and run the DE algorithm.
    
    usage = "usage: %prog [options] MODE"
    vers = '%prog v' + '%s' % (__version__)
    parser = OptionParser(usage, version=vers)
    # mode is a *required* argument, so process it as a positional arg.
    ##### Use 'metavar="FILE"' to change the metavariable listed in
    ##### the optparse help message.
##    parser.add_option("-a", "--adaptive-search", action="store_true",
##                      default=False, help="Use adaptive search parameters.")
    parser.add_option("-c", "--cmd", action="store", type="string",
                      dest="cmd", help="Add a command for function"+
                      " evaluation; should take the trial parameter vector as"+
                      " positional args and write the cost values to a file"+
                      " (de_fcn.dat).  Do not include the positional args in the"+
                      " command.")
    parser.add_option("-d", "--dim", action="store", type="int",
                      dest="nDim", default=1,
                      help="Parameter vector dimensionality [default: %default].")
    parser.add_option("-D", "--delta-tol", action="store", type="float",
                      dest="dtol", metavar="DTOL", default=0.5E-4,
                      help="Maximum change in convergence measure over"+
                      " the previous M generations [default: %default].")
    parser.add_option("-e", "--email", action="store", type="string",
                      dest="email", default=os.environ['USER']+'@fnal.gov',
                      help="Email address to use for error reports"+
                      " [default: %default].")
    parser.add_option("-f", "--state-file", action="store", type="string",
                      dest="stateFile", metavar="FILE", default="de_state.dat.gz",
                      help="Differential Evloution driver state filename"+
                      " [default: %default].")
    parser.add_option("-F", "--fcn-eval-file", action="store", type="string",
                      dest="fcnEvalFile", metavar="FILE", default="de_fcn.dat",
                      help="Use FILE for communicating function"+
                      " evaluation outputs to the DeDriver instance"+
                      " [default: %default]")
    parser.add_option("-g", "--group", action="store", type="string",
                      dest="group", default="common",
                      help="CAF group [default: %default].")
    parser.add_option("-G", "--gen-max", action="store", type="int",
                      dest="Gmax", metavar="GMAX", default=100,
                      help="Maximum number of generations [default: %default].")
    parser.add_option("-H", "--add-hard-constr", action="append",
                      type="string", dest="hardConstr", metavar="CONSTR",
                      default=[],
                      help="Add a *hard* constraint g(x[]) (e.g."+
                      " \"-[0] + 10\"""), where it is assumed that g(x[]) <= 0."+
                      " Constraints can be linear or nonlinear.  Multiple"+
                      " constraints are possible via passing this option"+
                      " multiple times.  In particular, lower and upper bounds"+
                      " are taken from the list of hard constraints.  Hard"+
                      " constraints will be *strictly* enforced; see the docs"+
                      " for the repair rule.")
    parser.add_option("-i", "--add-gauss-init", action="append",
                      type="string", dest="gaussParms", metavar="GAUSSPARMS",
                      default=[], help="Add a set of gaussian parameters for"+
                      "initialization in the form 'i:mu:sigma', where 'i'"+
                      "is the decision vector index, and mu and sigma are the"+
                      "gaussian mean and width, respectively.")
    parser.add_option("-j", "--segment", action="store", type="int",
                      dest="segment", default=-1,
                      help="CAF job segment number [default: %default].")
    parser.add_option("-l", "--log-file", action="store", type="string",
                      dest="logFile", default="de.log", metavar="LOGFILE",
                      help="Log output to file LOGFILE [default: %default].")
##     parser.add_option("-L", "--local", action="store_true", dest="local",
##                       default=False, help="Run minimization locally [default: %default].")
##     parser.add_option("-L", "--local", action="store", type="string",
##                       dest="runLocalDir", default="./de_local", metavar="DIR"
##                       help="""Run minimization locally, using DIR as the 
##                       working directory for function eval [default: %default].""")
    parser.add_option("-m", "--convg-history", action="store", type="int",
                      dest="m", default=5, help="Number of previous"+
                      " generations to look at for convergence criteria"+
                      " [default: %default].")
    parser.add_option("-n", "--num-evals-per-seg", action="store", type="int",
                      dest="nEvals", metavar="NEVALS", default=1,
                      help="Number of function evaluations per CAF segment"+
                      " [default: %default].")
    parser.add_option("-o", "--out-location", action="store", type="string",
                      dest="outLoc", metavar="OUTLOC",
                      default=os.environ['USER']+'@'+os.uname()[1]+
                      ':'+os.getcwd(),
                      help="CAF Output location [default: %default]")
    parser.add_option("-p", "--pareto", action="store_true", dest="pareto",
                      default=False, help="Use a Pareto algorithm to"+
                      " minimize multiple objective functions [default: %default].")
    parser.add_option("-q", "--queue", action="store", type="string",
                      dest="queue", default="short",
                      help="CAF queue [default: %default].")
    parser.add_option("-R", "--run-local-dir", action="store", type="string",
                      dest="runLocalDir", default="", metavar="DIR",
                      help="Use DIR as the working directory for"+
                      " local function evaluation [default: %default].")
    parser.add_option("-r", "--add-init-range", action="append", type="string",
                      dest="initRange", metavar="RANGESTR", default=[],
                      help="Add an initialization range"+
                      " for parameter i of the form, 'i:lo:hi', where 'lo' and"+
                      " 'hi' are the low and high bounds (floats), respectively."+
                      " These bounds will be used solely for initialization.")
    parser.add_option("-s", "--pop-size", action="store", type="int",
                      dest="popSize", metavar="NP", default=5,
                      help="Population size (>= 5) [default: %default].")
    parser.add_option("-S", "--add-soft-constr", action="append",
                      type="string", dest="softConstr", metavar="CONSTR",
                      default=[],
                      help="Add a *soft* constraint g(x[]) (e.g. \"-x[0] +"+
                      " 1.0E-1 * x[1]\"), where it is assumed that g(x[]) <= 0."+
                      " Constraints can be linear or nonlinear.  Multiple"+
                      " constraints are possible via passing this option"+
                      " multiple times.  Soft constraints may be violated, but"+
                      " these trial solutions will be de-weighted; in particular,"+
                      " non-viable trial solutions will not participate in the"+
                      " convergence test.  See the docs for the exact algorithm.")
    parser.add_option("-t", "--tarfile", action="store", type="string",
                      dest="tarfile", metavar="FILE",
                      help="The tarball to send with the job (monolithic CAF"+
                      " job submission mode)")
    parser.add_option("-T", "--tol", action="store", type="float",
                      dest="tol", default=1.0E-3,
                      help="Fractional tolerance for determining convergence"+
                      " [default: %default].")
    parser.add_option("-u","--tarball-url", action="store", type="string",
                      dest="url", metavar="URL",
                      help="URL for the objective function's"+
                      " tarball (split tarball CAF job submission mode).")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", default=False, help="Be verbose"+
                      " [default: %default].")
    parser.add_option("-x", "--explr-frac", action="store", type="float",
                      dest="xi", metavar="XI", default=0.75,
                      help="Fraction xi of Gmax"+
                      " generations which use an exploratory strategy; (1-xi)"+
                      " *Gmax generations use Elitist strategy [default: %default].")
    parser.add_option("", "--complt-seg-frac", action="store", type="float",
                      dest="cmpltSegFrac", metavar="FRAC", default=0.2,
                      help="Fraction of complete job segments (for each"+
                      " generation) required for algorithm to continue after"+
                      " elapsed job time exceeds queue real-time limit"+
                      " [default: %default].")
    parser.add_option("", "--phase1DEStrategy", action="store", type="string",
                      dest="ph1Strat", metavar="STRAT", default="rand-trig1",
                      help="Set the phase 1 DE strategy to STRAT [default = %default].")
    parser.add_option("", "--phase1LoDimCr", action="store", type="float",
                      dest="ph1LDCr", metavar="CR", default=0.5,
                      help="Phase 1 crossover"+
                      " probability for a low-dim parameter space"+
                      " [default: %default].")
    parser.add_option("", "--phase1HiDimCr", action="store", type="float",
                      dest="ph1HDCr", metavar="CR", default=0.7,
                      help="Phase 1 crossover"+
                      " probability for a high-dim parameter space"+
                      " [default: %default].")
    parser.add_option("", "--phase1F", action="store", type="float",
                      dest="ph1F", metavar="F", default=0.8,
                      help="Phase 1 scale factor F"+
                      " [default: %default].")
    parser.add_option("", "--phase2DEStrategy", action="store", type="string",
                      dest="ph2Strat", metavar="STRAT", default="best1",
                      help="Set the phase 2 DE strategy to STRAT [default = %default].")
    parser.add_option("", "--phase2LoDimCr", action="store", type="float",
                      dest="ph2LDCr", metavar="CR", default=0.0,
                      help="Phase 2 crossover"+
                      " probability for a low-dim parameter space"+
                      " [default: %default].")
    parser.add_option("", "--phase2HiDimCr", action="store", type="float",
                      dest="ph2HDCr", metavar="CR", default=0.1,
                      help="Phase 2 crossover"+
                      " probability for a high-dim parameter space"+
                      " [default: %default].")
    parser.add_option("", "--phase2F", action="store", type="float",
                      dest="ph2F", metavar="F", default=0.5,
                      help="Phase 2 scale factor F"+
                      " [default: %default].")
    parser.add_option("","--sade", action="store_true", dest="sade",
                      default=False, help="Self-Adaptive DE: adapt the"+
                      "DE control parameters (F, Cr) during the run.")
    parser.add_option("","--zeta", action="store", type="float", dest="zeta",
                      default=0.5, help="Percentage of the Pareto front to"+
                      "use when selecting a DE control parameter vector (F, Cr)."+
                      "  SADE parameter.")
    parser.add_option("", "--dataSrc", action="store", type="string",
                      dest="dataSrc", default="None",
                      help="Data source for the CAF [default: %default].")
    parser.add_option("", "--convg-strategy", action="store", type="string",
                      dest="convgStrat", metavar="STRAT", default="fractol",
                      help="Set the convergence strategy for the DE algorithm"+
                      " (either 'fractol' or 'chisq') [default: %default].")
    parser.add_option("", "--Gmin", action="store", type="int",
                      dest="Gmin", default=10,
                      help="Set the minimum number of generations to complete"+
                      " [default: %default].")
    parser.add_option("", "--dim-thresh", action="store", type="int",
                      dest="dimThresh", metavar="THRESH", default=6,
                      help="Set the dimensional threshold for the DE crossover"+
                      " parameter [default: %default].")
    parser.add_option("", "--trunc-frac", action="store", type="float",
                      dest="truncFrac", metavar="FRAC", default=0.15,
                      help="Set the truncation fraction for chi-square"+
                      " computation to FRAC [default: %default].")
    parser.add_option("", "--trig-prob", action="store", type="float",
                      dest="mt", metavar="MT", default=0.05,
                      help="Set the probability for trigonometric mutation"+
                      "to MT; the probability for standard DE mutation is then"+
                      "(1-MT).  NOTE: This option is only applicable to the"+
                      "'rand-trig' and 'best-trig' strategies.")
    parser.add_option("", "--save-init-pop-file", action="store", type="string",
                      dest="saveInitPopFile", metavar="FILE", default="",
                      help="Save the initial population of to file FILE."+
                      "  Should have a '.gz' extension.")
    parser.add_option("", "--load-init-pop-file", action="store", type="string",
                      dest="loadInitPopFile", metavar="FILE", default="",
                      help="Load the initial population from file FILE."+
                      "  Should have a '.gz' extension.")
    # DO NOT SET THE LOCK FILE NAME MANUALLY! THERE'S NO WAY FOR THE cronMode()
    # FUNCTION TO KNOW WHAT IT IS!
##    parser.add_option("", "--lock-file", action="store", type="string",
##                      dest="lockFile", metavar="FILE", default="dedriver.loc",
##                      help="Use LOCKFILE as the lock file for pausing"+
##                      " the algorithm [default: %default].")
    parser.add_option("", "--build-dir", action="store", type="string",
                      dest="buildDir", metavar="DIR", default="./build_de",
                      help="Use BUILDDIR as the tarball build directory"+
                      " for CAF job submission [default: %default].")
    parser.add_option("","--debug", action="store_true", dest="debug",
                      default=False, help="Debug flag.")
##    parser.add_option("", "--monitor-url", action="store", type="string",
##                      dest="monUrl", help="""Set the base URL (node, port, path) of the monitor.""")
##    parser.add_option("-M", "--monitor-location", action="store", type="string",
##                      dest="monitorLoc",
##                      help="""Location of intermediate results for web monitoring (e.g. user@node.edu:/path""")
    
    (options, args) = parser.parse_args()

    # Error handling for options
    availModes = ('setup', 'cron', 'caf', 'monitor')
    sep = ', '
    if len(args) != 1:
        parser.error("Invalid mode; available modes: %s" % sep.join(availModes))

    if args[0] == availModes[0]:
        setupMode(parser, options)
    elif args[0] == availModes[1]:
        cronMode(parser, options)
    elif args[0] == availModes[2]:
        cafMode(parser, options)
    elif args[0] == availModes[3]:
        monitorMode(parser, options)
    else:
        parser.error("Invalid mode; available modes: %s" % sep.join(availModes))
    



