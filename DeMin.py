#!/usr/bin/env python
#
# File: DeMin.py
# Author: Jason Galyardt
# Description:
# Differential Evolution minimizer (direct search evolutionary algorithm)
import os, sys, array, math, cPickle, copy, random, re


class Member:
    """Population member; holds one trial vector + (tuple of) fcn eval(s)"""
    _LARGEVAL = 1.0E20
    _defLoBnd = 0.0
    _defHiBnd = 1.0E3
    _constrEvalMax = 500
    _constrCallMax = 10

    def __init__(self, dim=1, gen=-1, ipop=-1):
        self.nDim = dim
        self.generation = gen     # Negative gen --> trial vect
        self.popIndex = ipop      # index in the population (< 0 --> trial)
        self.x = [0.0] * self.nDim
        self.y = Member._LARGEVAL
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
        if self.y > other.y:
            return 1
        elif self.y == other.y:
            return 0
        else:
            return -1


    def setDim(self, dim):
        if type(dim) != type(1):
            raise TypeError, 'Member.nDim must be an int!'
        if self.nDim < dim:
            for i in range(dim - self.nDim):
                self.x.append(0.0)
            y = Member._LARGEVAL
        elif self.nDim > dim:
            for i in range(self.nDim - dim):
                self.x.pop()
            y = Member._LARGEVAL
        self.nDim = dim


    def makeTrial(self):
        # Make a copy of self, but tag it as a trial;
        # return the copy
        trialMember = copy.deepcopy(self)
        trialMember.generation = -1
        trialMember.popIndex = -1
        trialMember.y = Member._LARGEVAL
        return trialMember

    
    def __add__(self, other):
        """ Overload the '+' operator; add the parameter values for
        each index rather than concatenating the underlying list objects. """
        if not isinstance(other, Member):
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
        if not isinstance(other, Member):
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
        if not isinstance(other, Member):
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
        if not isinstance(other, Member):
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
        if not isinstance(mutant, Member):
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
            while delta > 0.0:
                if evalCtr < Member._constrEvalMax:
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
            if callCtr < Member._constrCallMax:
                if evalCtr < Member._constrEvalMax:
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
    _strategies = (('rand',1), ('best',1), ('best',2), ('rand-sort',1),
                   ('rand-trig',1), ('best-trig',1))

    def __init__(self, dim=1, popSize=4, gen=-1, prob=0.1, f=0.5, mt=0.05):
        self.nDim = dim
        self.generation = gen  ## May want to take 'gen' out of constructor arglist
        self.np = popSize
        self.crossProb = prob
        self.F = f
        self.mt = mt       ## only used for rand-trig, best-trig strategies
        self.ibest = ''
        self.iworst = ''
        self.strategy = Population._strategies[0]
        self.hardConstraints = ()
        self.softConstraints = ()
        self.initRange = [(Population._defInitRange[0],
                           Population._defInitRange[1])]*self.nDim
        self.memberColl = {}        
        for i in range(self.np):
            tmp = Member(self.nDim, self.generation, i)
            tmpkey = tmp.getKey()
            for j in range(tmp.nDim):
                tmp[j] = 0.0
            self.memberColl[tmpkey] = tmp
            if i == 0:
                self.ibest = tmpkey
                self.iworst = tmpkey


    def __getitem__(self, key):
        return self.memberColl[key]


    def __setitem__(self, key, value):
        if not isinstance(value, Member):
            raise TypeError, (
                  'Objects of type %s cannot be stored by Population instances.' 
                  % (value.__class__))
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
        indexpatt = re.compile("\[([0-9]+)\]")
        # This pattern may need a bit of work -- it can't handle the
        # unary minus operator (e.g. '-x[0] + 1' --> leading '-' will
        # cause the pattern to not match) ; for now, just specify a unary minus
        # like so, '-1 * x[0] + 1'
        bndpatt = re.compile("(?:([-+]?)((?:\d+(?:\.\d*)?|\d*\.\d+)(?:[eE][-+]?\d+)?)\s*\*\s*)*x\[([0-9]+)\]\s*(?:([-+])\s*((?:\d+(?:\.\d*)?|\d*\.\d+)(?:[eE][-+]?\d+)?))?")
        # Loop through hard constraints, look for parameter bounds
        for constr in self.hardConstraints:
            j = -1
            strlist = indexpatt.findall(constr)
            loBnd = Member._defLoBnd
            hiBnd = Member._defHiBnd
            if len(strlist) == 1:
                # This constraint only involves one parameter,
                # so it is a boundary candidate.
                j = int(strlist[0])
                if j < self.nDim and j > 0:
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
                            hiBnd = -1.0 * s1 * c1 / c0
                        elif mgroups[0] == '-':
                            # s0 < 0
                            # we have a lower bound
                            loBnd = s1 * c1 / c0
                    for key in self.keys():
                        self[key].setParmBounds(j, loBnd, hiBnd)


    def initialize(self, costFcn, args=None):
        """ Initialize the population with random trial vectors;
        should only be done for generation -1!"""
        # Set the parameter upper and lower bounds from the hard constraints
        # given by the user.
        self.setPopParmBounds()
        for key in self.keys():
            # Initialize the Member objects' parameter vectors
            for i in range(self[key].nDim):
                self[key][i] = random.uniform(self.initRange[i][0], self.initRange[i][1])
            # Repair the hard and soft constraints
            self[key].repairHardConstr(self.hardConstraints)
            self[key].repairSoftConstr(self.softConstraints)
            # Evaluate the initial population's cost functions
            try:
                if args is None:
                    self[key].y = costFcn(self[key].x)
                else:
                    self[key].y = costFcn(self[key].x, args)
            except:
                self[key].y = Population._LARGEVAL
        self.setGeneration(1)


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
            self.iworst = random.choice(orderedKeys)
        bestCost = self[self.ibest].y
        worstCost = self[self.iworst].y
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
                sum += self[key].y
                sumsq += self[key].y**2
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
            for key in truncKeys:            
                # We've already checked every Member key in truncKeys[] for
                # viability and discarded the unviable ones
                sumTrunc += self[key].y
                sumsqTrunc += self[key].y**2
            muTrunc = sumTrunc / float(npTrunc)
            musqTrunc = sumsqTrunc / float(npTrunc)
            diffTrunc = musqTrunc - muTrunc**2
            if diffTrunc > 0:
                varTrunc = npTrunc * diffTrunc / float(npTrunc - 1)
            chisq = 0.0
            for key in truncKeys[1:]:
                chisq += (self[key].y - bestCost)**2
            chisq /= (musqTrunc)
            # 
            range = self[truncKeys[-1]].y - bestCost
            fracTol = 2.0 * abs(range) / (abs(self[truncKeys[-1]].y) + abs(bestCost) +
                                          Population._TINY)
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
            rndKeys = self.getRndMembers(2, parentKey, self.ibest)
            if random.random() <= self.mt:
                # 'mutant' is biased toward the region of lower function
                # value in the available parameter space (defined by
                # the 'rndKeys' list)
                pprime = (abs(self[rndKeys[0]].y) + abs(self[rndKeys[1]].y) +
                          abs(self[self.ibest].y))
                p = [abs(self[self.ibest].y) / pprime,
                     abs(self[rndKeys[0]].y) / pprime,
                     abs(self[rndKeys[1]].y) / pprime]
                mutant = ((self[self.ibest] + self[rndKeys[0]] +
                           self[rndKeys[1]]) / 3.0 +
                          (p[1] - p[0]) * (self[self.ibest] - self[rndKeys[0]]) +
                          (p[2] - p[1]) * (self[rndKeys[0]] - self[rndKeys[1]]) +
                          (p[0] - p[2]) * (self[rndKeys[1]] - self[self.ibest]))
            else:                
                mutant = self[self.ibest] + self.F * (self[rndKeys[0]] -
                                                      self[rndKeys[1]])
        elif self.strategy == ('rand-sort', 1):
            rndKeys = list(self.getRndMembers(3, parentKey))
            rndKeys.sort(lambda a,b: Member.__cmp__(self[a], self[b]))
            shuffledKeys = rndKeys[1:]
            random.shuffle(shuffledKeys)
            mutant = self[rndKeys[0]] + self.F * (self[shuffledKeys[0]] -
                                                  self[shuffledKeys[1]])
        elif self.strategy == ('rand', 1):
            rndKeys = self.getRndMembers(3, parentKey)
            mutant = self[rndKeys[0]] + self.F * (self[rndKeys[1]] -
                                                  self[rndKeys[2]])
        elif self.strategy == ('rand-trig',1):
            rndKeys = self.getRndMembers(3, parentKey)
            if random.random() <= self.mt:
                # 'mutant' is biased toward the region of lower function
                # value in the available parameter space (defined by
                # the 'rndKeys' list)
                pprime = (abs(self[rndKeys[0]].y) + abs(self[rndKeys[1]].y) +
                          abs(self[rndKeys[2]].y))
                p = [abs(self[rndKeys[0]].y) / pprime,
                     abs(self[rndKeys[1]].y) / pprime,
                     abs(self[rndKeys[2]].y) / pprime]
                mutant = ((self[rndKeys[0]] + self[rndKeys[1]] +
                           self[rndKeys[2]]) / 3.0 +
                          (p[1] - p[0]) * (self[rndKeys[0]] - self[rndKeys[1]]) +
                          (p[2] - p[1]) * (self[rndKeys[1]] - self[rndKeys[2]]) +
                          (p[0] - p[2]) * (self[rndKeys[2]] - self[rndKeys[0]]))
            else:                
                # use standard DE/rand/1 strategy
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


    def getTrialMember(self, parentKey):
        """ Generate a new trial vector to possibly replace the ith member
        from members of the current population."""
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
        sigma = 0.1
##         scale = random.gauss(mu,sigma)
        # Breit-Wigner / Cauchy pdf
        scale = mu + sigma * math.tan(math.pi * (random.random() - 0.5))
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


    def update(self):
        """ Update the population; remove all dominated individuals.
        """
        if self.generation % 2 == 0:
            nonDomKeys = self.nonDomSet()
##             while len(nonDomKeys) < 4:
            minkeys = max(int(math.ceil(0.10*self.np)), 4)
            maxtrials = len(self.memberColl)
            ntrials = 0
            while (len(nonDomKeys) < minkeys and
                   ntrials < maxtrials):
                # need to make sure we have enough members to
                # perform a mutation for the control parms.
                tmp = self.getRndMembers(1,*nonDomKeys)
                uniq = True
                for key in nonDomKeys:
                    if not self.memberColl[tmp[0]].isUnique(self.memberColl[key]):
                        uniq = False
                if uniq:
                    nonDomKeys.extend(tmp)
                ntrials += 1
            for key in self.keys():
                if key not in nonDomKeys:
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


#---------------------------------------------------------------------------------


class DeMin:
    """ Driver class for the serial Differential Evolution algorithm.
    USAGE:
    def myCost(vect):
        # example cost function
        return sum(vect)
    ir = [(0, 0.6, 0.8), (1, 0.4, 0.6)] # for dim=2
    de = DeMin(dim, np, -1, ir, myCost)
    de.setConstraints(softConstr, hardConstr)
    de.setDEStrategy(xi, dimThresh, ph1LDCr, ph1HDCr, ph1F,
                     ph2LDCr, ph2HDCr, ph2F, mt, ph1StratStr, ph2StratStr)
    
    de.setConvergenceStrategy('fractol',gmin,gmax,tol,dtol,m,truncFrac)
    de.initialize()
    bestCost = de.minimize()
    """
    # Dictionary of available farms and their crucial info:
    # entries are like this: (group, available queues)
    # Standard queues and their limits (for now just realtime limit, in seconds)
    _availConvgStrats = ('chisq', 'fractol')
    #
    def __init__(self, costFcn, dim=1, np=4, gen=-1, ir=[], args=None):
        """ Constructor for the DeMin class.
        """
        self.costFcn = costFcn
        self.args = args
        self.sade = False
        try:
            if args is None:
                costFcn([0.0]*dim)
            else:
                costFcn([0.0]*dim, args)
        except:
            self.args = None
            self.costFcn = DeMin._costFcn
        # Set the initial range for each parameter
        self.ir = ir
        # Create the population(s)
        self.population = Population(dim, np, gen)
        # Create the control parameter population
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
        # Differential Evolution Parameters
        self.pareto = False
        #self.xi = 0.6 # phase 1 strategy used for first (xi*Gmax) generations
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
        self.phase2LoDimCr = 0.1
        self.phase2HiDimCr = 0.2
        self.phase2F = 0.6
        #self.phase2Strategy = ('best', 1)             # 'DE/best/1/bin'
        self.phase2Strategy = ('best', 2)             # 'DE/best/1/bin'
        # number of previous generations to look at for convergence
        self.m = 5
        # List to hold the statistics of the previous 'm' generations
        # Use integer keys, value should be a tuple:
        # (bestCost, worstCost, meanCost, stdDev, fracTol, chisq, ndf)
        # For multi-objective optimization, each tuple entry should
        # be a tuple whose elements correspond to the various
        # functions.
        #
        # CODEMOD
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
        self.convgStrategy = DeMin._availConvgStrats[1]
        # Use 'tol' to define the tolerance for the convergence strategy,
        # use 'deltaTol' to define the tolerance for the time-stability of the
        # convergence.
        self.tol = 1.0E-3
        self.deltaTol = 1.0E-4
        self.truncFrac = 0.15
        # logistical files
        self.verbose = False
        self.debug = False
        

    def __repr__(self):
        """ Overload __repr__() so that `deDriver` will print the algorithm
        state to STDOUT.
        """
        stats = self.population.getStats(self.truncFrac)
        outStr = ("""DeMin state:
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
        """ % (self.population.generation,
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
                   self.population.getCovMatRepr(stats[8])))
        return outStr


    def setDEStrategy(self, xi=0.7, dimThresh=6,
                      ph1LDCr=0.5, ph1HDCr=0.7, ph1F=0.8,
                      ph2LDCr=0.1, ph2HDCr=0.2, ph2F=0.6, mt=0.05,
                      ph1StratStr='rand-trig1', ph2StratStr='best-trig1',
                      sade=False, zeta=0.5):
        """ Set the parameters for the DE strategy.
        """
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
        if strat not in DeMin._availConvgStrats:
            self.convgStrategy = DeMin._availConvgStrats[1]
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
            # reaching GMax generations.
            self.m = self.GMax
        if truncFrac > 0.0 and truncFrac <= 1.0:
            self.truncFrac = truncFrac
        #
        # CODEMOD
        self.prevMGenStats = [(Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL,
                               Population._LARGEVAL)] * self.m
        self.maxDeltaStats = [Population._LARGEVAL] * len(self.prevMGenStats[0])


    def initialize(self):
        """ Initialize the Population.
        """
        # Expect initRange to be a list of tuples, each tuple
        # of the form (index, lobnd, hibnd); assume that
        # the tuple values are valid.
        # NOTE: if len(ir) < self.population.nDim, the remainder
        # of the parameters will take their initial ranges from
        # the default range, defined in Population.__init__()
        #if len(self.popFile) == 0:        
        for item in self.ir:
            self.population.initRange[item[0]] = tuple(item[1:])
        self.population.initialize(self.costFcn, self.args)
        self.ctrlPop.initialize()
        

    def saveStatistics(self):
        """ Save the stats of the current population.
        """
        # Find the best and worst population members
        # Best and worst members are found in population.getStats(),
        # which is called below -- should happen for every generation.
        # Shift the current values back one index
        del self.prevMGenStats[0:1]
        #self.prevMGenStats.append(self.population.getStats(self.truncFrac)[:7])
        #
        # CODEMOD
        self.prevMGenStats.append(self.population.getStats(self.truncFrac)[:len(self.prevMGenStats[0])])
        # Loop over the statistics defined in Population::getStats()
        for j in range(len(self.prevMGenStats[0])):
            maxDelta = -Population._LARGEVAL
            delta = -Population._LARGEVAL
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


    def minimize(self):
        """ The minimization driver.
        """
        ctrlParent = None
        ctrlChild = None
        self.initialize()
        while not self.converged():
            self.updateStrategy()
            nextGenPop = copy.deepcopy(self.population)
            nextGen = self.population.generation + 1
            nextGenPop.setGeneration(nextGen)
            if self.sade:
                zeta = self.zeta
                if self.population.generation < 2:
                    # Eval all the ctrlPop members
                    zeta = 1.0
                ctrlParentKeys = self.ctrlPop.getPopFrac(zeta, self.fi)
                ctrlParentIndex = -1
            for parentIndex in self.population.keys():
                parent = copy.deepcopy(self.population[parentIndex])
                # child = self.population.getTrialMember(parent.getKey())
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
##             nextGenPop.setGeneration(nextGen)
            # self.population = nextGenPop
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
        # We've converged, so return the best cost value
        return self.population[self.population.ibest].y



## def testDeMin():
##     """ Test the DeMin class minimizer with a 3-D sphere problem.
##     DeMin USAGE:
##     def myCost(vect):
##         # example cost function
##         return sum(vect)
##     ir = [(0, 0.6, 0.8), (1, 0.4, 0.6)] # for dim=2
##     de = DeMin(dim, np, -1, ir, myCost)
##     de.setConstraints(softConstr, hardConstr)
##     de.setDEStrategy(xi, dimThresh, ph1LDCr, ph1HDCr, ph1F,
##                      ph2LDCr, ph2HDCr, ph2F, mt, ph1StratStr, ph2StratStr)
    
##     de.setConvergenceStrategy('fractol',gmin,gmax,tol,dtol,m,truncFrac)
##     de.initialize()
##     bestCost = de.minimize()
##     """
##     # Setup the cost function as a 3-D sphere
##     def SphereFcn(vect):
##         S = sphere(3)
##         if not S.setVect(vect):
##             return sphere._LARGEVAL
##         S.calc()
##         return S.f
##     # Set the initial parameters for the DeMin object
##     dim = 3
##     np = 15
## ##     ir = [(0, 0.6, 0.8), (1, 0.4, 0.6), (2, 0.3, 0.5),
## ##           (3, 0.2, 0.4), (4, 0.1, 0.3)]
##     ir = [(0, -10.0, 10.0), (1, -10.0, 10.0), (2, -10.0, 10.0)]
##     softConstr = []
##     # Constrain the parameters to lie on [-10.0, 10.0]
##     hardConstr = ['-x[0]-10.0', 'x[0]-10.0',
##                   '-x[1]-10.0', 'x[1]-10.0',
##                   '-x[2]-10.0', 'x[2]-10.0']
##     xi = 0.7
##     gmin = 5
##     gmax = 300
##     tol = 1e-4
##     dtol = 5e-3
##     m = 5
##     truncFrac = 0.4
##     de = DeMin(SphereFcn, dim, np, -1, ir)
##     de.setConstraints(softConstr, hardConstr)
##     de.setDEStrategy(xi)
##     de.setConvergenceStrategy('fractol', gmin, gmax, tol, dtol, m, truncFrac)
##     de.initialize()
##     minCost = de.minimize()
##     print 'prevMGenStats:'
##     for i in range(len(de.prevMGenStats)):
##         print
##         for j in range(len(de.prevMGenStats[0])):
##             print '%.6E ' % de.prevMGenStats[i][j],
##     print 'maxDeltaFracTol = %.6E' % de.maxDeltaStats[4]
##     print 
##     print `de`

