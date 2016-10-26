from gurobipy import *
import random
import numpy 
import math

numscen = 50  # n
numbatch = 15  # M
numeval = 500 # N


### define some data for this simple problem ###
rhs = [7.0, 4.0]
penalty = 5.0

candx1 = 0.0
candx2 = 0.0

### function for building and solving the extensive form for a given sample omega1, omega2
### returns the objective value and the solution
def solveEF(omega1, omega2):

    ef = Model("extensiveform")
    ef.params.logtoconsole=0
    x1 = ef.addVar(obj=1.0)
    x2 = ef.addVar(obj=1.0)
    y1 = {}
    y2 = {}
    for j in range(numscen):
        y1[j] = ef.addVar(obj=penalty/float(numscen))
        y2[j] = ef.addVar(obj=penalty/float(numscen))

    ef.update()

    for j in range(numscen):
        ef.addConstr(omega1[j]*x1 + x2 + y1[j] >= rhs[0])
        ef.addConstr(omega2[j]*x1 + x2 + y2[j] >= rhs[1])

    ef.modelSense = GRB.MINIMIZE
    ef.update()
    ef.optimize()
    return [ef.objval, x1.x, x2.x]

### function for building and solving the second stage problem for 
#### a given value omega1,omega2 and first-stage solution x1,x2
### returns the objective value
def solveSSP(x1val,x2val,omega1,omega2):
    y1val = 0.0
    y2val = 0.0
    y1val = max(rhs[0] - omega1*x1val- x2val,0.0)
    y2val = max(rhs[1] - omega2*x1val- x2val,0.0)
    return x1val+x2val+penalty*(y1val+y2val)

objvals = numpy.zeros(numbatch)

### First separately estimate lower and upper bounds #### 
for k in range(numbatch):

    ### set up and solve extensive form for a sample average approximation problem
    omega1 = [random.uniform(1.0,4.0) for j in range(numscen)]  # creates numscen scenarios
    omega2 = [random.uniform(0.333333,0.666667) for j in range(numscen)]  
    [objvals[k],candx1,candx2] = solveEF(omega1,omega2) 
    

print 'mean SAA objval = ', numpy.mean(objvals)
print 'stdev of SAA objval = ', numpy.std(objvals)
lbmean = numpy.mean(objvals)
lbwidth = numpy.std(objvals)/math.sqrt(numbatch)*2.145

### here, candx1 and candx2 are the solution obtained in last batch above
### more generally, might want to choose "most promising" solution
print 'candidate solution = (', candx1, ',', candx2, ')'
evalvals = numpy.zeros(numeval)
omega1 = [random.uniform(1.0,4.0) for j in range(numeval)]
omega2 = [random.uniform(0.333333,0.666667) for j in range(numeval)]
for k in range(numeval):
    evalvals[k] = solveSSP(candx1,candx2,omega1[k],omega2[k])
  
print 'mean solution objval = ', numpy.mean(evalvals)
print 'stdev of solution objval = ', numpy.std(evalvals)
ubmean = numpy.mean(evalvals)
ubwidth = numpy.std(evalvals)/math.sqrt(numeval)*2.145

print 'ci on lower bound = [', lbmean-lbwidth, ',', lbmean+lbwidth, ']'
print 'ci on upper bound = [', ubmean-ubwidth, ',', ubmean+ubwidth, ']'


#### Now let's try the method which directly estimates optimality gap

omega1 = [random.uniform(1.0,4.0) for j in range(numscen)]  # creates numscen scenarios
omega2 = [random.uniform(0.333333,0.666667) for j in range(numscen)]  
[objval,candx1,candx2] = solveEF(omega1,omega2)
print 'candidate solution = (', candx1, ',', candx2, ')'

gapvals = numpy.zeros(numbatch)
for k in range(numbatch):

    ### set up and solve extensive form for a sample average approximation problem
    omega1 = [random.uniform(1.0,4.0) for j in range(numscen)]  # creates numscen scenarios
    omega2 = [random.uniform(0.333333,0.666667) for j in range(numscen)]  
    [sampleopt,curx1,curx2] = solveEF(omega1,omega2) 

    ### now evaluate the candidate solution using the same scenarios
    curvals = numpy.zeros(numscen)
    for j in range(numscen):
        curvals[j] = solveSSP(candx1,candx2,omega1[j],omega2[j])

    gapvals[k] = numpy.mean(curvals) - sampleopt 

print 'mean gap estimate = ', numpy.mean(gapvals)
print '95% c.i. on gap = [0, ', numpy.mean(gapvals) + numpy.std(gapvals)/math.sqrt(numbatch)*1.761, ']'

