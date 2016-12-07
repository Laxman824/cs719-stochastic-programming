import cPickle
from gurobipy import *
import time
import numpy
import math

### Read data from file you choose: commont/uncomment to choose the different files
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('nd8-4-8.pdat','r')

Fset = cPickle.load(dfile)  # set of facilities (list of strings)
Hset = cPickle.load(dfile)  # set of warehouses (list of strings)
Cset = cPickle.load(dfile)  # set of customers (list of strings)
arcExpCost = cPickle.load(dfile)  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
facCap = cPickle.load(dfile)   # facility capacities (dictionary mapping F to floats) 
curArcCap = cPickle.load(dfile)  # current arc capacities (dictionary mapping (i,j) to floats, where either 
                                 # i is facility, j is warehouse, or i is warehouse and j is customer
unmetCost = cPickle.load(dfile)  # penalty for unment customer demands (dicationary mapping C to floats)
demmean = cPickle.load(dfile)  #  mean demand (conditional that it is positive)
demstdev = cPickle.load(dfile)   # stdev of demand (conditional that it is positive)
probnodem = cPickle.load(dfile)  # probability each demand is zero
dfile.close() 


### Define sets of arcs (used as keys to dictionaries)
FHArcs = [(i,j) for i in Fset for j in Hset]  ## arcs from facilities to warehouses
HCArcs = [(i,j) for i in Hset for j in Cset]   ## arcs from warehouses to customers
AllArcs = FHArcs + HCArcs
print AllArcs

### Make them Gurobi tuplelists
FHArcs = tuplelist(FHArcs)
HCArcs = tuplelist(HCArcs)
AllArcs = tuplelist(AllArcs)


def solveSSP(capexpvals, demvals):

    ssp = Model("secondstage")
    ssp.params.logtoconsole=0

    ### Second stage vars, flows on arcs 
    flow ={}
    for (i,j) in AllArcs:
        flow[i,j]=ssp.addVar(name='Flow%s_%s' % (i,j))

    ### Second stage vars, unmet demands
    unmet = {}
    for i in Cset:
        unmet[i] = ssp.addVar(obj=float(unmetCost[i]), name='Unmet%s' %(i))

    ssp.modelSense = GRB.MINIMIZE
    ssp.update()

    ### Arc Capacity constraints 
    flowbnds = {}
    for (i,j) in AllArcs:
        flowbnds[i,j] = ssp.addConstr(
                flow[i,j] <= curArcCap[i,j] + capexpvals[i,j], name='Capacity%s_%s' %(i,j))
    
    ### Facility capacity 
    faccap = {}
    for i in Fset:
        faccap[i] = ssp.addConstr(
                quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) <= facCap[i], name='FacCap%s' %i)
    
    ### Warehouse balance
    for i in Hset:
        ssp.addConstr(
                quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) - quicksum(flow[j,i] for (j,i) in
				    AllArcs.select('*',i)) == 0, name='WhBal%s' % i)

    ### Customer demands
    custdem = {}
    for i in Cset:
        custdem[i] = ssp.addConstr(
                quicksum(flow[j,i] for (j,i) in AllArcs.select('*',i)) + unmet[i] >= demvals[i], name='CustDem%s' % i)
    
    
    ssp.update()
    ssp.optimize()

    rhs = 0.0 
    custdemduals = ssp.getAttr('pi', custdem)
    faccapduals = ssp.getAttr('pi', faccap) 
    flowbndduals = ssp.getAttr('pi', flowbnds)
    rhs = sum(custdemduals[i]*demScen[i] for i in Cset) + sum(facCap[i]*faccapduals[i] for i in Fset) + sum(flowbndduals[i,j]*curArcCap[i,j] for (i,j) in AllArcs) 
 
    ## Solve
    return [ssp.objVal, rhs, flowbndduals]


#### Now let's try stochastic approximation
theta = 1.0
curx = {}
for (i,j) in AllArcs:
    curx[i,j] = 0.0

#### Use theta/sqrt(n) step size
for it in range(3000):
    demScen = {}
    for i in Cset:
        r = numpy.random.uniform(0.0,1.0)        
        if r < probnodem[i]:
            demScen[i] = 0.0
        else:
            demScen[i] = numpy.random.normal(demmean[i],demstdev[i])
    [objval, rhs, coefs] = solveSSP(curx,demScen)
    
    step = theta/numpy.power((float(it)+1.0),0.75)
    for (i,j) in AllArcs:
        curx[i,j] = min(max(0.0, curx[i,j] - step*(coefs[i,j] + arcExpCost[i,j])),200.0)

for (i,j) in AllArcs:
    print i, ',', j, ': ', curx[i,j] 

### estimate an upper bound, and build linear objective for estimating lower bound
numeval = 1000
evalvals = numpy.zeros(numeval)
firststagecost = sum(curx[i,j]*arcExpCost[i,j] for (i,j) in AllArcs)
const = 0.0
lincoefs = {}
for (i,j) in AllArcs:
    lincoefs[i,j] = 0.0
for k in range(numeval):
    demScen = {}
    for i in Cset:
        r = numpy.random.uniform(0.0,1.0)        
        if r < probnodem[i]:
            demScen[i] = 0.0
        else:
            demScen[i] = numpy.random.normal(demmean[i],demstdev[i])
    [evalvals[k],rhs,coefs] = solveSSP(curx,demScen)
    evalvals[k] = evalvals[k]+firststagecost
    for (i,j) in AllArcs:
        lincoefs[i,j] += (coefs[i,j] + arcExpCost[i,j])
    const += rhs
  
print 'mean solution objval = ', numpy.mean(evalvals)
print 'stdev of solution objval = ', numpy.std(evalvals)
ubmean = numpy.mean(evalvals)
ubwidth = numpy.std(evalvals)/math.sqrt(numeval)*2.145
print 'ci on upper bound = [', ubmean-ubwidth, ',', ubmean+ubwidth, ']'

### estimate a lower bound by optimizing linear function over the given bounds
### optimal solution is trivial: if coef is positive, set to lower bound
### if coef is negative, set to upper bound 
lbobj = 0.0
for (i,j) in AllArcs:
    print 'coef', i, ',', j, '= ', lincoefs[i,j]
    if lincoefs[i,j] < 0.0:
        lbobj += lincoefs[i,j]*200.0

print 'estimated lower bound = ', (lbobj + const)/1000.0






