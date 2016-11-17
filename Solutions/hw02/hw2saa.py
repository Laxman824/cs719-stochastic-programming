import cPickle
from gurobipy import *
import time
import numpy
import math

numpy.random.seed(423)

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

print demmean
print demstdev

### Define sets of arcs (used as keys to dictionaries)
FHArcs = [(i,j) for i in Fset for j in Hset]  ## arcs from facilities to warehouses
HCArcs = [(i,j) for i in Hset for j in Cset]   ## arcs from warehouses to customers
AllArcs = FHArcs + HCArcs
print AllArcs

### Make them Gurobi tuplelists
FHArcs = tuplelist(FHArcs)
HCArcs = tuplelist(HCArcs)
AllArcs = tuplelist(AllArcs)


## implement SAA below here
numscen = 50
numbatch = 15
numeval = 1000

## this function builds and solves the extensive form for a given set of scenarios, returns objval and solution
def solveScenModel(demScens,Sset):

    m = Model("netdesign")
    m.params.logtoconsole=0
    
    ### Arc expansion variables 
    capAdd = {} 
    for (i,j) in AllArcs:
        capAdd[i,j]=m.addVar(obj=arcExpCost[i,j], name='CapExp%s_%s' % (i,j))
    
    ### Second stage vars, flows on arcs 
    flow ={}
    for (i,j) in AllArcs:
        for k in Sset:
            flow[i,j,k]=m.addVar(name='Flow%s_%s_%s' % (i,j,k))
    
    ### Second stage vars, unmet demands
    unmet = {}
    for i in Cset:
        for k in Sset:
            unmet[i,k] = m.addVar(obj=float(unmetCost[i])/float(len(Sset)), name='Unmet%s_%s' %(i,k))
    
    m.modelSense = GRB.MINIMIZE
    
    m.update()
    
    
    ### Arc Capacity constraints 
    for (i,j) in AllArcs:
        for k in Sset:
            m.addConstr(
                flow[i,j,k] <= curArcCap[i,j] + capAdd[i,j], name='Capacity%s_%s_%s' %(i,j,k))
    
    ### Facility capacity 
    for i in Fset:
        for k in Sset:
            m.addConstr(
                quicksum(flow[i,j,k] for (i,j) in AllArcs.select(i,'*')) <= facCap[i], name='FacCap%s_%s' %(i,k))
    
    ### Warehouse balance
    for i in Hset:
        for k in Sset:
            m.addConstr(
                quicksum(flow[i,j,k] for (i,j) in AllArcs.select(i,'*')) - quicksum(flow[j,i,k] for (j,i) in
				    AllArcs.select('*',i)) == 0, name='WhBal%s_%s' % (i,k))
    
    ### Customer demands
    for i in Cset:
        for k in Sset:
            m.addConstr(
                quicksum(flow[j,i,k] for (j,i) in AllArcs.select('*',i)) + unmet[i,k] >= demScens[i,k],
				    name='CustDem%s_%s' % (i,k))
    
    
    m.update()
    
    ## Solve
    m.optimize()
    capexpvals = {}
    for (i,j) in AllArcs:
        capexpvals[i,j] = capAdd[i,j].x
    return [m.objVal, capexpvals]

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
    for (i,j) in AllArcs:
        ssp.addConstr(
                flow[i,j] <= curArcCap[i,j] + capexpvals[i,j], name='Capacity%s_%s' %(i,j))
    
    ### Facility capacity 
    for i in Fset:
        ssp.addConstr(
                quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) <= facCap[i], name='FacCap%s' %i)
    
    ### Warehouse balance
    for i in Hset:
        ssp.addConstr(
                quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) - quicksum(flow[j,i] for (j,i) in
				    AllArcs.select('*',i)) == 0, name='WhBal%s' % i)

    ### Customer demands
    for i in Cset:
        ssp.addConstr(
                quicksum(flow[j,i] for (j,i) in AllArcs.select('*',i)) + unmet[i] >= demvals[i], name='CustDem%s' % i)
    
    
    ssp.update()
    
    ## Solve
    ssp.optimize()
    return ssp.objVal

### Here is where the SAA, independent lower and upper bound occurs

objvals = numpy.zeros(numbatch)
bestval = 1e10 
bestx = {}

for k in range(numbatch):

    ### set up and solve extensive form for a sample average approximation problem
    demScens = {}
    Sset = range(numscen)
    for s in range(numscen):
        for i in Cset:
            r = numpy.random.uniform(0.0,1.0)        
            if r < probnodem[i]:
				    demScens[(i,s)] = 0.0
            else:
                demScens[(i,s)] = numpy.random.normal(demmean[i],demstdev[i])
    [objvals[k],candx] = solveScenModel(demScens,Sset) 
    if objvals[k] < bestval:
        bestval = objvals[k]
        bestx = candx

### Note, critical value is for 95% 2-sided confidence interval (tau, 14 degrees of freedom, p=0.975)
print 'mean SAA objval = ', numpy.mean(objvals)
print 'stdev of SAA objval = ', numpy.std(objvals)
lbmean = numpy.mean(objvals)
lbwidth = numpy.std(objvals)/math.sqrt(numbatch)*2.145

### now independently estimate upper bound for the "bestx" candidate solution

evalvals = numpy.zeros(numeval)
firststagecost = sum(candx[i,j]*arcExpCost[i,j] for (i,j) in AllArcs)
for k in range(numeval):
    demScen = {}
    for i in Cset:
        r = numpy.random.uniform(0.0,1.0)        
        if r < probnodem[i]:
            demScen[i] = 0.0
        else:
            demScen[i] = numpy.random.normal(demmean[i],demstdev[i])
    evalvals[k] = solveSSP(candx,demScen)+firststagecost
    
  
print 'mean solution objval = ', numpy.mean(evalvals)
print 'stdev of solution objval = ', numpy.std(evalvals)
ubmean = numpy.mean(evalvals)
ubwidth = numpy.std(evalvals)/math.sqrt(numeval)*2.145

print 'ci on lower bound = [', lbmean-lbwidth, ',', lbmean+lbwidth, ']'
print 'ci on upper bound = [', ubmean-ubwidth, ',', ubmean+ubwidth, ']'


### Finally, let's do the common random numbers approach to estimate gap

### first generate a candidate x
demScens = {}
Sset = range(numscen)
for s in range(numscen):
    for i in Cset:
        r = numpy.random.uniform(0.0,1.0)        
        if r < probnodem[i]:
            demScens[(i,s)] = 0.0
        else:
            demScens[(i,s)] = numpy.random.normal(demmean[i],demstdev[i])
[objval,candx] = solveScenModel(demScens,Sset) 
firststagecost = sum(candx[i,j]*arcExpCost[i,j] for (i,j) in AllArcs)

### now estimate the gap
gapvals = numpy.zeros(numbatch)
objvals = numpy.zeros(numbatch*numscen)
for k in range(numbatch):

    ### set up and solve extensive form for a sample average approximation problem
    demScens = {}
    Sset = range(numscen)
    for s in range(numscen):
        for i in Cset:
            r = numpy.random.uniform(0.0,1.0)        
            if r < probnodem[i]:
				    demScens[(i,s)] = 0.0
            else:
                demScens[(i,s)] = numpy.random.normal(demmean[i],demstdev[i])
    [objval,curcand] = solveScenModel(demScens,Sset) 

    curvals = numpy.zeros(numscen)
    for s in range(numscen):
        demScen = {}
        for i in Cset:
            demScen[i] = demScens[(i,s)]
        curvals[s] = solveSSP(candx,demScen)+firststagecost
        objvals[k*numscen+s]=curvals[s]

    gapvals[k] = numpy.mean(curvals) - objval 

### Note, critical value is for 95% 1-sided confidence interval (tau, 14 degrees of freedom, p=0.95)
print 'estimated objval of solution', numpy.mean(objvals)
print '95% c.i. on solution obj val: [', numpy.mean(objvals)-numpy.std(evalvals)/math.sqrt(numeval)*2.145, ',', numpy.mean(objvals)+numpy.std(evalvals)/math.sqrt(numeval)*2.145, ']'
print 'mean gap estimate = ', numpy.mean(gapvals)
print '95% c.i. on gap = [0, ', numpy.mean(gapvals) + numpy.std(gapvals)/math.sqrt(numbatch)*1.761, ']'

