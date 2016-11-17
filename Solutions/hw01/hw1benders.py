import cPickle
from gurobipy import *
import time

### Read data from file you choose: commont/uncomment to choose the different files
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

#dfile = open('nd10-4-10-15.pdat','r')
dfile = open('nd15-10-20-3000.pdat','r')

Fset = cPickle.load(dfile)  # set of facilities (list of strings)
Hset = cPickle.load(dfile)  # set of warehouses (list of strings)
Cset = cPickle.load(dfile)  # set of customers (list of strings)
Sset = cPickle.load(dfile)  # set of scenarios (list of strings)
arcExpCost = cPickle.load(dfile)  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
facCap = cPickle.load(dfile)   # facility capacities (dictionary mapping F to floats) 
curArcCap = cPickle.load(dfile)  # current arc capacities (dictionary mapping (i,j) to floats, where either 
                                 # i is facility, j is warehouse, or i is warehouse and j is customer
unmetCost = cPickle.load(dfile)  # penalty for unment customer demands (dicationary mapping C to floats)
demScens = cPickle.load(dfile)  # demand scenarios (dictionary mapping (i,k) tuples to floats, where i is customer, k is
                                #scenario
dfile.close() 

### This is just a check of the data. Probably you want to comment/delete these lines once you see the structure 

#print Fset
#print Hset
#print Cset
#print Sset 
#print arcExpCost
#print facCap
#print curArcCap
#print unmetCost
#print demScens


### Define sets of arcs (used as keys to dictionaries)
FHArcs = [(i,j) for i in Fset for j in Hset]  ## arcs from facilities to warehouses
HCArcs = [(i,j) for i in Hset for j in Cset]   ## arcs from warehouses to customers
AllArcs = FHArcs + HCArcs
print AllArcs

### Make them Gurobi tuplelists
FHArcs = tuplelist(FHArcs)
HCArcs = tuplelist(HCArcs)
AllArcs = tuplelist(AllArcs)


##### Build the Master Model #####
start = time.time()
master = Model("ben-master")

### Arc expansion variables 
capAdd = {} 
for (i,j) in AllArcs:
    capAdd[i,j]=master.addVar(obj=arcExpCost[i,j], name='CapExp%s_%s' % (i,j))

### Scenario cost variables
theta = {}
for k in Sset:
    theta[k]=master.addVar(obj=1.0/float(len(Sset)),name='Theta%s' % k)



### Build the Subproblem Model ####
subprob = Model("scen-prob")
subprob.params.logtoconsole=0  ## turns off display of Gurobi output when solving subproblems

flow ={}
for (i,j) in AllArcs:
    flow[i,j]=subprob.addVar(name='Flow%s_%s' % (i,j))

###  unmet demands
unmet = {}
for i in Cset:
    unmet[i] = subprob.addVar(obj=unmetCost[i], name='Unmet%s' %i)

subprob.modelSense = GRB.MINIMIZE

subprob.update()


### Arc Capacity constraints. For now just use current capacity (will be updated according to master problem solution)
flowbndsc = {}
for (i,j) in AllArcs:
    flowbndsc[i,j] = subprob.addConstr(flow[i,j] <= curArcCap[i,j], name='Capacity%s_%s' %(i,j))

### Facility capacity 
faccapc = {}
for i in Fset:
    faccapc[i] = subprob.addConstr(
            quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) <= facCap[i], name='FacCap%s' %i)

### Warehouse balance (no need to save these, as they do not contribute to Benders cuts b/c RHS = 0)
for i in Hset:
    subprob.addConstr(
            quicksum(flow[i,j] for (i,j) in AllArcs.select(i,'*')) - quicksum(flow[j,i] for (j,i) in
				AllArcs.select('*',i)) == 0, name='WhBal%s' % i)

### Customer demands. For now use avg demand. Will be updated in algorithm
custdemc = {}
for i in Cset:
    custdemc[i] = subprob.addConstr(
            quicksum(flow[j,i] for (j,i) in AllArcs.select('*',i)) + unmet[i] >= quicksum(demScens[i,k] for k in
				Sset)/float(len(Sset)), name='CustDem%s' % i)

subprob.update()

cutfound = 1
iter = 1

while cutfound:
    print '================ Iteration ', iter, ' ==================='
    iter = iter+1
    ncuts = 0
    # Solve current master problem
    cutfound = 0 
    master.update()
    master.optimize() 
 
    print 'current optimal solution:'
    for (i,j) in AllArcs:
        if capAdd[i,j].x > 0.00001:
            print('Arc %s,%s expanded by %g' % (i,j,capAdd[i,j].x))
    print 'objval = ', master.objVal
	 # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
    for (i,j) in AllArcs:
        flowbndsc[i,j].RHS = curArcCap[i,j] + capAdd[i,j].x

    for k in Sset:
        for i in Cset:
		      custdemc[i].RHS = demScens[i,k]
        subprob.update()	
        subprob.optimize()

        if subprob.objVal > theta[k].x + 0.00001:  ### violation tolerance
            rhs = 0.0 
            custdemduals = subprob.getAttr('pi', custdemc)
            faccapduals = subprob.getAttr('pi', faccapc) 
            flowbndduals = subprob.getAttr('pi', flowbndsc)
            rhs = quicksum(custdemduals[i]*demScens[i,k] for i in Cset) - quicksum(facCap[i]*faccapduals[i] for i in Fset) + quicksum(flowbndduals[i,j]*curArcCap[i,j] for (i,j) in AllArcs) 
            master.addConstr(theta[k] - quicksum(flowbndduals[i,j]*capAdd[i,j] for (i,j) in AllArcs) >= rhs)
            cutfound = 1 
            ncuts = ncuts + 1
    print 'num cuts found in this iteration = ', ncuts

print 'total time = ', time.time() - start
