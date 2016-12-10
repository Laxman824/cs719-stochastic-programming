import cPickle
from gurobipy import *
import time

### Read data from file you choose: commont/uncomment to choose the different files
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('nd15-10-20-500.pdat','r')

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

avarweight = 0.5
alpha = 0.95
dfile.close() 


### Define sets of arcs (used as keys to dictionaries)
FHArcs = [(i,j) for i in Fset for j in Hset]  ## arcs from facilities to warehouses
HCArcs = [(i,j) for i in Hset for j in Cset]   ## arcs from warehouses to customers
AllArcs = FHArcs + HCArcs

### Make them Gurobi tuplelists
FHArcs = tuplelist(FHArcs)
HCArcs = tuplelist(HCArcs)
AllArcs = tuplelist(AllArcs)


##### Start building the Model #####

start = time.time()
m = Model("netdesign")

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


### gamma var, for defining Average value at risk
gamma = m.addVar(obj=avarweight, lb=-GRB.INFINITY, name='Gamma')

### total cost vars (one per scenario, but includes first-stage costs)
costs = {}
for k in Sset:
    costs[k] = m.addVar(obj=0.0, name='cost%s' %k)

### shortfall vars for defining average value at risk
shortfall = {}
for k in Sset:
    shortfall[k] = m.addVar(obj=avarweight/((1.0-alpha)*float(len(Sset))), name='Short%s' %k)

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

### Total cost definitions
for k in Sset:
    m.addConstr(costs[k] == quicksum(unmet[i,k]*unmetCost[i] for i in Cset) + quicksum(capAdd[i,j]*arcExpCost[i,j] for
	 (i,j) in AllArcs))

### Shortfall definitions
for k in Sset:
    m.addConstr(shortfall[k] >= costs[k] - gamma)

m.update()

## Solve
m.optimize()
stochobjval = m.objVal

print('\nOBJECTIVE: %g' % m.objVal)
expandcost = sum(capAdd[i,j].x*arcExpCost[i,j] for (i,j) in AllArcs) 
print('\nEXPANSION COST: %g' % expandcost)
expectsscost = sum(unmetCost[i]*unmet[i,k].x for k in Sset for i in Cset)/float(len(Sset))
print('\nEXPECTED UNMET COST: %g' % expectsscost)
print('\nVaR (Gamma): %g' % gamma.x)
avar = gamma.x + 1.0/(float(len(Sset))*(1.0-alpha))*sum(shortfall[k].x for k in Sset)
print('\nAVaR: %g' % avar)


print('SOLUTION:')
for (i,j) in AllArcs:
    if capAdd[i,j].x > 0.00001:
        print('Arc %s,%s expanded by %g' % (i,j,capAdd[i,j].x))
print('AVERAGE UNMET DEMAND:') 
for i in Cset:
    avgunmet = sum(unmet[i,k].x for k in Sset)/(len(Sset))
    print('   Customer %s: %g' % (i, avgunmet))


print 'total time for building and solving extensive form = ', time.time() - start



