import cPickle
from gurobipy import *
import time

### Read data from file you choose: commont/uncomment to choose the different files
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('nd10-4-10-15.pdat','r')
#dfile = open('nd15-10-20-3000.pdat','r')

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

print Fset
print Hset
print Cset
print Sset 
#print arcExpCost
print facCap
#print curArcCap
print unmetCost
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
stochobjval = m.objVal

print('\nEXPECTED COST : %g' % m.objVal)
print('SOLUTION:')
for (i,j) in AllArcs:
    if capAdd[i,j].x > 0.00001:
        print('Arc %s,%s expanded by %g' % (i,j,capAdd[i,j].x))
print('AVERAGE UNMET DEMAND:') 
for i in Cset:
    avgunmet = quicksum(unmet[i,k].x for k in Sset)/(len(Sset))
    print('   Customer %s: %g' % (i, avgunmet))


print 'total time for building and solving extensive form = ', time.time() - start



### Now build mean value problem

mvm = Model("mvnetdesign")

### Arc expansion variables 
mvcapAdd = {} 
for (i,j) in AllArcs:
    mvcapAdd[i,j]=mvm.addVar(obj=arcExpCost[i,j], name='CapExp%s_%s' % (i,j))

### Second stage vars, flows on arcs 
mvflow ={}
for (i,j) in AllArcs:
    mvflow[i,j]=mvm.addVar(name='Flow%s_%s' % (i,j))

### Second stage vars, unmet demands
mvunmet = {}
for i in Cset:
    mvunmet[i] = mvm.addVar(obj=float(unmetCost[i]), name='Unmet%s' %(i))

mvm.modelSense = GRB.MINIMIZE

mvm.update()


### Arc Capacity constraints 
for (i,j) in AllArcs:
    mvm.addConstr(
            mvflow[i,j] <= curArcCap[i,j] + mvcapAdd[i,j], name='Capacity%s_%s' %(i,j))

### Facility capacity 
for i in Fset:
    mvm.addConstr(
            quicksum(mvflow[i,j] for (i,j) in AllArcs.select(i,'*')) <= facCap[i], name='FacCap%s' %i)

### Warehouse balance
for i in Hset:
    mvm.addConstr(
            quicksum(mvflow[i,j] for (i,j) in AllArcs.select(i,'*')) - quicksum(mvflow[j,i] for (j,i) in
				AllArcs.select('*',i)) == 0, name='WhBal%s' % i)

### Customer demands
for i in Cset:
    mvm.addConstr(
            quicksum(mvflow[j,i] for (j,i) in AllArcs.select('*',i)) + mvunmet[i] >= quicksum(demScens[i,k] for k in
				Sset)/float(len(Sset)), name='CustDem%s' % i)


mvm.update()

## Solve
mvm.optimize()

print('\nMEAN VALUE PROBLEM COST : %g' % mvm.objVal)
print('SOLUTION:')
for (i,j) in AllArcs:
    if mvcapAdd[i,j].x > 0.00001:
        print('Arc %s,%s expanded by %g' % (i,j,mvcapAdd[i,j].x))


##### Finally, fix mena value solution in the stochastic model and re-solve stochastic model
for (i,j) in AllArcs:
    capAdd[i,j].lb = capAdd[i,j].ub = mvcapAdd[i,j].x

m.update()

m.optimize()

print('\nEXPECTED COST OF MEAN VALUE SOLUTION: %g' % m.objVal)
print('\nVALUE OF STOCHASTIC SOLUTION: %g' % (m.objVal - stochobjval))


