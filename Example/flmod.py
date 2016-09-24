import cPickle
from gurobipy import *

### Read data from file `fl.dat'
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

# import os
# os.chdir("/Users/shenzhenyuan/PycharmProjects/cs719/Example")

f = open('fl.dat','r')
cap = cPickle.load(f)
fixedc = cPickle.load(f)
costs = cPickle.load(f)
demscens = cPickle.load(f)
pen = cPickle.load(f)
f.close()

## get sizes of data (8 facilities to 10 customers in 100 scenarios)
fac = range(len(cap))
cust = range(len(demscens[0]))
scen = range(len(demscens))

print cap
print fixedc
print costs
print demscens
print pen
nscen = len(demscens)

##### Start building the Model #####

m = Model("facility")

### First stage vars, facility open decisions (x_i)
fopen = {} 
for i in fac:
    fopen[i]=m.addVar(vtype=GRB.BINARY, obj=fixedc[i], name='Open%d' % i)

### Second stage vars, ship from/to decisions (y_ijk)
### Assume all scenarios have probability 1/nscen, so scale objective coefficients by this
transport={}
for i in fac:
    for c in cust:
        for s in scen:
            transport[i,c,s] = m.addVar(obj=costs[i][c]/nscen, name='Transport%d_%d_%d' % (i, c, s))

### Second stage vars, unmet demands (z_jk)
unmet = {}
for c in cust:
    for s in scen:
        unmet[c,s] = m.addVar(obj=float(pen[c])/nscen, name='Unmet%d_%d' %(c,s))

m.modelSense = GRB.MINIMIZE

m.update()


### Production constraints (sum_j {y_ijk} - u_i * x_i <= 0)
for i in fac:
    for s in scen:
        m.addConstr(
            quicksum(transport[i,c,s] for c in cust) <= cap[i] * fopen[i], name='Capacity%d_%d' %(i,s))

### Demand constraints (sum_i {y_ijk} + z_jk = d_jk)
for c in cust:
    for s in scen:
        m.addConstr(
            quicksum(transport[i,c,s] for i in fac) + unmet[c,s] == demscens[s][c], name='Demand%d_%d' %(c,s))

m.update()

## Solve
m.optimize()

print('\nEXPECTED COST : %g' % m.objVal)
print('SOLUTION:')
for i in fac:
    if fopen[i].x > 0.99:
        print('Plant %s open' % i)
        for c in cust:
            if transport[i,c,0].x > 0.00001:
                print('  Transport %g units to customer %s in scenario 0' % \
                      (transport[i,c,0].x, c))
print('AVERAGE UNMET DEMAND:') 
for c in cust:
    avgunmet = quicksum(unmet[c,s].x for s in scen)/nscen

    print('   Customer %s: %g' % (c, avgunmet.getValue()))






############ Build and solve the mean value problem  ############

mvm = Model("avgfacility")

### First stage vars, facility open decisions
mvfopen = {} 
for i in fac:
    mvfopen[i]=mvm.addVar(vtype=GRB.BINARY, obj=fixedc[i], name='Open%d' % i)

### Ship from/to decisions 
### Assume all scenarios have probability 1/nscen, so scale objective coefficients by this
mvtransport={}
for i in fac:
    for c in cust:
        mvtransport[i,c] = mvm.addVar(obj=costs[i][c], name='Transport%d_%d' % (i, c))

### Unmet demands
mvunmet = {}
for c in cust:
    mvunmet[c] = mvm.addVar(obj=float(pen[c]), name='Unmet%d' %c)

mvm.modelSense = GRB.MINIMIZE

mvm.update()


### Production constraints
for i in fac:
    mvm.addConstr(
        quicksum(mvtransport[i,c] for c in cust) <= cap[i] * mvfopen[i], name='Capacity%d' % i )

### Demand constraints
for c in cust:
    mvm.addConstr(
        quicksum(mvtransport[i,c] for i in fac) + mvunmet[c] == quicksum(demscens[s][c] for s in scen)/nscen, name='Demand%d' % c)

mvm.update()

## Solve
mvm.optimize()
print('\nTOTAL COST OF MEAN VALUE MODEL: %g' % mvm.objVal)
print('SOLUTION:')
for i in fac:
    if mvfopen[i].x > 0.99:
        print('Plant %s open' % i)
        for c in cust:
            if mvtransport[i,c].x > 0.00001:
                print('  Transport %g units to customer %s 0' % \
                      (mvtransport[i,c].x, c))

print('UNMET DEMAND:') 
for c in cust:
    print('   Customer %s: %g' % (c, mvunmet[c].x))




############# Fix this as first-stage solution in stochastic model    ##########
for i in fac:
    if mvfopen[i].x > 0.99:
        fopen[i].lb = 1.0
    else:
        fopen[i].ub = 0.0

m.update()
m.optimize()
print('\nEXPECTED COST of MV SOLUTION: %g' % m.objVal)

