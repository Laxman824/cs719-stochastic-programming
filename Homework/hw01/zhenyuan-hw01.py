'''
Homework 1 of cs719 Stochastic Programmin
Zhenyuan Shen (zhenyuan@cs.wisc.edu)
2016.09.24-26
Three models: (1) Stochastic Programming with Extensive Form
              (2) Stochastic Programming with Bender Decomposition
              (3) Mean Value model
'''

import cPickle
from gurobipy import *

### Read data from file you choose: commont/uncomment to choose the different files
### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

# import os
# os.chdir("/Users/shenzhenyuan/PycharmProjects/cs719/Homework/hw01")

dfile = open('nd1041015.pdat','r')
# dfile = open('nd1510203000.pdat','r')

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

''' ================================  Model 1: SP with Extensive Form ================================='''
m = Model("ArcExpansionExtensiveForm")

''' Variables '''
# First stage vars, arc expansion amounts
arcExpAmount = {}
for arc in AllArcs:
    arcExpAmount[arc]=m.addVar(obj=arcExpCost[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

# Second stage vars, unmet demands (z_jk)
nscen = len(Sset)
unmet = {}
for c in Cset:
    for k in Sset:
        unmet[c,k] = m.addVar(obj=float(unmetCost[c])/nscen, name='Unmet_%s_%s' % (c,k))

# Intermediate vars, transport
transport = {}
for arc in AllArcs:
    for k in Sset:
        transport[arc[0], arc[1], k] = m.addVar(obj=0, name="Transport_%s_%s_%s" % (arc[0], arc[1], k))

m.modelSense = GRB.MINIMIZE
m.update()


''' Constraints '''
# Production constraints (sum_h {x_fhk} - b_f <= 0)
for f in Fset:
    for k in Sset:
        m.addConstr(
            quicksum(transport[f,h,k] for h in Hset) <= facCap[f], name='Capacity_%s' % (f))

# Demand constraints (sum_i {y_ijk} + z_jk = d_jk)
for c in Cset:
    for k in Sset:
        m.addConstr(
            quicksum(transport[h,c,k] for h in Hset) + unmet[c,k] == demScens[(c,k)],name='Demand_%s_%s' % (c,k))


# Line capacity constraints
for arc in AllArcs:
    for k in Sset:
        m.addConstr(transport[arc[0], arc[1], k] <= curArcCap[arc] + arcExpAmount[arc],
                    name='LineCapacity_%s_%s_%s' % (arc[0], arc[1], k))

# Transportation contraints
for h in Hset:
    for k in Sset:
        m.addConstr(
            # quicksum(transport[f,h,k] for f in Fset) >= quicksum(transport[h,c,k] for c in Cset),
            # name="Transportation_%s_%s" % (h, k))
            quicksum(transport[f, h, k] for f, h in FHArcs.select('*', h)) >= quicksum(
                transport[h, c, k] for h, c in HCArcs.select(h, '*')),
            name="Transportation_%s_%s" % (h, k))

m.update()


print('\n-------------------------- Extensive Form  ------------------------------')
m.optimize()

if m.status == GRB.Status.OPTIMAL:

    print('EXPECTED COST : %g' % m.objVal)
    # print('SOLUTION:')
    # for k in Sset:
    #     print('  Senario: %s' % k)
    #     for arc in AllArcs:
    #         print('    [%s -> %s] (Transport, Current, Expanded) = (%g, %g, %g)'
    #               % (arc[0], arc[1],
    #                  transport[arc[0], arc[1], k].x,
    #                  arcExpAmount[arc[0], arc[1]].x,
    #                  curArcCap[arc]))
    #
    #
    # print('  AVERAGE UNMET DEMAND:')
    # for c in Cset:
    #     avgunmet = quicksum(unmet[c,s].x for s in Sset)/nscen
    #     print('   Customer %s: %g' % (c, avgunmet.getValue()))

spObjValue = m.objVal

''' ================================  Model 2: SP with Bender Decomposition ================================='''
# Master problem and Value function decision variables
master = Model("master")
master.params.logtoconsole=0

arcExpAmount_mp = {}
for arc in AllArcs:
    arcExpAmount_mp[arc] = master.addVar(obj=arcExpCost[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

theta = {}
for k in Sset:
    theta[k] = master.addVar(vtype=GRB.CONTINUOUS, obj=1.0/nscen, name="Theta_%s" % k)

master.modelSense = GRB.MINIMIZE
master.update()


# Subproblem and subproblem decision variables
sub = Model("subproblem")
sub.params.logtoconsole=0

unmet_sub = {}
for c in Cset:
    unmet_sub[c] = sub.addVar(obj=float(unmetCost[c]), name='Unmet_%s' % (c))

# Intermediate vars, transport
transport_sub = {}
for arc in AllArcs:
    transport_sub[arc] = sub.addVar(obj=0, name="Transport_%s_%s" % (arc[0], arc[1]))

sub.modelSense = GRB.MINIMIZE
sub.update()

''' Subproblem production constraints '''
# Production constraints (sum_h {x_fh} - b_f <= 0)
prodcon = {}
for f in Fset:
    prodcon[f] = sub.addConstr(
        quicksum(transport_sub[f,h] for h in Hset) <= facCap[f], name='Capacity_%s' % (f))

# Demand constraints (sum_h {y_hc} + z_c = d_c)
# For now just use scenario 0 data -- it will be reset during algorithm as it loops through scnearios
demcon = {}
for c in Cset:
    demcon[c] = sub.addConstr(
        quicksum(transport_sub[h,c] for h in Hset) + unmet_sub[c] == demScens[(c,Sset[0])],name='Demand_%s' % (c))


# Line capacity constraints
# For now just set right-hand side to capacity[p] -- it will be reset during algorithm
# based on master problem solution
capcon = {}
for arc in AllArcs:
    capcon[arc] = sub.addConstr(transport_sub[arc] <= curArcCap[arc],
                name='LineCapacity_%s_%s' % (arc[0], arc[1]))

# Transportation contraints
transcon = {}
for h in Hset:
    transcon[h] = sub.addConstr(
        quicksum(transport_sub[f,h] for f in Fset) >= quicksum(transport_sub[h,c] for c in Cset),
        name="Transportation_%s" % (h))

sub.update()


''' Begin the cutting plane loop '''
cutfound = 1  ## keep track if any violated cuts were found
iter = 1
lb = 0
ub = float("inf")
nCuts = 0
print('\n-------------------------- Bender Decomposition ------------------------------')
while cutfound:

    print '--- Iteration ', iter, ' ---'
    iter = iter+1
    # Solve current master problem
    cutfound = 0
    master.update()
    master.optimize()

    assert master.status == GRB.Status.OPTIMAL
    # print 'current optimal solution:'
    # for arc in AllArcs:
    #     print 'x[', arc[0], ' ', arc[1], ']=', arcExpAmount_mp[arc].x
    # for k in Sset:
    #     print 'theta[', k, ']=', theta[k].x
    # print 'objval = ', master.objVal
    lb = master.objVal

	 # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
    for arc in AllArcs:
        capcon[arc].RHS = curArcCap[arc] + arcExpAmount_mp[arc].x # Not arcExpAmount_mp[arc]

    Q = {}
    for k in Sset:
        for c in Cset:
            demcon[c].RHS = demScens[(c, k)]
        sub.update()
        sub.optimize()

        assert sub.status == GRB.Status.OPTIMAL

        # Display info, compute Benders cut, display, add to master
        # print 'sub[', k, '] objval = ', sub.objVal
        # for c in Cset:
        #     print 'unmet_sub[', c, ']=', unmet_sub[c].x

        if sub.objVal > theta[k].x + 0.000001:  ### violation tolerance
            rhs = 0.0
            for c in Cset:
                rhs += demScens[(c, k)] * demcon[c].Pi

            for f in Fset:
                rhs += facCap[f] * prodcon[f].Pi

            xcoef = {}
            for arc in AllArcs:
                xcoef[arc] = capcon[arc].Pi
                rhs += capcon[arc].Pi * curArcCap[arc]
                # print 'xcoef[', arc[0], ' ', arc[1], ']=', xcoef[arc]
                # print

            # print 'rhs = ', rhs
            master.addConstr(theta[k] - quicksum(xcoef[arc] * arcExpAmount_mp[arc] for arc in AllArcs) >= rhs)
            nCuts += 1
            cutfound = 1
        Q[k] = sub.objVal

    ub = quicksum(arcExpCost[arc] * arcExpAmount_mp[arc].x for arc in AllArcs).getValue() \
         + quicksum(1.0/nscen * Q[k] for k in Sset).getValue()
    print("    [lowerBound, upperBound] = [%f, %f],  nCut = %d" % (lb, ub, nCuts))

print('\nEXPECTED COST : %g' % master.objVal)

''' ================================  Model 3: Mean Value Model ================================='''

print('\n-------------------------- Mean Value Model ------------------------------')
mvm = Model("mvmodel")
mvm.params.logtoconsole=0
meanDemScens = {}
for c in Cset:
    meanDemScens[c] = quicksum(demScens[c,i] for i in Sset).getValue() / nscen

''' Variables '''
# First stage vars, arc expansion amounts
arcExpAmount_mv = {}
for arc in AllArcs:
    arcExpAmount_mv[arc]=mvm.addVar(obj=arcExpCost[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

# unmet demands (z_j)
unmet_mv = {}
for c in Cset:
    unmet_mv[c] = mvm.addVar(obj=unmetCost[c], name='Unmet_%s' % (c))

# Intermediate vars, transport
transport_mv = {}
for arc in AllArcs:
    transport_mv[arc] = mvm.addVar(obj=0, name="Transport_%s_%s" % (arc[0], arc[1]))

mvm.modelSense = GRB.MINIMIZE
mvm.update()


''' Constraints '''
# Production constraints (sum_h {x_fh} - b_f <= 0)
for f in Fset:
    mvm.addConstr(
        quicksum(transport_mv[f,h] for h in Hset) <= facCap[f], name='Capacity_%s' % (f))

# Demand constraints (sum_i {y_ij} + z_j = meand_j)
for c in Cset:
    mvm.addConstr(
        quicksum(transport_mv[h,c] for h in Hset) + unmet_mv[c] == meanDemScens[c],name='Demand_%s' % (c))


# Line capacity constraints
for arc in AllArcs:
    mvm.addConstr(transport_mv[arc] <= curArcCap[arc] + arcExpAmount_mv[arc],
                name='LineCapacity_%s_%s' % (arc[0], arc[1]))

# Transportation contraints
for h in Hset:
    mvm.addConstr(
        quicksum(transport_mv[f, h] for f, h in FHArcs.select('*', h)) >= quicksum(
            transport_mv[h, c] for h, c in HCArcs.select(h, '*')),
        name="Transportation_%s" % (h))

mvm.update()

''' Solve '''
mvm.optimize()


if mvm.status == GRB.Status.OPTIMAL:
    print('EXPECTED COST of deterministic model: %g' % mvm.objVal)
    # print('SOLUTION:')
    # for arc in AllArcs:
    #     print('    [%s -> %s] (Transport, Current, Expanded) = (%g, %g, %g)'
    #           % (arc[0], arc[1], transport_mv[arc].x, arcExpAmount_mv[arc].x, curArcCap[arc]))
    #
    #
    # print('  AVERAGE UNMET DEMAND:')
    # for c in Cset:
    #     print('   Customer %s: %g' % (c, unmet_mv[c].x))


############# Fix this as first-stage solution in stochastic model    ##########
for arc in AllArcs:
    arcExpAmount[arc].ub = arcExpAmount_mv[arc].x
    arcExpAmount[arc].lb = arcExpAmount_mv[arc].x

m.update()
m.optimize()
print('\nEXPECTED COST of MV SOLUTION: %g' % m.objVal)
print('\nThe value of stochastic solution: %g' % (m.objVal - spObjValue))
