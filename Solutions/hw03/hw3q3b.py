import cPickle
from gurobipy import *
import time

### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('forest-rev50-50-775.pdat','r')

B = cPickle.load(dfile)  # set of base locations (list of strings) 
F = cPickle.load(dfile)  # set of fire districts (list of strings)
S = cPickle.load(dfile)  # set of scenarios (list of strings0
c = cPickle.load(dfile)  # dictionary of cost parameters c[i,j] is unit cost to relocate from i in B to j in B 
h = cPickle.load(dfile)  # dictionary or purchase params, h[i] is unit cost to add new resource at i in B
init = cPickle.load(dfile)  # dictionary of initial resource placements, init[i] is amt of resource at i in B
closesets = cPickle.load(dfile)  # dictionary of "close enough sets". closesets[f] is a list of base locations (subset of B)
demscens = cPickle.load(dfile)  # dictionary of demand scnearios. demscens[s,f] is demand of resource in scenario s at district f 
costscens = cPickle.load(dfile) # dictionary of cost scenarios. costscnes[s,f] is cost of resource shortage in scenario s at district f
budget = 470.0   ### hard coded here!!!
dfile.close() 

### NoTE: demscens and costscens are very sparse. In particular, for a given s, there are only one or two districts f
### that have nonzero demscens vals (and only these have costscens vals, since without demand there is no need for a cost param)
### Here I define the set of keys which exist in these sets

SFkeys = demscens.keys()
SFkeys = tuplelist(SFkeys)  ## make tuplelist for easy selection

### it may also be useful to have a "reverse" version of the dictionary closesets, which provides for each facility i in
### B, the set of districts it can serve. This is constructed here
closedists = {}
for i in B:
    closedists[i] = []

for f in F:
    for i in closesets[f]:
        closedists[i].append(f)



##### Build the Master Model #####
m = Model("ben-master")
#m.params.logtoconsole=0
start = time.time()

### repositioning variables 
repos = {} 
for i in B:
    for j in B:
        repos[i,j]=m.addVar(vtype=GRB.CONTINUOUS, name='Repos%s_%s' %(i,j))

### capacity addition variables
add = {}
for i in B:
    add[i]=m.addVar(vtype=GRB.CONTINUOUS, name='Add%s' %i)

### total avail capacity variables (not necessary, but simplifies cuts)
totavail = {}
for i in B:
    totavail[i]=m.addVar(vtype=GRB.CONTINUOUS, name='TotAval%s' %i)

### Variable for representing second-stage cost
theta = {}
for s in S:
    theta[s] = m.addVar(name='theta%s' %s, obj=1.0/float(len(S)))

m.modelSense = GRB.MINIMIZE

m.update()

## only reposition what is available
for i in B:
    m.addConstr(quicksum(repos[i,j] for j in B) == init[i])

### budget on repoisitioning and purchasing
m.addConstr(quicksum(c[i,j]*repos[i,j] for i in B for j in B) +quicksum(h[i]*add[i] for i in B) <= budget)

### define the totavail vars
for i in B:
    m.addConstr(totavail[i] == add[i] + quicksum(repos[j,i] for j in B))


#### Now build the subproblem

subprob = Model("subprob")
subprob.params.logtoconsole=0

### Second stage vars, service  
service ={}
for f in F:
    for i in B:
        service[i,f]=subprob.addVar(name='Service%s_%s' % (i,f))

### Second stage vars, unmet demands
unmet = {}
for f in F:
    ### objective coefficient needs to be updated when re-solving problem
    unmet[f] = subprob.addVar(obj=1.0, name='Unmet%s' %f)

subprob.update()
subprob.modelSense = GRB.MINIMIZE

### service only what is available in each scenario
servlim = {}
for i in B:
    ## rhs needs to be updated
    servlim[i] = subprob.addConstr(
           quicksum(service[i,f] for f in closedists[i]) <= 1.0)

### Meet or record unmet demand in each scenario
demc = {}
for f in F:
    ## rhs needs to be updated
    demc[f] = subprob.addConstr(
            quicksum(service[i,f] for i in closesets[f])+unmet[f] == 0.0)

subprob.update()

cutfound = 0   ### when set = 1, does the phase 0, otherwise skips it
iter = 1

while cutfound:
    print '================ Iteration ', iter, ' ==================='
    iter = iter+1
    ncuts = 0
    # Solve current master problem
    cutfound = 0 
    m.update()
    m.optimize() 

    print 'lower bound = ', m.objval

    for i in B:
        servlim[i].RHS = totavail[i].x

    cutconst = 0.0
    sscost = 0.0
    ub = 0.0
    for k in S:
        for f in F:
             demc[f].RHS = 0.0 
             unmet[f].obj = 0.0
        for (k,f) in SFkeys.select(k,'*'):
             demc[f].RHS = demscens[k,f]
             unmet[f].obj = costscens[k,f]

        subprob.update()
        subprob.optimize()
        if subprob.objval > theta[k].x + 0.000001:
            cutfound = 1
            cutconst = sum(demc[f].pi*demscens[k,f] for (k,f) in SFkeys.select(k,'*'))
            m.addConstr(theta[k] >= cutconst + quicksum(totavail[i]*servlim[i].pi for i in B)) 
        
        ub += subprob.objval
  
    print 'upper bound = ', ub/float(len(S)) 


### Now change variable types to integer
for i in B:
    for j in B:
        repos[i,j].vtype = GRB.INTEGER

### capacity addition variables
for i in B:
    add[i].vtype = GRB.INTEGER

### total avail capacity variables (not necessary, but simplifies cuts)
for i in B:
    totavail[i].vtype = GRB.INTEGER


### Define the callback function
def BendersCallback(model, where):
    if where == GRB.Callback.MIPSOL:
        for i in B:
            servlim[i].RHS = model.cbGetSolution(totavail[i])

        cutconst = 0.0
        sscost = 0.0
        ub = 0.0
        for k in S:
            for f in F:
                demc[f].RHS = 0.0 
                unmet[f].obj = 0.0
            for (k,f) in SFkeys.select(k,'*'):
                demc[f].RHS = demscens[k,f]
                unmet[f].obj = costscens[k,f]

            subprob.update()
            subprob.optimize()
            if subprob.objval > model.cbGetSolution(theta[k]) + 0.000001:
                cutfound = 1
                cutconst = sum(demc[f].pi*demscens[k,f] for (k,f) in SFkeys.select(k,'*'))
                model.cbLazy(theta[k] >= cutconst + quicksum(totavail[i]*servlim[i].pi for i in B)) 

m.Params.lazyConstraints = 1

m.optimize(BendersCallback)


print('\nEXPECTED COST : %g' % m.objVal)
print('SOLUTION:')
for i in B:
    for j in B:
        if repos[i,j].x > 0.00001:
             print('Reposition from %s to %s: %g' % (i,j,repos[i,j].x))
for i in B:
    if add[i].x > 0.00001:
        print('Add at base %s: %g' % (i, add[i].x))

print 'total time = ', time.time() - start

