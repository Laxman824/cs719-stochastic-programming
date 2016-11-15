import cPickle
from gurobipy import *

### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('forest4-10-31.pdat','r')

B = cPickle.load(dfile)  # set of base locations (list of strings) 
F = cPickle.load(dfile)  # set of fire districts (list of strings)
S = cPickle.load(dfile)  # set of scenarios (list of strings0
c = cPickle.load(dfile)  # dictionary of cost parameters c[i,j] is unit cost to relocate from i in B to j in B 
h = cPickle.load(dfile)  # dictionary or purchase params, h[i] is unit cost to add new resource at i in B
init = cPickle.load(dfile)  # dictionary of initial resource placements, init[i] is amt of resource at i in B
closesets = cPickle.load(dfile)  # dictionary of "close enough sets". closesets[f] is a list of base locations (subset of B)
demscens = cPickle.load(dfile)  # dictionary of demand scnearios. demscens[s,f] is demand of resource in scenario s at district f 
costscens = cPickle.load(dfile) # dictionary of cost scenarios. costscnes[s,f] is cost of resource shortage in scenario s at district f
budget = 500.0   ### hard coded here!!!
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

### This is just a check of the data. Probably you want to comment/delete these lines once you see the structure 

print B 
print F 
print S 
print c 
print h
print init
print closesets
print demscens
print costscens
print budget

##### Build the Master Model #####
m = Model("master")
m.params.logtoconsole=0

### repositioning variables 
repos = {} 
for i in B:
    for j in B:
        repos[i,j]=m.addVar(name='Repos%s_%s' %(i,j))

### capacity addition variables
add = {}
for i in B:
    add[i]=m.addVar(name='Add%s' %i)

### total avail capacity variables (not necessary, but simplifies cuts)
totavail = {}
for i in B:
    totavail[i]=m.addVar(name='TotAval%s' %i)

### Variable for representing second-stage cost
theta = m.addVar(name="theta", obj=1.0)

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
           quicksum(service[i,f] for f in closedists[i]) <= 0.0)

### Meet or record unmet demand in each scenario
demc = {}
for f in F:
    ## rhs needs to be updated
    demc[f] = subprob.addConstr(
            quicksum(service[i,f] for i in closesets[f])+unmet[f] == 0.0)

subprob.update()

iter = 1
lambdap = 0.2929

### Set as initial solution, totavail = what is currently there
curavail = {}
for i in B:
    curavail[i] = init[i]
currepos = {}
for i in B:
    for j in B:
        if i == j:
            currepos[i,j] = init[i]
        else:
            currepos[i,j] = 0.0
curadd = {}
for i in B:
    curadd[i] = 0.0

### solve subproblems and calculate objective value and cut for initial solution
sscost = 0.0
cutconst = 0.0
cutcoefs = {} ### only coefs on the totavail vars
for i in B:
    cutcoefs[i] = 0.0
for k in S:
    for f in F:
         demc[f].RHS = 0.0 
         unmet[f].obj = 0.0
    for (k,f) in SFkeys.select(k,'*'):
         demc[f].RHS = demscens[k,f]
         unmet[f].obj = costscens[k,f]

    subprob.update()
    subprob.optimize()

    sscost += subprob.objval
    for i in B:        
        cutcoefs[i] += servlim[i].pi 
    cutconst += sum(demc[f].pi*demscens[k,f] for (k,f) in SFkeys.select(k,'*'))
  
curub = sscost/float(len(S))
print 'upper bound = ', sscost/float(len(S)) 
m.addConstr(float(len(S))*theta >= cutconst + quicksum(totavail[i]*cutcoefs[i] for i in B))
curlb = curub - 1.0 # just to get loop started

while curub > curlb + 0.00001:

    print '================ Iteration ', iter, ' ==================='
    iter = iter+1
    ncuts = 0
    # Solve current master lower bounding problem
    theta.ub = curub ### (really there should be no upper bound on theta here, but this is a safe upper bound)
    m.setObjective(theta,GRB.MINIMIZE)
    m.update()
    m.optimize() 

    print 'lower bound = ', m.objval
    curlb = m.objval

    level = curlb + lambdap*(curub - curlb)

	 ### Now change master problem objective and set upper bound on second-stage cost
    theta.ub = level
    m.setObjective(quicksum((curavail[i] - totavail[i])*(curavail[i] - totavail[i]) for i in B)
	                +quicksum((curadd[i] - add[i])*(curadd[i]-add[i]) for i in B)
						 +quicksum((currepos[i,j]-repos[i,j])*(currepos[i,j]-repos[i,j]) for i in B for j in B),GRB.MINIMIZE)
    m.update()
    m.optimize()

	 ### now solve subproblems to eval at solution of modified master and find upper bound

    for i in B:
        servlim[i].RHS = totavail[i].x

    cutconst = 0.0
    sscost = 0.0
    cutcoefs = {} ### only coefs on the totavail vars
    for i in B:
        cutcoefs[i] = 0.0
    for k in S:
        for f in F:
             demc[f].RHS = 0.0 
             unmet[f].obj = 0.0
        for (k,f) in SFkeys.select(k,'*'):
             demc[f].RHS = demscens[k,f]
             unmet[f].obj = costscens[k,f]

        subprob.update()
        subprob.optimize()

        sscost += subprob.objval
        for i in B:        
            cutcoefs[i] += servlim[i].pi 
        cutconst += sum(demc[f].pi*demscens[k,f] for (k,f) in SFkeys.select(k,'*'))
 
    if sscost/float(len(S)) < curub:
        curub = sscost/float(len(S))
        for i in B:
            curavail[i] = totavail[i].x
            curadd[i] = add[i].x
            for j in B:
                currepos[i,j] = repos[i,j].x
    print 'upper bound = ', curub 
    if sscost > float(len(S))*theta.x + 0.000001:
         m.addConstr(float(len(S))*theta >= cutconst + quicksum(totavail[i]*cutcoefs[i] for i in B))

print('SOLUTION:')
for i in B:
    for j in B:
        if repos[i,j].x > 0.00001:
             print('Reposition from %s to %s: %g' % (i,j,repos[i,j].x))
for i in B:
    if add[i].x > 0.00001:
        print('Add at base %s: %g' % (i, add[i].x))

