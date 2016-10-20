import cPickle
from gurobipy import *

### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('forest41031.pdat','r')

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


demscensDict = {}
for (i,j) in demscens:
    try:
        demscensDict[i].append(j)
    except KeyError:
        demscensDict[i] = [j]
for i in demscensDict:
    print i, demscensDict[i]

AllArcs = [(i,j) for i in B for j in B if j != i]
AllArcs = tuplelist(AllArcs)

# record all close relations
BF = [(b,f) for b in B for f in closedists[b]]
BF = tuplelist(BF)

nscen = len(S)



''' ================================  single-cut Benders decomposition ================================='''
''' Master problem's variable '''
master = Model("master")
master.params.logtoconsole=0

moveAmount_mp = {}
for arc in AllArcs:
    moveAmount_mp[arc] = master.addVar(obj=0.0, name='move_%s_%s' %(arc[0], arc[1]))

buyAmount_mp = {}
nunit_mp = {}
for b in B:
    buyAmount_mp[b] = master.addVar(obj=0.0, name='purchase_%s' %(b))
    nunit_mp[b] = master.addVar(obj=0.0, name='nunit_%s' % (b))

theta = master.addVar(vtype=GRB.CONTINUOUS, obj=1.0, name="Theta")


master.modelSense = GRB.MINIMIZE
master.update()

''' Master problem's constraints '''
# Budget constraints (sum_b {h_b x_b} + sum_arc c_ij z_ij = C)
budgetcon = master.addConstr(
    quicksum(buyAmount_mp[b]*h[b] for b in B) + quicksum(moveAmount_mp[arc]*c[arc] for arc in AllArcs) <= budget, name='Budget')

capacitycon = {}
for b in B:
    capacitycon[b] = master.addConstr(init[b] + buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in AllArcs.select('*', b))
                       - quicksum(moveAmount_mp[arc] for arc in AllArcs.select(b, '*')) == nunit_mp[b], name='Capacity_%s' %(b))

master.update()

''' Subproblem's variables '''
sub = Model("subproblem")
sub.params.logtoconsole=0

unmet_sub = {}
for f in F:
    unmet_sub[f] = sub.addVar(obj=0, name='Unmet_%s' % (f))

respond_sub = {}
for b, f in BF:
    respond_sub[b, f] = sub.addVar(obj=0, name="Respond_%s_%s" % (b, f))

sub.modelSense = GRB.MINIMIZE
sub.update()

''' Subproblem's constraints '''
# Capacity constraints (sum_fk {y_bfk} <= I_b + sum_j c_jb z_jb)
capcon = {}
for b in B:
    capcon[b] = sub.addConstr(
        quicksum(respond_sub[arc] for arc in BF.select(b, '*')) <= init[b], name='Capacity_%s' % (b))


# Demand constraints (sum_b {y_bfk} + u_fk = d_fk)
demcon = {}
for f in F:
    demcon[f] = sub.addConstr(
        quicksum(respond_sub[arc] for arc in BF.select('*', f)) + unmet_sub[f] == 0, name='Demand_%s' % (f))

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
    print 'current optimal solution:'
    for b in B:
        print 'X[' ,b ,'] = ', buyAmount_mp[b].x
    for arc in AllArcs:
        print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
    print 'objval = ', master.objVal
    lb = master.objVal

	 # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
    for b in B:
        capcon[b].RHS = nunit_mp[b].x

    Q = 0
    rhs = 0.0
    xzcoef = {}
    for k in S:
        for f in F:
            if (k, f) in costscens:
                unmet_sub[f].obj = costscens[k, f]
                demcon[f].RHS = demscens[k, f]
            else:
                unmet_sub[f].obj = 0
                demcon[f].RHS = 0
        sub.update()
        sub.optimize()

        assert sub.status == GRB.Status.OPTIMAL
        Q += sub.objVal

        # print " ----------------------- "
        # print k, sub.objVal
        #
        # # Display info, compute Benders cut, display, add to master
        # print 'sub[', k, '] objval = ', sub.objVal
        # for f in F:
        #     print 'unmet_sub[', f, '] = ', unmet_sub[f].x
        # for b, f in BF:
        #     print 'respond_sub[', b, ', ', f, '] = ', respond_sub[b,f].x
        # print len(respond_sub)

        for b in B:
            xzcoef[k, b] = capcon[b].Pi
            rhs += xzcoef[k,b] * init[b]

        for f in F:
            if (k, f) in demscens:
                rhs += demscens[k,f] * demcon[f].Pi

    assert lb == theta.x

    ub = Q
    if ub > lb + 0.000001:  ### violation tolerance
        master.addConstr(
            theta - 1.0 * quicksum(
                xzcoef[k, b] * (buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in AllArcs.select('*', b))
                                - quicksum(moveAmount_mp[arc] for arc in AllArcs.select(b, '*')))
                             for (k,b) in xzcoef) >= rhs)
        nCuts += 1
        cutfound = 1

    print("    [lowerBound, upperBound] = [%f, %f],  nCut = %d" % (lb, ub, nCuts))

print('\nEXPECTED COST : %g' % master.objVal)