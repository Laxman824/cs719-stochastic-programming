from gurobipy import *

Cset = ['C1','C2','C3','C4','C5']
Sset = ['S1','S2']
Fset = ['F1','F2','F3','F4','F5']

#demand values in each scenario
demand = { ('C1','S1'): 13,    ### 13?
           ('C2','S1'): 8,  
           ('C3','S1'): 6,  
           ('C4','S1'): 11,  
           ('C5','S1'): 15,  
           ('C1','S2'): 8,  
           ('C2','S2'): 12,   ### 12?
           ('C3','S2'): 7,  
           ('C4','S2'): 6,
           ('C5','S2'): 4 }

# capacity of each facility, if opened
capacity = { 'F1': 26,
             'F2': 25,
             'F3': 16,
				 'F4': 8,
				 'F5': 10}

# fixed cost of opening each facility
fixedCost = {'F1': 120,
             'F2': 100,
             'F3': 90,
				 'F4': 40,
				 'F5': 60}


# Transportation costs per thousand units
transCost = { ('F1','C1'): 2,
              ('F1','C2'): 7,
              ('F1','C3'): 6,
              ('F1','C4'): 9,
              ('F1','C5'): 6,
              ('F2','C1'): 7,
              ('F2','C2'): 3,
              ('F2','C3'): 8,
              ('F2','C4'): 5,
              ('F2','C5'): 9,
              ('F3','C1'): 6,
              ('F3','C2'): 8,
              ('F3','C3'): 4,
              ('F3','C4'): 10,
              ('F3','C5'): 3,
				  ('F4','C1'): 2,
              ('F4','C2'): 7,
              ('F4','C3'): 6,
              ('F4','C4'): 2,
              ('F4','C5'): 5,
              ('F5','C1'): 2,
              ('F5','C2'): 9,
              ('F5','C3'): 6,
              ('F5','C4'): 7,
              ('F5','C5'): 12 }

# penalties per unit of unmet customer demand
penalties = {'C1': 30, 'C2': 30, 'C3':30, 'C4': 30, 'C5': 30 }


# Master Problem  and master decision variables
master = Model("master")

# Plant open decision variables: open[p] == 1 if plant p is open.
### Here we start with LP relaxation, so varibles are declared as continuous, with bounds [0,1]
fopen = {} 
for p in Fset:
    fopen[p] = master.addVar(ub=1.0, obj=fixedCost[p], name="Open%s" % p)

# Value function decision variables
theta = {}
for k in Sset:
    theta[k] = master.addVar(vtype=GRB.CONTINUOUS, obj=0.5, name="Theta%s" % k)


master.modelSense = GRB.MINIMIZE
master.update()


## Subproblem and subproblem decision variables
sub = Model("facility")
sub.params.logtoconsole=0  ## turns off display of Gurobi output when solving subproblems

# Transportation decision variables: how much to transport from
# a plant p to a customer w
transport = {}
for w in Cset:
    for p in Fset:
        transport[w,p] = sub.addVar(obj=transCost[p,w],
                                     name="Trans%s.%s" % (p, w))

unmet = {}
for w in Cset:
    unmet[w] = sub.addVar(obj=penalties[w],name="Unmet%s" % w)

# The objective is to minimize the total fixed and variable costs
sub.modelSense = GRB.MINIMIZE 
sub.update() 

# Subproblem production constraints
# For now just set right-hand side to capacity[p] -- it will be reset during algorithm
# based on master problem solution
capcon = {}
for p in Fset:
    capcon[p] = sub.addConstr(
        quicksum(transport[w,p] for w in Cset) <= capacity[p],
           "Capacity%s" % p)

# Demand constraints
# For now just use scenario 0 data -- it will be reset during algorithm as it loops through scnearios
demcon = {}
for w in Cset:
    demcon[w] = sub.addConstr(quicksum(transport[w,p] for p in Fset) + unmet[w] >= demand[w,'S1'],
                    "Demand%s" % w)

## Upper bound constraints (optional -- these strengthen the subproblem relaxation)
## upper bounds need to be updated to demand[w,p]*fopen[p].x
usesubprobcuts = 1
if usesubprobcuts:
    trub = {}
    for w in Cset:
        for p in Fset:
            trub[w,p] = sub.addConstr(transport[w,p] <= demand[w,'S1'])

sub.update()

# Begin the LP cutting plane loop
cutfound = 1  ## keep track if any violated cuts were found  (set cutfound = 0 to skip this loop)
iter = 1
totcuts = 0
while cutfound:

    print '================ Iteration ', iter, ' ==================='
    iter = iter+1
    # Solve current master problem
    cutfound = 0 
    master.update()
    master.optimize() 
 
    print 'current LP master objval = ', master.objVal

	 # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
    for k in Sset:
        for p in Fset:
            capcon[p].RHS = fopen[p].x*capacity[p]
        for w in Cset:
		      demcon[w].RHS = demand[w,k]
        if usesubprobcuts:
            for w in Cset:
                for p in Fset:
                    trub[w,p].RHS = demand[w,k]*fopen[p].x
        sub.update()	
        sub.optimize()

        if sub.objVal > theta[k].x + 0.000001:  ### violation tolerance
            totcuts += 1
            xcoef = {} 
            for p in Fset:
                xcoef[p] = capacity[p]*capcon[p].Pi
            rhs = 0.0 
            for w in Cset:
                rhs += demand[w,k]*demcon[w].Pi
            if usesubprobcuts:
                for p in Fset:
                    for w in Cset:
                        xcoef[p] += demand[w,k]*trub[w,p].Pi
            master.addConstr(theta[k] - quicksum(xcoef[p]*fopen[p] for p in Fset) >= rhs)
            cutfound = 1 

print 'Benders cuts in LP master: ', totcuts

### Now declare the variables to be binary
for p in Fset:
    fopen[p].vType = GRB.BINARY


### Define the callback function that Gurobi will call when it finds an integer feasible solution 
### This is where you need to search for more Benders cuts and add them if any violated
### See Gurobi's "callback.py" and "tsp.py" examples for more details


def BendersCallback(model, where):
    if where == GRB.Callback.MIPSOL:

        ## Set up and solve the Benders subproblems, just like in cutting plane loop

        # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
        for p in Fset:
            capcon[p].RHS = model.cbGetSolution(fopen[p])*capacity[p] 

        for k in Sset:
            for w in Cset:
		          demcon[w].RHS = demand[w,k]
            if usesubprobcuts:
                for w in Cset:
                    for p in Fset:
                        trub[w,p].RHS = demand[w,k]*model.cbGetSolution(fopen[p])
            sub.update()	
            sub.optimize()

            if sub.objVal > model.cbGetSolution(theta[k]) + 0.000001:  ### violation tolerance
                xcoef = {} 
                for p in Fset:
                    xcoef[p] = capacity[p]*capcon[p].Pi
                rhs = 0.0 
                for w in Cset:
                    rhs += demand[w,k]*demcon[w].Pi
                if usesubprobcuts:
                    for p in Fset:
                        for w in Cset:
                            xcoef[p] += demand[w,k]*trub[w,p].Pi
                model.cbLazy(theta[k] - quicksum(xcoef[p]*fopen[p] for p in Fset) >= rhs)


### Pass BendersCallback as argument to the optimize function on the master
master.Params.lazyConstraints = 1

### Unfortunately, 
master.optimize(BendersCallback)

print 'current optimal solution:'
for p in Fset:
    print 'fopen[', p, ']=', fopen[p].x
for k in Sset:
    print 'theta[', k, ']=', theta[k].x
print 'objval = ', master.objVal

