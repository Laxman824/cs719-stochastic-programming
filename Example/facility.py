from gurobipy import *

Cset = ['C1','C2','C3','C4']
Sset = ['S1','S2']
Fset = ['F1','F2','F3']

#demand values in each scenario
demand = { ('C1','S1'): 13,    ### 13?
           ('C2','S1'): 8,  
           ('C3','S1'): 6,  
           ('C4','S1'): 11,  
           ('C1','S2'): 8,  
           ('C2','S2'): 12,   ### 12?
           ('C3','S2'): 7,  
           ('C4','S2'): 6 }

# capacity of each facility, if opened
capacity = { 'F1': 26,
             'F2': 25,
             'F3': 18 }

# fixed cost of opening each facility
fixedCost = {'F1': 120,
             'F2': 100,
             'F3': 90 }


# Transportation costs per thousand units
transCost = { ('F1','C1'): 2,
              ('F1','C2'): 7,
              ('F1','C3'): 6,
              ('F1','C4'): 9,
              ('F2','C1'): 7,
              ('F2','C2'): 3,
              ('F2','C3'): 8,
              ('F2','C4'): 5,
              ('F3','C1'): 6,
              ('F3','C2'): 8,
              ('F3','C3'): 4,
              ('F3','C4'): 10 }

# penalties per unit of unmet customer demand
penalties = {'C1': 30, 'C2': 30, 'C3':30, 'C4': 30 }


# Master Problem  and master decision variables
master = Model("master")

# Plant open decision variables: open[p] == 1 if plant p is open.
fopen = {} 
for p in Fset:
    fopen[p] = master.addVar(vtype=GRB.BINARY, obj=fixedCost[p], name="Open%s" % p)

# Value function decision variables
theta = {}
for k in Sset:
    theta[k] = master.addVar(vtype=GRB.CONTINUOUS, obj=0.5, name="Theta%s" % k)


master.modelSense = GRB.MINIMIZE
master.update()


## Subproblem and subproblem decision variables
sub = Model("facility")

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
    demcon[w] = sub.addConstr(quicksum(transport[w,p] for p in Fset) + unmet[w] == demand[w,'S1'],
                    "Demand%s" % w)

sub.update()

# Begin the cutting plane loop
cutfound = 1  ## keep track if any violated cuts were found
iter = 1
while cutfound:

    print '================ Iteration ', iter, ' ==================='
    iter = iter+1
    # Solve current master problem
    cutfound = 0 
    master.update()
    master.optimize() 
 
    print 'current optimal solution:'
    for p in Fset:
        print 'fopen[', p, ']=', fopen[p].x
    for k in Sset:
        print 'theta[', k, ']=', theta[k].x
    print 'objval = ', master.objVal

	 # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
    for k in Sset:
        for p in Fset:
            capcon[p].RHS = fopen[p].x*capacity[p]
        for w in Cset:
		      demcon[w].RHS = demand[w,k]
        sub.update()	
        sub.optimize()

        # Display info, compute Benders cut, display, add to master 
        print 'sub[', k, '] objval = ', sub.objVal
        for w in Cset:
            print 'unmet[', w, ']=', unmet[w].x
            #for p in facil:
            #    print transport[w][p].x
        if sub.objVal > theta[k].x + 0.000001:  ### violation tolerance
            xcoef = {} 
            for p in Fset:
                xcoef[p] = capacity[p]*capcon[p].Pi
                print 'xcoef[',p,']=', xcoef[p]
                print
            rhs = 0.0 
            for w in Cset:
                rhs += demand[w,k]*demcon[w].Pi
            print 'rhs = ', rhs
            master.addConstr(theta[k] - quicksum(xcoef[p]*fopen[p] for p in Fset) >= rhs)
            cutfound = 1 


