import cPickle
from gurobipy import *
import time;
class hw03q3:
    def __init__(self):
        self.budget = 470.0   ### hard coded here!!!

    # def read_data(self, fileName="../hw02/forest41031.pdat"):

    def read_data(self, fileName="forestrev5050775.pdat"):
        ### This file was generated from a separate python file, using the `cPickle' module
        ### This is just for convenience -- data can be read in many ways in Python

        dfile = open(fileName,'r')

        self.B = cPickle.load(dfile)  # set of base locations (list of strings)
        self.F = cPickle.load(dfile)  # set of fire districts (list of strings)
        self.S = cPickle.load(dfile)  # set of scenarios (list of strings0
        self.c = cPickle.load(dfile)  # dictionary of cost parameters c[i,j] is unit cost to relocate from i in B to j in B
        self.h = cPickle.load(dfile)  # dictionary or purchase params, h[i] is unit cost to add new resource at i in B
        self.init = cPickle.load(dfile)  # dictionary of initial resource placements, init[i] is amt of resource at i in B
        self.closesets = cPickle.load(dfile)  # dictionary of "close enough sets". closesets[f] is a list of base locations (subset of B)
        self.demscens = cPickle.load(dfile)  # dictionary of demand scnearios. demscens[s,f] is demand of resource in scenario s at district f
        self.costscens = cPickle.load(dfile) # dictionary of cost scenarios. costscnes[s,f] is cost of resource shortage in scenario s at district f
        dfile.close()

        ### NoTE: demscens and costscens are very sparse. In particular, for a given s, there are only one or two districts f
        ### that have nonzero demscens vals (and only these have costscens vals, since without demand there is no need for a cost param)
        ### Here I define the set of keys which exist in these sets

        self.SFkeys = self.demscens.keys()
        self.SFkeys = tuplelist(self.SFkeys)  ## make tuplelist for easy selection

        ### it may also be useful to have a "reverse" version of the dictionary closesets, which provides for each facility i in
        ### B, the set of districts it can serve. This is constructed here
        self.closedists = {}
        for i in self.B:
            self.closedists[i] = []

        for f in self.F:
            for i in self.closesets[f]:
                self.closedists[i].append(f)

        self.AllArcs = [(i, j) for i in self.B for j in self.B if j != i]
        self.AllArcs = tuplelist(self.AllArcs)

        # record all close relations
        self.BF = [(b, f) for b in self.B for f in self.closedists[b]]
        self.BF = tuplelist(self.BF)
        self.nscen = len(self.S)

    def MIPwBenders(self):
        '''
        Standard sequence of MIPs with multi-cut Benders Decomposition
        :return:
        '''

        ''' Master problem's variable '''
        master = Model("master")
        master.params.logtoconsole = 0

        moveAmount_mp = {}
        for arc in self.AllArcs:
            moveAmount_mp[arc] = master.addVar(obj=0.0, vtype=GRB.INTEGER, name='move_%s_%s' % (arc[0], arc[1]))

        buyAmount_mp = {}
        nunit_mp = {}
        left_mp = {}
        for b in self.B:
            buyAmount_mp[b] = master.addVar(obj=0.0, vtype=GRB.INTEGER, name='purchase_%s' % (b))
            nunit_mp[b] = master.addVar(obj=0.0, name='nunit_%s' % (b))
            left_mp[b]  = master.addVar(name='left_%s' % (b))

        theta = {}
        for k in self.S:
            theta[k] = master.addVar(obj=1.0/self.nscen, name="Theta_%s" %k)


        master.modelSense = GRB.MINIMIZE
        master.update()

        ''' Master problem's constraints '''
        # Budget constraints (sum_b {h_b x_b} + sum_arc c_ij z_ij <= C)
        budgetcon = master.addConstr(
            quicksum(buyAmount_mp[b] * self.h[b] for b in self.B) + quicksum(moveAmount_mp[arc] * self.c[arc]
                                                            for arc in self.AllArcs) <= self.budget, name='Budget')

        capacitycon = {}
        for b in self.B:
            capacitycon[b] = master.addConstr(
                buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b)) + left_mp[b] == nunit_mp[b],
                name='Capacity_%s' % (b))

        for b in self.B:
            master.addConstr(self.init[b] - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')) == left_mp[b], name='Capacity2_%s' %(b))

        master.update()

        ''' Subproblem's variables '''
        sub = Model("subproblem")
        sub.params.logtoconsole = 0

        unmet_sub = {}
        for f in self.F:
            unmet_sub[f] = sub.addVar(obj=0, name='Unmet_%s' % (f))

        respond_sub = {}
        for b, f in self.BF:
            respond_sub[b, f] = sub.addVar(obj=0, name="Respond_%s_%s" % (b, f))

        sub.modelSense = GRB.MINIMIZE
        sub.update()

        ''' Subproblem's constraints '''
        # Capacity constraints (sum_f {y_bfk} <= V_b)
        capcon = {}
        for b in self.B:
            capcon[b] = sub.addConstr(
                quicksum(respond_sub[arc] for arc in self.BF.select(b, '*')) <= 0, name='Capacity_%s' % (b))

        # Demand constraints (sum_b {y_bfk} + u_fk = d_fk)
        demcon = {}
        for f in self.F:
            demcon[f] = sub.addConstr(
                quicksum(respond_sub[arc] for arc in self.BF.select('*', f)) + unmet_sub[f] == 0, name='Demand_%s' % (f))

        sub.update()

        ''' Begin the cutting plane loop '''
        startTime = time.time()
        cutfound = 1
        iter = 1
        nCuts = 0
        print('\n-------------------------- Bender Decomposition ------------------------------')
        while cutfound:
            # Solve current master problem
            cutfound = 0
            master.update()
            master.optimize()

            assert master.status == GRB.Status.OPTIMAL
            lb = master.objVal
            ub = float("inf")
            # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
            for b in self.B:
                capcon[b].RHS = nunit_mp[b].x

            Q = {}
            for k in self.S:
                rhs = 0.0
                for f in self.F:
                    if (k, f) in self.SFkeys:
                        unmet_sub[f].obj = self.costscens[k, f]
                        demcon[f].RHS = self.demscens[k, f]
                    else:
                        unmet_sub[f].obj = 0.0
                        demcon[f].RHS = 0.0
                sub.update()
                sub.optimize()

                assert sub.status == GRB.Status.OPTIMAL
                Q[k] = sub.objVal

                # fast implementation with quicksum rather than for-loop
                rhs += quicksum(
                    capcon[b].Pi * nunit_mp[b] for b in self.B)

                # rhs += quicksum(self.demscens[k, f] * demcon[f].Pi for (k, f) in self.SFkeys) # BUG!
                rhs += quicksum(self.demscens[k, f] * demcon[f].Pi for (k, f) in self.SFkeys.select(k,'*'))

                if Q[k] > theta[k].x + 0.000001: ### violation tolerance
                    master.addConstr(theta[k] >= rhs)
                    nCuts += 1
                    cutfound = 1

            assert abs(lb - quicksum(theta[k].x for k in self.S).getValue() / float(self.nscen)) < 0.00000001
            # update upperbound when finding a feasible MIP solution
            ub = quicksum(Q[k] for k in self.S).getValue() / float(self.nscen)
            print("Iter %d: [lowerBound, upperBound] = [%f, %f],  nCut = %d" % (iter, lb, ub, nCuts))
            iter = iter + 1

        timeCost = time.time() - startTime

        print 'current optimal solution:'
        for b in self.B:
            if buyAmount_mp[b].x:
                print 'X[', b, '] = ', buyAmount_mp[b].x
        for arc in self.AllArcs:
            if moveAmount_mp[arc].x:
                print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
        print 'EXPECTED COST: ', master.objVal
        print 'TIME COST: ', timeCost


    def BranchAndCut(self, cutfound):
        '''
        Branch-and-Cut with Benders cuts added within a lazy constraint callback
        :param cutfound: flag to indicate whether to use phase 0 to keep the Benders cuts for the 1-st stage LP relaxation.
        :return:
        '''
        ''' Master problem's variable '''
        master = Model("master")
        # master.params.logtoconsole = 0

        moveAmount_mp = {}
        for arc in self.AllArcs:
            moveAmount_mp[arc] = master.addVar(name='move_%s_%s' % (arc[0], arc[1]))

        buyAmount_mp = {}
        nunit_mp = {}
        for b in self.B:
            buyAmount_mp[b] = master.addVar(name='purchase_%s' % (b))
            nunit_mp[b] = master.addVar(name='nunit_%s' % (b))

        theta = {}
        for k in self.S:
            theta[k] = master.addVar(obj=1.0/self.nscen, name="Theta_%s" %k)

        master.modelSense = GRB.MINIMIZE
        master.update()

        ''' Master problem's constraints '''
        # Budget constraints (sum_b {h_b x_b} + sum_arc c_ij z_ij <= C)
        budgetcon = master.addConstr(
            quicksum(buyAmount_mp[b] * self.h[b] for b in self.B) + quicksum(moveAmount_mp[arc] * self.c[arc]
                                                            for arc in self.AllArcs) <= self.budget, name='Budget')

        capacitycon = {}
        for b in self.B:
            capacitycon[b] = master.addConstr(
                self.init[b] + buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
                - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')) == nunit_mp[b],
                name='Capacity_%s' % (b))

        for b in self.B:
            master.addConstr(self.init[b] - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')) >= 0, name='Capacity2_%s' %(b))

        master.update()

        ''' Subproblem's variables '''
        sub = Model("subproblem")
        sub.params.logtoconsole = 0

        unmet_sub = {}
        for f in self.F:
            unmet_sub[f] = sub.addVar(name='Unmet_%s' % (f))

        respond_sub = {}
        for b, f in self.BF:
            respond_sub[b, f] = sub.addVar(name="Respond_%s_%s" % (b, f))

        sub.modelSense = GRB.MINIMIZE
        sub.update()

        ''' Subproblem's constraints '''
        # Capacity constraints (sum_f {y_bfk} <= V_b)
        capcon = {}
        for b in self.B:
            capcon[b] = sub.addConstr(
                quicksum(respond_sub[arc] for arc in self.BF.select(b, '*')) <= 0, name='Capacity_%s' % (b))

        # Demand constraints (sum_b {y_bfk} + u_fk = d_fk)
        demcon = {}
        for f in self.F:
            demcon[f] = sub.addConstr(
                quicksum(respond_sub[arc] for arc in self.BF.select('*', f)) + unmet_sub[f] == 0, name='Demand_%s' % (f))

        sub.update()

        ''' Begin the cutting plane loop '''
        startTime = time.time()
        iter = 1
        nCuts = 0
        print('\n-------------------------- Branch-and-Cut ------------------------------')
        while cutfound:
            # Solve current master problem
            cutfound = False
            master.update()
            master.optimize()

            assert master.status == GRB.Status.OPTIMAL
            lb = master.objVal
            # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
            for b in self.B:
                capcon[b].RHS = nunit_mp[b].x

            Q = {}
            for k in self.S:
                rhs = 0.0
                for f in self.F:
                    if (k, f) in self.SFkeys:
                        unmet_sub[f].obj = self.costscens[k, f]
                        demcon[f].RHS = self.demscens[k, f]
                    else:
                        unmet_sub[f].obj = 0
                        demcon[f].RHS = 0
                sub.update()
                sub.optimize()

                assert sub.status == GRB.Status.OPTIMAL
                Q[k] = sub.objVal

                rhs += quicksum(
                    capcon[b].Pi * (self.init[b] + buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
                                - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*'))) for b in self.B)

                rhs += quicksum(self.demscens[k, f] * demcon[f].Pi for (k, f) in self.SFkeys.select(k,'*'))

                if Q[k] > theta[k].x + 0.000001: ### violation tolerance
                    master.addConstr(theta[k] >= rhs)
                    nCuts += 1
                    cutfound = True

            assert abs(lb - quicksum(theta[k].x for k in self.S).getValue() / float(self.nscen)) < 0.00000001
            print("Iter %d: lowerBound = %f,  nCut = %d" % (iter, lb, nCuts))
            iter = iter + 1

        ### Now declare the variables to be integral
        for arc in self.AllArcs:
            moveAmount_mp[arc].vType = GRB.INTEGER
        for b in self.B:
            buyAmount_mp[b].vType = GRB.INTEGER
        master._vars = sub
        master.update()

        def BendersCallback(model, where):
            if where == GRB.Callback.MIPSOL:
                for b in self.B:
                    capcon[b].RHS = model.cbGetSolution(nunit_mp[b])

                for k in self.S:
                    rhs = 0.0
                    for f in self.F:
                        if (k, f) in self.SFkeys:
                            unmet_sub[f].obj = self.costscens[k, f]
                            demcon[f].RHS = self.demscens[k, f]
                        else:
                            unmet_sub[f].obj = 0
                            demcon[f].RHS = 0
                    sub = model._vars
                    sub.update()
                    sub.optimize()

                    assert sub.status == GRB.Status.OPTIMAL

                    rhs += quicksum(
                        capcon[b].Pi * (self.init[b] + buyAmount_mp[b] + quicksum(
                            moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
                                        - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*'))) for b in
                        self.B)

                    rhs += quicksum(self.demscens[k, f] * demcon[f].Pi for (k, f) in self.SFkeys.select(k,'*'))

                    if sub.objVal > model.cbGetSolution(theta[k]) + 0.000001:  ### violation tolerance
                        model.cbLazy(theta[k] >= rhs)

        master.Params.lazyConstraints = 1
        master.optimize(BendersCallback)

        timeCost = time.time() - startTime
        print 'current optimal solution:'
        for b in self.B:
            if buyAmount_mp[b].x > 0:
                print 'X[', b, '] = ', buyAmount_mp[b].x
        for arc in self.AllArcs:
            if moveAmount_mp[arc].x > 0:
                print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
        print('EXPECTED COST: ', master.objVal)
        print('TIME COST: ', timeCost)


if __name__ == "__main__":
    mysolver = hw03q3()
    mysolver.read_data()
    mysolver.MIPwBenders()
    mysolver.BranchAndCut(False)
    mysolver.BranchAndCut(True)