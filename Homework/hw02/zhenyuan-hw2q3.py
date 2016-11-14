import cPickle
from gurobipy import *

class hw02q3:
    def __init__(self):
        pass
    def read_data(self, fileName = 'forest41031.pdat'):
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
        # self.budget = 500.0   ### hard coded here!!!
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

        ### This is just a check of the data. Probably you want to comment/delete these lines once you see the structure

        print self.B
        print self.F
        print self.S
        print self.c
        print self.h
        print self.init
        print self.closesets
        print self.demscens
        print self.costscens


        self.demscensDict = {}
        for (i,j) in self.demscens:
            try:
                self.demscensDict[i].append(j)
            except KeyError:
                self.demscensDict[i] = [j]
        for i in self.demscensDict:
            print i, self.demscensDict[i]

        self.AllArcs = [(i,j) for i in self.B for j in self.B if j != i]
        self.AllArcs = tuplelist(self.AllArcs)

        # record all close relations
        self.BF = [(b,f) for b in self.B for f in self.closedists[b]]
        self.BF = tuplelist(self.BF)

        self.nscen = len(self.S)

    def runSingleCutBender(self, budget = 500.0):
        '''
        Apply single-cut bender decomposition
        :param budget: total budget (500.0 by default)
        :return:
        '''
        ''' Master problem's variable '''
        master = Model("master")
        master.params.logtoconsole=0

        moveAmount_mp = {}
        for arc in self.AllArcs:
            moveAmount_mp[arc] = master.addVar(obj=0.0, name='move_%s_%s' %(arc[0], arc[1]))

        buyAmount_mp = {}
        nunit_mp = {}
        for b in self.B:
            buyAmount_mp[b] = master.addVar(obj=0.0, name='purchase_%s' %(b))
            nunit_mp[b] = master.addVar(obj=0.0, name='nunit_%s' % (b))

        theta = master.addVar(vtype=GRB.CONTINUOUS, obj=1.0, name="Theta")


        master.modelSense = GRB.MINIMIZE
        master.update()

        ''' Master problem's constraints '''
        # Budget constraints (sum_b {h_b x_b} + sum_arc c_ij z_ij = C)
        budgetcon = master.addConstr(
            quicksum(buyAmount_mp[b]*self.h[b] for b in self.B) + quicksum(moveAmount_mp[arc]*self.c[arc]
                                        for arc in self.AllArcs) <= budget, name='Budget')

        capacitycon = {}
        for b in self.B:
            capacitycon[b] = master.addConstr(self.init[b] + buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
             - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')) == nunit_mp[b], name='Capacity_%s' %(b))

        master.update()

        ''' Subproblem's variables '''
        sub = Model("subproblem")
        sub.params.logtoconsole=0

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
        cutfound = 1  ## keep track if any violated cuts were found
        iter = 1
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
            # for b in self.B:
            #     print 'X[' ,b ,'] = ', buyAmount_mp[b].x
            # for arc in self.AllArcs:
            #     print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
            # print 'objval = ', master.objVal
            lb = master.objVal

             # Fix the right-hand side in subproblem constraints according to each scenario and master solution, then solve
            for b in self.B:
                capcon[b].RHS = nunit_mp[b].x

            Q = 0
            rhs = 0.0
            xzcoef = {}
            for k in self.S:
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

                for b in self.B:
                    xzcoef[k, b] = capcon[b].Pi
                    rhs += xzcoef[k,b] * self.init[b]

                for f in self.F:
                    if (k, f) in self.SFkeys:
                        rhs += self.demscens[k,f] * demcon[f].Pi

            assert lb == theta.x

            ub = Q/float(self.nscen)
            print("    [lowerBound, upperBound] = [%f, %f],  nCut = %d" % (lb, ub, nCuts))
            if ub > lb + 0.000001:  ### violation tolerance
                master.addConstr(
                    theta - 1.0/self.nscen * quicksum(
                        xzcoef[k, b] * (buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
                                        - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')))
                                     for (k,b) in xzcoef) >= rhs/float(self.nscen))
                nCuts += 1
                cutfound = 1



        print('\nEXPECTED COST : %g' % master.objVal)
        print 'current optimal solution:'
        for b in self.B:
            print 'X[', b, '] = ', buyAmount_mp[b].x
        for arc in self.AllArcs:
            print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
        print 'objval = ', master.objVal


    def runLevelMethod(self, lamb = 0.2929, budget = 500.0):
        '''
        Apply level method
        :param lamb: lambda (0.2929 by default)
        :return:
        '''

        nunit = {}
        for b in self.B:
            nunit[b] = self.init[b]

        ''' Master problem's variable '''
        master = Model("master")
        master.params.logtoconsole=0

        moveAmount_mp = {}
        for arc in self.AllArcs:
            moveAmount_mp[arc] = master.addVar(obj=0.0, name='move_%s_%s' %(arc[0], arc[1]))

        buyAmount_mp = {}
        nunit_mp = {}
        for b in self.B:
            buyAmount_mp[b] = master.addVar(obj=0.0, name='purchase_%s' %(b))
            nunit_mp[b] = master.addVar(obj=0.0, name='nunit_%s' % (b))

        theta = master.addVar(vtype=GRB.CONTINUOUS, obj=1.0, name="Theta")

        phi = master.addVar(vtype=GRB.CONTINUOUS, obj=0.0, name="Phi")

        master.modelSense = GRB.MINIMIZE
        master.update()

        ''' Master problem's constraints '''
        # Budget constraints (sum_b {h_b x_b} + sum_arc c_ij z_ij = C)
        budgetcon = master.addConstr(
            quicksum(buyAmount_mp[b]*self.h[b] for b in self.B) + quicksum(moveAmount_mp[arc]*self.c[arc]
                                        for arc in self.AllArcs) <= budget, name='Budget')

        capacitycon = {}
        for b in self.B:
            capacitycon[b] = master.addConstr(self.init[b] + buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
             - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')) == nunit_mp[b], name='Capacity_%s' %(b))

        master.update()

        ''' Subproblem's variables '''
        sub = Model("subproblem")
        sub.params.logtoconsole=0

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
        iter = 1
        fxt = []
        print('\n-------------------------- Level Method ------------------------------')
        while 1:
            print '--- Iteration ', iter, ' ---'
            iter = iter + 1
            for b in self.B:
                capcon[b].RHS = nunit[b]  # nunit[b] is the # units at location b for t-th iteration, update iteratively.
            Q = 0
            rhs = 0.0
            xzcoef = {}
            for k in self.S:
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

                for b in self.B:
                    xzcoef[k, b] = capcon[b].Pi
                    rhs += xzcoef[k,b] * self.init[b]

                for f in self.F:
                    if (k, f) in self.SFkeys:
                        rhs += self.demscens[k,f] * demcon[f].Pi

            fxt.append(Q/float(self.nscen))

            ## compute z_lb
            # add one cut to m_t: \theta >= \pi_cap * v + \pi_dem * d; note that v = init + x + z_in - z_out (v is intermediate var)
            master.addConstr(
                theta - 1.0 / self.nscen * quicksum(
                    xzcoef[k, b] * (
                    buyAmount_mp[b] + quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select('*', b))
                    - quicksum(moveAmount_mp[arc] for arc in self.AllArcs.select(b, '*')))
                    for (k, b) in xzcoef) >= rhs / float(self.nscen))
            master.update()
            master.optimize()

            assert master.status == GRB.Status.OPTIMAL
            # print 'current optimal solution:'
            # for b in self.B:
            #     print 'X[' ,b ,'] = ', buyAmount_mp[b].x
            # for arc in self.AllArcs:
            #     print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
            # print 'objval = ', master.objVal
            zlb = master.objVal

            assert zlb == theta.x

            ## compute z_ub
            zub = min(fxt)
            print("    [lowerBound, upperBound] = [%f, %f]" % (zlb, zub))
            if zub > zlb + 0.000001:
                lt = zlb + (zub-zlb) * lamb

                # compute x_{t+1}
                theta.obj = 0.0
                phi.obj = 1.0
                tempUb = master.addConstr(theta <= lt)
                tempDist = master.addQConstr(quicksum((buyAmount_mp[b] - buyAmount_mp[b].x)*(buyAmount_mp[b] - buyAmount_mp[b].x) for b in self.B)
                                            + quicksum((moveAmount_mp[arc] - moveAmount_mp[arc])*(moveAmount_mp[arc] - moveAmount_mp[arc]) for arc in self.AllArcs)
                                            <= phi)
                master.update()
                master.optimize()
                assert master.status == GRB.Status.OPTIMAL

                # update for next iteration t+1
                for b in self.B:
                    nunit[b] = nunit_mp[b].x

                # reset the constraints and variables
                theta.obj = 1.0
                phi.obj = 0.0
                master.remove(tempUb)
                master.remove(tempDist)
            else:
                break;



        print('\nEXPECTED COST : %g' % master.objVal)
        print 'current optimal solution:'
        for b in self.B:
            print 'X[', b, '] = ', buyAmount_mp[b].x
        for arc in self.AllArcs:
            print 'Z[', arc[0], ' ', arc[1], ']=', moveAmount_mp[arc].x
        print 'objval = ', master.objVal


if __name__ == "__main__":
    mysolver = hw02q3()
    mysolver.read_data()
    mysolver.runSingleCutBender()
    mysolver.runLevelMethod()
