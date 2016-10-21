import cPickle
from gurobipy import *
import time
import numpy as np
import math

class hw02q4:
    def __init__(self):
        self.M = 15 # number of batches
        self.n = 50 # number of scenarios
        self.n1 = 50 # number of sample size
        self.expUB = 200 # upper bound of expansion(to make the 1st stage feasible region compact)
        self.N = 1000 # number of samples for estimation
        self.tvalue_14 = 2.145
        self.tvalue_1000 = 1.962 # should use tvalue_999 however.

    def read_data(self, fileName="nd848.pdat"):
        ### Read data from file you choose: commont/uncomment to choose the different files
        ### This file was generated from a separate python file, using the `cPickle' module
        ### This is just for convenience -- data can be read in many ways in Python

        dfile = open(fileName,'r')

        self.Fset = cPickle.load(dfile)  # set of facilities (list of strings)
        self.Hset = cPickle.load(dfile)  # set of warehouses (list of strings)
        self.Cset = cPickle.load(dfile)  # set of customers (list of strings)
        self.arcExpCost = cPickle.load(dfile)  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
        self.facCap = cPickle.load(dfile)   # facility capacities (dictionary mapping F to floats)
        self.curArcCap = cPickle.load(dfile)  # current arc capacities (dictionary mapping (i,j) to floats, where either
                                         # i is facility, j is warehouse, or i is warehouse and j is customer
        self.unmetCost = cPickle.load(dfile)  # penalty for unment customer demands (dicationary mapping C to floats)
        self.demmean = cPickle.load(dfile)  #  mean demand (conditional that it is positive)
        self.demstdev = cPickle.load(dfile)   # stdev of demand (conditional that it is positive)
        self.probnodem = cPickle.load(dfile)  # probability each demand is zero
        dfile.close()


        ### Define sets of arcs (used as keys to dictionaries)
        self.FHArcs = [(i,j) for i in self.Fset for j in self.Hset]  ## arcs from facilities to warehouses
        self.HCArcs = [(i,j) for i in self.Hset for j in self.Cset]   ## arcs from warehouses to customers
        self.AllArcs = self.FHArcs + self.HCArcs
        print self.AllArcs

        ### Make them Gurobi tuplelists
        self.FHArcs = tuplelist(self.FHArcs)
        self.HCArcs = tuplelist(self.HCArcs)
        self.AllArcs = tuplelist(self.AllArcs)


    def solveEF(self, sampledDemScens, nscen):
        '''
        function for building and solving the stochastic programming via EF formulation
        :param sampledDemScens: n sets of sampled demands of customers (n is the # scenarios)
        :return: [total cost(stage 1 + stage 2), optimal solution]
        '''
        m = Model("ArcExpansionExtensiveForm")
        m.params.logtoconsole=0

        ''' Variables '''
        # First stage vars, arc expansion amounts
        arcExpAmount = {}
        for arc in self.AllArcs:
            arcExpAmount[arc] = m.addVar(obj=self.arcExpCost[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

        # Second stage vars, unmet demands (z_jk)
        unmet = {}
        for c in self.Cset:
            for k in range(nscen):
                unmet[c, k] = m.addVar(obj=float(self.unmetCost[c]) / nscen, name='Unmet_%s_%s' % (c, k))

        # Intermediate vars, transport
        transport = {}
        for arc in self.AllArcs:
            for k in range(nscen):
                transport[arc[0], arc[1], k] = m.addVar(obj=0, name="Transport_%s_%s_%s" % (arc[0], arc[1], k))

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' Constraints '''
        # Expansion upper bound (to make the 1st stage feasible region compact)
        for arc in self.AllArcs:
            m.addConstr(arcExpAmount[arc] <= self.expUB, name='ExpansionLimit_%s_%s' % (arc[0], arc[1]))

        # Production constraints (sum_h {x_fhk} - b_f <= 0)
        for f in self.Fset:
            for k in range(nscen):
                m.addConstr(
                    quicksum(transport[f, h, k] for h in self.Hset) <= self.facCap[f], name='Capacity_%s' % (f))

        # Demand constraints (sum_i {y_ijk} + z_jk = d_jk)
        for c in self.Cset:
            for k in range(nscen):
                m.addConstr(
                    quicksum(transport[h, c, k] for h in self.Hset) + unmet[c, k] >= sampledDemScens[c][k],
                    name='Demand_%s_%s' % (c, k))

        # Line capacity constraints
        for arc in self.AllArcs:
            for k in range(nscen):
                m.addConstr(transport[arc[0], arc[1], k] <= self.curArcCap[arc] + arcExpAmount[arc],
                            name='LineCapacity_%s_%s_%s' % (arc[0], arc[1], k))

        # Transportation contraints
        for h in self.Hset:
            for k in range(nscen):
                m.addConstr(
                    # quicksum(transport[f,h,k] for f in Fset) >= quicksum(transport[h,c,k] for c in Cset),
                    # name="Transportation_%s_%s" % (h, k))
                    quicksum(transport[f, h, k] for f, h in self.FHArcs.select('*', h)) >= quicksum(
                        transport[h, c, k] for h, c in self.HCArcs.select(h, '*')),
                    name="Transportation_%s_%s" % (h, k))

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL
        candidate = {}
        for arc in self.AllArcs:
            candidate[arc] = arcExpAmount[arc].x
        return [m.objval, candidate]

    def solveSSP(self, arcExp, sampledDemScens):
        '''
        function for building and solving the second stage problem given first stage solution
        :param arcExp: first stage solution
        :param sampledDemScens: ONE set of sampled demands of customers
        :return: second stage cost
        '''
        m = Model("ArcExpansionSecondStage")
        m.params.logtoconsole=0

        ''' Variables '''
        # Second stage vars, unmet demands (z_jk)
        unmet = {}
        for c in self.Cset:
            unmet[c] = m.addVar(obj=self.unmetCost[c], name='Unmet_%s' % (c))

        # Intermediate vars, transport
        transport = {}
        for arc in self.AllArcs:
            transport[arc] = m.addVar(obj=0, name="Transport_%s_%s" % (arc[0], arc[1]))

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' Constraints '''
        # Production constraints (sum_h {x_fhk} - b_f <= 0)
        for f in self.Fset:
            m.addConstr(
                quicksum(transport[f, h] for h in self.Hset) <= self.facCap[f], name='Capacity_%s' % (f))

        # Demand constraints (sum_i {y_ijk} + z_jk = d_jk)
        for c in self.Cset:
            m.addConstr(
                quicksum(transport[h, c] for h in self.Hset) + unmet[c] >= sampledDemScens[c], name='Demand_%s' % (c))
                # !! WARN: sampledDemScens[c] maybe negative, so use '>=' not '=='(probably results in infeasibility)

        # Line capacity constraints
        for arc in self.AllArcs:
            m.addConstr(transport[arc] <= self.curArcCap[arc] + arcExp[arc], name='LineCapacity_%s_%s' % (arc[0], arc[1]))

        # Transportation contraints
        for h in self.Hset:
            m.addConstr(
                # quicksum(transport[f,h] for f in Fset) >= quicksum(transport[h,c] for c in Cset),
                # name="Transportation_%s_%s" % (h))
                quicksum(transport[f, h] for f, h in self.FHArcs.select('*', h)) >= quicksum(
                    transport[h, c] for h, c in self.HCArcs.select(h, '*')), name="Transportation_%s" % (h))

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL
        return m.objval

    def saa_1(self):
        '''
        First separately estimate lower and upper bounds
        :return:
        '''
        objvals = np.zeros([self.M, 1])
        arcExpansion = {}
        # compute lower bound E[Zn^*]
        for k in range(self.M):
            samples = {}
            ### set up and solve extensive form for a sample average approximation problem
            for c in self.Cset:
                identity = np.random.rand(self.n, 1)
                samp = self.demstdev[c] * np.random.randn(self.n, 1) + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                samples[c] = samp

            [objvals[k], arcExpansion[k]] = self.solveEF(samples, self.n)

        # print 'mean SAA objval = ', np.mean(objvals)
        # print 'stdev of SAA objval = ', np.std(objvals)
        lbmean = np.mean(objvals)
        lbwidth = np.std(objvals) / math.sqrt(self.M) * self.tvalue_14

        # compute upper bound Ef(\hat{x}, \ksi)
        candArcExp = arcExpansion[objvals.argmin()]
        evalvals = np.zeros([self.N, 1])
        for k in range(self.N):
            samples = {}
            for c in self.Cset:
                identity = np.random.rand()
                samp = self.demstdev[c] * np.random.randn() + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                samples[c] = samp
            evalvals[k] = quicksum(candArcExp[arc]*self.arcExpCost[arc] for arc in self.AllArcs).getValue() \
                          + self.solveSSP(candArcExp, samples)

        # print 'mean solution objval = ', np.mean(evalvals)
        # print 'stdev of solution objval = ', np.std(evalvals)
        ubmean = np.mean(evalvals)
        ubwidth = np.std(evalvals)/math.sqrt(self.N)*self.tvalue_1000
        print('\n-------------------------- separately estimate lower and upper bounds ------------------------------')
        print 'ci on lower bound = [', lbmean - lbwidth, ',', lbmean + lbwidth, ']'
        print 'ci on upper bound = [', ubmean - ubwidth, ',', ubmean + ubwidth, ']'
        print 'The width of the confidence intervial is ', (ubmean + ubwidth) - (lbmean - lbwidth)

    def saa_2(self):
        '''
        Now let's try the method which directly estimates optimality gap
        :return:
        '''

        # Generate the candidate solution by solving a SAA problem with n = 50 scenarios.
        samples = {}
        for c in self.Cset:
            identity = np.random.rand(self.n, 1)
            samp = self.demstdev[c] * np.random.randn(self.n, 1) + self.demmean[c]
            samp = samp * (identity > self.probnodem[c])
            samples[c] = samp
        [objvals, candArcExp] = self.solveEF(samples, self.n)

        # Then use M = 15 batches of size n1 = 50 to estimate the optimality gap of the solution.
        gapvals = np.zeros(self.M)
        for k in range(self.M):
            samples = {}
            for c in self.Cset:
                identity = np.random.rand(self.n1, 1)
                samp = self.demstdev[c] * np.random.randn(self.n1, 1) + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                samples[c] = samp

            [sampleopt, curArcExp] = self.solveEF(samples, self.n1)

            # now evaluate the candidate solution using the same scenarios
            curvals = np.zeros(self.n1)
            for j in range(self.n1):
                sample = {}
                for c in samples:
                    sample[c] = samples[c][j]
                curvals[j] = quicksum(candArcExp[arc]*self.arcExpCost[arc] for arc in self.AllArcs).getValue()\
                             + self.solveSSP(candArcExp, sample) # !! SHOULD include 1-st stage cost

            gapvals[k] = np.mean(curvals) - sampleopt
        print('\n-------------------------- directly estimates optimality gap ------------------------------')
        print 'mean gap estimate = ', np.mean(gapvals)
        print '95% c.i. on gap = [0, ', np.mean(gapvals) + np.std(gapvals)/math.sqrt(self.M)*self.tvalue_14, ']'

        # Also estimate the objective value of the candidate solution using an independent sample of size N = 1000.
        evalvals = np.zeros([self.N, 1])
        for k in range(self.N):
            samples = {}
            for c in self.Cset:
                identity = np.random.rand()
                samp = self.demstdev[c] * np.random.randn() + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                samples[c] = samp
            evalvals[k] = quicksum(candArcExp[arc]*self.arcExpCost[arc] for arc in self.AllArcs).getValue() \
                          + self.solveSSP(candArcExp, samples)
        ubmean = np.mean(evalvals)
        ubwidth = np.std(evalvals) / math.sqrt(self.N) * self.tvalue_1000
        print 'ci on upper bound = [', ubmean - ubwidth, ',', ubmean + ubwidth, ']'

if __name__ == "__main__":
    mysolver = hw02q4()
    mysolver.read_data()
    mysolver.saa_1()
    mysolver.saa_2()