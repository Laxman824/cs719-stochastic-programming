import cPickle
from gurobipy import *
import time
import numpy as np
import math

class hw03q1:
    def __init__(self):
        self.T = 3000 # number of iterations
        self.expUB = 200.0 # upper bound of expansion(to make the 1st stage feasible region compact)
        self.N = 1000 # number of samples for estimation
        self.M = 15 # number of batches
        self.tvalue_14 = 2.145
        self.tvalue_14_oneside = 1.761
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

    def SA(self):
        ## initialize 1-st stage decision var
        arcExp = {arc: self.expUB*np.random.rand() for arc in self.AllArcs}  # candArcExp := X
        t = 0
        while t < self.T:
            ## compute subgradient by solving 2-nd stage SP
            t += 1
            sample = {}
            for c in self.Cset:
                identity = np.random.rand()
                samp = self.demstdev[c] * np.random.randn() + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                sample[c] = samp
            # compute subgradient of Q(x,\xi)
            Q, subQ = self.solveSSP(arcExp, sample)

            objVal = quicksum(self.arcExpCost[arc] * arcExp[arc] for arc in self.AllArcs).getValue() + Q

            gamma_t = 1.0 / (t ** 0.75)
            y = {arc: arcExp[arc] - gamma_t * (self.arcExpCost[arc] + subQ[arc]) for arc in self.AllArcs}
            arcExp = self.project(y)

            # print("Objective value at t={0} is {1}".format(t, objVal))

        ## 1) extimate 95% C.I. with 1000 scenarios on the objective value (also is the upperbound of the SA) Ef(\hat{x}, \ksi)
        evalvals = np.zeros([self.N, 1])

        samples = []
        coefx = {arc: 0.0 for arc in self.AllArcs}
        for k in range(self.N):
            sample = {}
            for c in self.Cset:
                identity = np.random.rand()
                samp = self.demstdev[c] * np.random.randn() + self.demmean[c]
                samp = samp * (identity > self.probnodem[c])
                sample[c] = samp
            samples.append(sample)
            Q, subQ = self.solveSSP(arcExp, sample)
            assert max(subQ.values()) <= 0

            # compute f(\hat{x}, \ksi)
            evalvals[k] = quicksum(self.arcExpCost[arc] * arcExp[arc] for arc in self.AllArcs).getValue() + Q

            for arc in self.AllArcs:
                coefx[arc] += (self.arcExpCost[arc] + subQ[arc])/float(self.N)

        ubmean = np.mean(evalvals)
        ubwidth = np.std(evalvals) / math.sqrt(self.N) * self.tvalue_1000
        print '95% c.i. on the objective value = [', ubmean - ubwidth, ',', ubmean + ubwidth, ']'

        ## 2) estimate the lower bound w/ the same samples
        const = ubmean - quicksum(coefx[arc]*arcExp[arc] for arc in self.AllArcs).getValue()
        lb = const + self.estLowerBound(coefx)
        print 'The estimated lower bound = ', lb

        ## 3) estimate 95% C.I. with 15 batches on the lower bound
        lbs = np.zeros([self.M, 1])
        for i in range(self.M):
            evalvals = np.zeros([self.N, 1])
            coefx = {arc: 0.0 for arc in self.AllArcs}
            for k in range(self.N):
                sample = {}
                for c in self.Cset:
                    identity = np.random.rand()
                    samp = self.demstdev[c] * np.random.randn() + self.demmean[c]
                    samp = samp * (identity > self.probnodem[c])
                    sample[c] = samp
                samples.append(sample)
                Q, subQ = self.solveSSP(arcExp, sample)

                # compute f(\hat{x}, \ksi)
                evalvals[k] = quicksum(self.arcExpCost[arc] * arcExp[arc] for arc in self.AllArcs).getValue() + Q

                for arc in self.AllArcs:
                    coefx[arc] += (self.arcExpCost[arc] + subQ[arc]) / float(self.N)

            ubmean = np.mean(evalvals)
            const = ubmean - quicksum(coefx[arc] * arcExp[arc] for arc in self.AllArcs).getValue()

            lbs[i] = const + self.estLowerBound(coefx)

        lbmean = np.mean(lbs)
        lbwidth = np.std(lbs) / math.sqrt(self.M) * self.tvalue_14
        print '95% c.i. on the lower bound = [', lbmean - lbwidth, ',', lbmean + lbwidth, ']'


    def estLowerBound(self, coefx):
        ''' estimate the lower bound
        min \hat{l}(x;\hat{x})
        s.t. 0 <= x <= 200
        '''
        m = Model("Estimate LowerBound")
        m.params.logtoconsole = 0

        ''' Variables '''
        # arc expension amount
        x = {}
        for arc in self.AllArcs:
            x[arc] = m.addVar(obj=coefx[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' Constraints '''
        # Expansion upper bound (to make the 1st stage feasible region compact)
        for arc in self.AllArcs:
            m.addConstr(x[arc] <= self.expUB, name='ExpansionLimit_%s_%s' % (arc[0], arc[1]))

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL
        return m.objval


    def project(self, y):
        ''' project y to the space of X by solving the optimization problem as follows
        min \phi
        s.t. (x-y)^2 <= \phi
             0 <= x <= 200
        '''
        m = Model("Projection")
        m.params.logtoconsole=0

        ''' Variables '''
        # arc expension amount
        x = {}
        for arc in self.AllArcs:
            x[arc] = m.addVar(obj=0.0, name='Expansion_%s_%s' % (arc[0], arc[1]))

        # phi
        phi = m.addVar(obj=1.0, name='Phi')

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' Constraints '''
        # Expansion upper bound (to make the 1st stage feasible region compact)
        for arc in self.AllArcs:
            m.addConstr(x[arc] <= self.expUB, name='ExpansionLimit_%s_%s' % (arc[0], arc[1]))

        m.addQConstr(quicksum((x[arc] - y[arc])*(x[arc] - y[arc]) for arc in self.AllArcs) <= phi)

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL

        return {arc: x[arc].x for arc in self.AllArcs}


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

        # Demand constraints (sum_i {y_ijk} + z_jk >= d_jk) (Note: should be >=, not =, since the sampled d_jk might be negative)
        for c in self.Cset:
            m.addConstr(
                quicksum(transport[h, c] for h in self.Hset) + unmet[c] >= sampledDemScens[c], name='Demand_%s' % (c))
                # !! WARN: sampledDemScens[c] maybe negative, so use '>=' not '=='(probably results in infeasibility)

        # Line capacity constraints
        capcon = {}
        for arc in self.AllArcs:
            capcon[arc] = m.addConstr(transport[arc] <= self.curArcCap[arc] + arcExp[arc], name='LineCapacity_%s_%s' % (arc[0], arc[1]))

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
        return m.objval, {arc: capcon[arc].Pi for arc in self.AllArcs}

if __name__ == "__main__":
    mysolver = hw03q1()
    mysolver.read_data()
    mysolver.SA()
