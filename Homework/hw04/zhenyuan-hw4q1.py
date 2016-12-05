from gurobipy import *
import numpy as np
import math

class hw04q1:
    def __init__(self):
        self.c = {1:2400, 2:3000, 3:2600, 4:2000}
        self.a = {1:0.5, 2:1.0, 3:0.75, 4:0.3}
        self.mu = {1:50, 2:35, 3:30, 4:10}
        self.N = 100
        self.epsilon = 0.05
        self.T = 10
        self.M = 100000
        self.D = [1,2,3,4]
        self.arc = tuplelist([(1,4), (2,4), (3,4), (2,1), (2,3)])
        self.K = list(range(self.N))
        self.tvalue_9 = 2.262
        self.tvalue_inf = 1.960

    def bigMStrengthen(self, sample):
        '''
        Big M Strengthening
        :param sample:
        :return:
        '''
        ''' variables '''
        m = Model("LP relaxation")
        m.params.logtoconsole = 0

        x = {}
        z = {}
        for i in self.D:
            x[i] = m.addVar(name='nurse_%s' % (i))
            for s in range(self.N):
                z[i, s] = m.addVar(name='eff_nurse_%s_%s' % (i, s))

        w = {}
        for s in self.K:
            w[s] = m.addVar(ub=1.0, name='isfeasible_%s' % (s))

        y = {}
        for arc in self.arc:
            for s in self.K:
                y[arc, s] = m.addVar(name='move_%s_%s_%s' % (arc[0], arc[1], s))

        m.modelSense = GRB.MAXIMIZE
        m.update()

        ''' constraints '''
        for i in self.D:
            for s in self.K:
                m.addConstr(
                    x[i] - quicksum(y[arc, s] for arc in self.arc.select(i, '*')) + 1 / 1.2 * quicksum(
                        y[arc, s] for arc in self.arc.select('*', i)) == z[i, s], name='nurse conservation')
                m.addConstr(12 * z[i, s] >= w[s] * self.a[i] * sample[i][s], name="demand constraint")
        m.addConstr(quicksum(w.values()) >= self.N * (1 - self.epsilon), name='probability constraint')

        m.update()

        M={}
        for r in self.D:  # m
            for k in self.K:  # N
                # modify objective function
                z[r, k].obj = -12.0
                m.update()
                m.optimize()
                assert m.status == GRB.Status.OPTIMAL
                M[r, k] = m.objVal + self.a[r]*sample[r][k]

                # recover objective function
                z[r, k].obj = 0.0
                m.update()
        return M

    def solveMIP(self, M, sample):
        '''
        Solve MIP
        :param M: coefficients
        :param sample:
        :return:
        '''
        ''' variables '''
        m = Model("MIP model")
        m.params.logtoconsole = 0
        x = {}
        z = {}
        for i in self.D:
            x[i] = m.addVar(obj=self.c[i], vtype=GRB.INTEGER, name='nurse_%s' % (i))
            for k in range(self.N):
                z[i, k] = m.addVar(obj=0.0, name='eff_nurse_%s_%s' % (i, k))
        w = {}
        for k in self.K:
            w[k] = m.addVar(obj=0.0, vtype=GRB.BINARY, name='isfeasible_%s' % (k))

        y = {}
        for arc in self.arc:
            for k in self.K:
                y[arc, k] = m.addVar(obj=0.0, name='move_%s_%s_%s' % (arc[0], arc[1], k))

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' constraints '''
        for i in self.D:
            for k in self.K:
                m.addConstr(
                    x[i] - quicksum(y[arc, k] for arc in self.arc.select(i, '*')) + 1 / 1.2 * quicksum(
                        y[arc, k] for arc in self.arc.select('*', i)) == z[i, k],
                    name='nurse conservation')
                m.addConstr(12 * z[i, k] + M[i, k] * (1 - w[k]) >= self.a[i] * sample[i][k], name="demand constraint")
        m.addConstr(quicksum(w.values()) >= self.N * (1 - self.epsilon), name='probability constraint')

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL
        x_res = {}
        y_res = {}
        z_res = {}
        w_res = {}
        for i in self.D:
            x_res[i] = x[i].x
            for k in range(self.N):
                z_res[i, k] = z[i, k].x
        for k in self.K:
            w_res[k] = w[k].x
        for arc in self.arc:
            for k in self.K:
                y_res[arc, k] = y[arc, k]

        return m.objVal, (x_res, y_res, z_res, w_res)

    def evalSol(self, soln, sample):
        ''' Evaluate the feasibility probability of the candidate solution '''
        x = soln[0]
        nFeasible = 0.0

        ''' variables '''
        m = Model("MIP model")
        m.params.logtoconsole = 0
        z = {}
        for i in self.D:
            z[i] = m.addVar(obj=0.0, name='eff_nurse_%s' % (i))

        y = {}
        for arc in self.arc:
            y[arc] = m.addVar(obj=0.0, name='move_%s_%s' % (arc[0], arc[1]))

        m.modelSense = GRB.MINIMIZE
        m.update()

        ''' constraints '''
        for i in self.D:
            m.addConstr(
                x[i] - quicksum(y[arc] for arc in self.arc.select(i, '*')) + 1 / 1.2 * quicksum(
                    y[arc] for arc in self.arc.select('*', i)) == z[i], name='nurse conservation')
        for k in self.K:
            demConstr = {}
            for i in self.D:
                demConstr[i] = m.addConstr(12 * z[i] >= self.a[i] * sample[i][k], name="demand constraint")
            m.update()
            m.optimize()
            for i in self.D:
                m.remove(demConstr[i])
            m.update()
            if m.status == GRB.Status.OPTIMAL:
                nFeasible += 1

        return nFeasible/self.N

    def solve(self):
        ## sampling
        samples = {i: np.random.poisson(self.mu[i], (self.T, self.N)) for i in self.D}

        objvals = np.zeros([self.T, 1])
        soln = {}
        ## run the algorithm T times
        minCost = float("inf")
        for t in range(self.T):
            sample = {i: samples[i][t] for i in samples}
            sample[2] += 0.2*sample[1]
            sample[3] += 0.2*sample[1] + 0.1*sample[2]
            sample[4] += 0.3*sample[1] + 0.6*sample[2] + 0.2*sample[3]
            ## compute big M strengthening
            M = self.bigMStrengthen(sample)
            ## solve MIP
            objvals[t], soln[t] = self.solveMIP(M, sample)
            if objvals[t] < minCost:
                minCost = objvals[t]
                candSoln = soln[t]
        ## statistically estimate the lower bound of the optimal value
        lbmean = np.mean(objvals)
        lbwidth = np.std(objvals) / math.sqrt(self.T) * self.tvalue_9
        print 'extimated lower bound = ', lbmean
        print '95% ci on lower bound = [', lbmean - lbwidth, ',', lbmean + lbwidth, ']'

        ## statistically estimate the chance constraint of the candidate solution
        samples = {i: np.random.poisson(self.mu[i], (self.M, self.N)) for i in self.D}
        prob = np.zeros([self.M, 1])
        solId = np.random.choice(self.N)
        for m in range(self.M):
            sample = {i: samples[i][m] for i in samples}
            sample[2] += 0.2*sample[1]
            sample[3] += 0.2*sample[1] + 0.1*sample[2]
            sample[4] += 0.3*sample[1] + 0.6*sample[2] + 0.2*sample[3]

            prob[m] = self.evalSol(candSoln, sample)
        meanProb = np.mean(prob)
        widthProb = np.std(prob) / math.sqrt(self.M) * self.tvalue_inf
        print 'extimated feasible probability = ', meanProb
        print '95% ci on feasible probability = [', meanProb - widthProb, ',', meanProb + widthProb, ']'

if __name__ == "__main__":
    mysolver = hw04q1()
    mysolver.solve()