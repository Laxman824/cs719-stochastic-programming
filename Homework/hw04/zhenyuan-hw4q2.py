from gurobipy import *
import time
import numpy
import math

class hw04q2:
    def __init__(self):
        ### Base data
        self.T = 5  ## (Number of stages)
        self.fracoutcomes = [0.05, 0.1, 0.15]  ## set of possible fraction destroyed outcomes in each period
        self.classes = [1,2,3,4]
        self.initArea = {1:8000, 2:10000, 3:20000, 4:60000}
        self.yields = {1:0, 2:250, 3:510, 4:710}

        ### The following code builds data structures for representing the scenario tree
        ### It's not necessarily the most efficient implementation -- just meant to be relatively simple

        ### This will be a list of indices to the nodes in the tree, will serve as keys to other data structures
        self.Nodes = []

        ### This dictionary maps Nodes to Stages. NodeStage[n] gives stage of node n
        self.NodeStage = {}

        ## This dictionary maps each stages to a list of nodes which are in that stage
        self.StageNodes = {}

        ### This dictionary maps Nodes to a list of Chilren nodes. NodeChildren[n] gives list of nodes, which are the children
        ### of node n, NodeChildren[n] is an empty list if NodeStage[n] == T
        self.NodeChildren = {}

        self.NodeParent = {}

        ### This dictionary maps Nodes to a list of numbers, representing the history of fractions lost up to that node in tree
        self.NodeFracHist = {}

        ### Initialize with root node
        self.curnodeid = 1
        self.Nodes.append(self.curnodeid)
        self.NodeStage[self.curnodeid] = 1
        self.StageNodes[1] = [self.curnodeid]
        self.NodeFracHist[self.curnodeid] = [0.0]  ### I am using convention that nothing is "lost" at stage t=1 since we already know how
        ###much is available
        self.NodeChildren[1] = []

        ### Build the tree moving forward in stages
        for t in range(2, self.T + 1):  ## 2...T
            self.StageNodes[t] = []
            for n in self.StageNodes[t - 1]:
                for f in self.fracoutcomes:
                    self.curnodeid += 1
                    self.NodeFracHist[self.curnodeid] = self.NodeFracHist[n] + [f]
                    self.NodeChildren[n] += [self.curnodeid]
                    self.NodeChildren[self.curnodeid] = []
                    self.NodeStage[self.curnodeid] = t
                    self.StageNodes[t] += [self.curnodeid]
                    self.Nodes += [self.curnodeid]

                    self.NodeParent[self.curnodeid] = n

        ### The following recursively defined function gets all leaf nodes (which correspond to scenarios) that originate from a node n
        def getLeafs(n):
            if self.NodeStage[n] == self.T:
                return [n]
            else:
                tmp = []
                for c in self.NodeChildren[n]:
                    tmp += getLeafs(c)
                return tmp

        ### The full set of scenarios corresponds to the set of node indices that are Leafs in the tree from node 1
        self.Scens = getLeafs(1)

        ### Let's see what we have!

        ### here is the full observation in each scenario (path through tree)
        for n in self.Scens:
            print(self.NodeFracHist[n])

        self.Nscen = len(self.Scens)
        self.N = len(self.Nodes)
        print("number of scenarios: %d" % self.Nscen)

        ### print the node indices in each stage
        for t in range(1, self.T + 1):
            print(self.StageNodes[t])

        ## as an example, let's pick out node 3 and look at key information:
        print("Node 3 time stage: %d" % self.NodeStage[3])
        print("Node 3 scenarios that share this same history:")
        print(getLeafs(3))

    def solveNodeForm(self):
        '''
        solve the multi-stage stochastic programming w/ node-based formulation
        :return:
        '''

        m = Model("Node Form")
        m.params.logtoconsole = 0
        ''' variables '''
        x = {}; y = {}; z = {}; v = {}
        for i in self.Nodes:
            v[i] = m.addVar(obj=1.0, name='profit_%s' % (i))
            y[i] = m.addVar(name='yield_%s' % (i))
            for c in self.classes:
                x[c,i] = m.addVar(name='harvest_acre_%s_%s' % (c, i))
                z[c,i] = m.addVar(name='left_acre_%s_%s' % (c, i))

        m.modelSense = GRB.MAXIMIZE
        m.update()

        ''' constraints '''
        for i in self.Nodes:
            if i == self.Nodes[0]:
                ## z_1
                for c in self.classes:
                    m.addConstr(z[c, i] == self.initArea[c] - x[c, i], name='z_{}_{}'.format(c,i))
            else:
                ## z_i
                ancestor = self.NodeParent[i]
                xi = self.NodeFracHist[i][-1]
                m.addConstr(z[1,i] == xi*(quicksum(z[c,ancestor] for c in self.classes))
                            + x[2,i] + x[3,i] + x[4,i], name='z_1_{}'.format(i))
                m.addConstr(z[2,i] == (1-xi) * z[1,ancestor] - x[2,i], name='z_2_{}'.format(i))
                m.addConstr(z[2,i] == (1-xi) * z[2,ancestor] - x[3,i], name='z_3_{}'.format(i))
                m.addConstr(z[2,i] == (1-xi) * (z[3,ancestor]+z[4,ancestor]) - x[4,i], name='z_4_{}'.format(i))

            m.addConstr(y[i] == quicksum(self.yields[c] * x[c,i] for c in self.classes), name='y_{}'.format(i))
            m.addConstr(v[i] <= 10*y[i], name='v1_{}'.format(i))
            m.addConstr(v[i] <= 10*98e6 + 7*(y[i]-98e6), name='v2_{}'.format(i))
            m.addConstr(v[i] <= 10*98e6 + 49e6 + 5*(y[i]-105e6), name='v3_{}'.format(i))

        m.update()
        m.optimize()
        assert m.status == GRB.Status.OPTIMAL
        print("Expected Income = ", m.objVal)

    def solveScenForm(self):
        '''
        solve the multi-stage stochastic programming w/ scenario-based formulation
        :return:
        '''
        pass

if __name__ == "__main__":
    mysolver = hw04q2()
    mysolver.solveNodeForm()