import cPickle
from gurobipy import *
import time

class hw03q2:
    def __init__(self):
        self.avarweight = 0.5
        self.alpha = 0.95

    def read_data(self, fileName="nd848.pdat"):
        ### Read data from file you choose: commont/uncomment to choose the different files
        ### This file was generated from a separate python file, using the `cPickle' module
        ### This is just for convenience -- data can be read in many ways in Python

        dfile = open('nd151020500.pdat','r')  # 'nd15-10-20-500.pdat'

        self.Fset = cPickle.load(dfile)  # set of facilities (list of strings)
        self.Hset = cPickle.load(dfile)  # set of warehouses (list of strings)
        self.Cset = cPickle.load(dfile)  # set of customers (list of strings)
        self.Sset = cPickle.load(dfile)  # set of scenarios (list of strings)
        self.arcExpCost = cPickle.load(dfile)  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
        self.facCap = cPickle.load(dfile)   # facility capacities (dictionary mapping F to floats)
        self.curArcCap = cPickle.load(dfile)  # current arc capacities (dictionary mapping (i,j) to floats, where either
                                         # i is facility, j is warehouse, or i is warehouse and j is customer
        self.unmetCost = cPickle.load(dfile)  # penalty for unment customer demands (dicationary mapping C to floats)
        self.demScens = cPickle.load(dfile)  # demand scenarios (dictionary mapping (i,k) tuples to floats, where i is customer, k is scenario
        dfile.close()

        ### Define sets of arcs (used as keys to dictionaries)
        self.FHArcs = [(i,j) for i in self.Fset for j in self.Hset]  ## arcs from facilities to warehouses
        self.HCArcs = [(i,j) for i in self.Hset for j in self.Cset]   ## arcs from warehouses to customers
        self.AllArcs = self.FHArcs + self.HCArcs

        ### Make them Gurobi tuplelists
        self.FHArcs = tuplelist(self.FHArcs)
        self.HCArcs = tuplelist(self.HCArcs)
        self.AllArcs = tuplelist(self.AllArcs)

        self.nscen = len(self.Sset)

    def RM(self):
        m = Model("Mean Risk Model")

        ''' Variables '''
        # First stage vars, arc expansion amounts
        arcExpAmount = {}
        for arc in self.AllArcs:
            arcExpAmount[arc] = m.addVar(obj=self.arcExpCost[arc], name='Expansion_%s_%s' % (arc[0], arc[1]))

        # gamma
        gamma = m.addVar(obj=self.avarweight, lb=-GRB.INFINITY, name='gamma')

        # w_k
        w = {}
        for k in self.Sset:
            w[k] = m.addVar(obj=self.avarweight/float((1-self.alpha) * self.nscen), name='w_%s' %k)

        # Second stage vars, unmet demands (y_ck)
        unmet = {}
        for c in self.Cset:
            for k in self.Sset:
                unmet[c,k] = m.addVar(obj=float(self.unmetCost[c])/self.nscen, name='Unmet_%s_%s' % (c,k))

        # Intermediate vars, transport (z_ijk)
        transport = {}
        for arc in self.AllArcs:
            for k in self.Sset:
                transport[arc[0], arc[1], k] = m.addVar(obj=0, name="Transport_%s_%s_%s" % (arc[0], arc[1], k))

        m.modelSense = GRB.MINIMIZE
        m.update()


        ''' Constraints '''
        # Production constraints (sum_h {z_fhk} - b_f <= 0)
        for f in self.Fset:
            for k in self.Sset:
                m.addConstr(
                    quicksum(transport[f,h,k] for h in self.Hset) <= self.facCap[f], name='Capacity_%s_%s' % (f,k))

        # Demand constraints (sum_h {z_hck} + y_ck = d_ck)
        for c in self.Cset:
            for k in self.Sset:
                m.addConstr(
                    quicksum(transport[h,c,k] for h in self.Hset) + unmet[c,k] == self.demScens[(c,k)],name='Demand_%s_%s' % (c,k))


        # Line capacity constraints (z_ijk <= x_ij + u_ij)
        for arc in self.AllArcs:
            for k in self.Sset:
                m.addConstr(transport[arc[0], arc[1], k] <= self.curArcCap[arc] + arcExpAmount[arc],
                            name='LineCapacity_%s_%s_%s' % (arc[0], arc[1], k))

        # Transportation contraints (sum_f {z_fhk} >= sum_c {z_hck})
        for h in self.Hset:
            for k in self.Sset:
                m.addConstr(
                    # quicksum(transport[f,h,k] for f in Fset) >= quicksum(transport[h,c,k] for c in Cset),
                    # name="Transportation_%s_%s" % (h, k))
                    quicksum(transport[f, h, k] for f, h in self.FHArcs.select('*', h)) >= quicksum(
                        transport[h, c, k] for h, c in self.HCArcs.select(h, '*')),
                    name="Transportation_%s_%s" % (h, k))

        # risk constraints
        for k in self.Sset:
            m.addConstr(w[k] >= quicksum(self.arcExpCost[arc] * arcExpAmount[arc] for arc in self.AllArcs)
                        + quicksum(self.unmetCost[c] * unmet[c,k] for c in self.Cset) - gamma, name="Risk_%s" %k)

        m.update()

        print('\n-------------------------- Extensive Form  ------------------------------')
        m.optimize()

        assert m.status == GRB.Status.OPTIMAL

        print('Total Cost: {}'.format(m.objVal))

        E_Z = quicksum(self.arcExpCost[arc] * arcExpAmount[arc].x for arc in self.AllArcs).getValue()\
              + 1.0/self.nscen * quicksum(self.unmetCost[c] * unmet[c,k] for c in self.Cset for k in self.Sset).getValue()
        AVaR_Z = gamma.x + 1.0/((1-self.alpha) * self.nscen)*quicksum(w[k] for k in self.Sset).getValue()
        print('Expected Cost: {}'.format(E_Z))
        print('Average Value at Risk of Cost: {}'.format(AVaR_Z))


if __name__ == "__main__":
    mysolver = hw03q2()
    mysolver.read_data()
    mysolver.RM()
