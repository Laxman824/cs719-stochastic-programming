from gurobipy import *
import time
import numpy
import math

D = [0,1,2,3]
cost = {}
cost[0] = 2400
cost[1] = 2600
cost[2] = 3000
cost[3] = 2000
time = {}
time[0] = 0.5
time[1] = 0.75 
time[2] = 1.0 
time[3] = 0.3
diffmult = 1.2
eps = 0.05
SampSize = 100
eligDepts = {}
eligDepts[0] = [0,3]
eligDepts[1] = [1,3]
eligDepts[2] = [0,1,2,3]
eligDepts[3] = [3]
revElig = {}
revElig[0]=[0,2]
revElig[1]=[1,2]
revElig[2]=[2]
revElig[3]=[0,1,2,3]


def genSample(N):
	basearr = numpy.random.poisson(lam=[50,35,30,10],size=(N,4))
	postarr = basearr.astype(float)
	postarr[:,1] = basearr[:,1]+0.2*basearr[:,0]
	postarr[:,2] = basearr[:,2]+0.2*basearr[:,0]+0.1*basearr[:,1]
	postarr[:,3] = basearr[:,3]+0.3*basearr[:,0]+0.6*basearr[:,1]+0.2*basearr[:,2]
	print postarr
	return postarr.astype(int)


def buildSolveExtensive(sample):
	
	m = Model("cc-ext")
	N = len(sample)

	staffing = {}
	for i in D:
		staffing[i] = m.addVar(name='Staff%d' % i,vtype=GRB.INTEGER,obj=cost[i])

	z = {}
	for k in range(N):
		z[k] = m.addVar(name='z%d' %k,vtype=GRB.BINARY)

	assign = {}
	for i in D:
		for j in eligDepts[i]:
			for k in range(N):
				assign[i,j,k] = m.addVar(name='assign%d_%d_%d' %(i,j,k))	

	m.modelSense = GRB.MINIMIZE
	m.update()

	m.addConstr(quicksum(z[k] for k in range(N)) >= (1.0-eps)*float(N))

	for k in range(N):
		for j in D:
			m.addConstr(quicksum(assign[i,j,k] for i in revElig[j]) >= sample[k][j]*z[k])
		for i in D:
			m.addConstr(time[i]*assign[i,i,k] + quicksum(time[j]*diffmult*assign[i,j,k] for j in eligDepts[i] if j != i) <= 12*staffing[i])

	m.update()
	m.optimize()

	xsol = []
	for j in D:
		print("staffing in department %d = %d" % (j,staffing[j].x))
		xsol.append(staffing[j].x)

	return [m.objval, xsol]

minobjval = 0.0
for k in range(10):	
	sample = genSample(SampSize)
	[objval, xsol] = buildSolveExtensive(sample)
	if k==0 or objval < minobjval:
		minobjval = objval

print("min over 10 solves (99 percent confidence lower bound): %f" % minobjval)

## arbitrarily chooose last sampled solution to estimate its objective value
print("objective value of last generated solution: %f" % objval)
largesamp = genSample(100000)
ninfeas = 0

for k in range(100000):
	f = Model("feascheck")
	f.params.logtoconsole=0

	assign = {}
	for i in D:
		for j in eligDepts[i]:
			assign[i,j] = f.addVar(name='assign%d_%d' %(i,j))	

	slack = {}
	for j in D:
		slack[j] = f.addVar(name='slack%d' %i, obj=1.0)

	f.modelSense = GRB.MINIMIZE
	f.update()

	for j in D:
		f.addConstr(quicksum(assign[i,j] for i in revElig[j]) + slack[j] >= largesamp[k][j])
	for i in D:
		f.addConstr(time[i]*assign[i,i] + quicksum(time[j]*diffmult*assign[i,j] for j in eligDepts[i] if j != i) <= 12*xsol[i])

	f.update()
	f.optimize()

	if f.objval > 0.00001:
		ninfeas = ninfeas + 1

phat = float(ninfeas)/100000.0
stderr = math.sqrt(phat*(1-phat)/100000.0)

print("estimated probability of violation: %f" % phat)
print("confidence interval: [%f, %f]" % (phat - 1.96*stderr, phat+1.96*stderr))

