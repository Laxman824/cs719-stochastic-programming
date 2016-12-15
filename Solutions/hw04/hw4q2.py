from gurobipy import *
import time
import numpy
import math

### Base data

T=5  ## (Number of stages)
fracoutcomes = [0.05,0.1,0.15] ## set of possible fraction destroyed outcomes in each period

### The following code builds data structures for representing the scenario tree 
### It's not necessarily the most efficient implementation -- just meant to be relatively simple

### This will be a list of indices to the nodes in the tree, will serve as keys to other data structures
Nodes = []

### This dictionary maps Nodes to Stages. NodeStage[n] gives stage of node n
NodeStage = {}

## This dictionary maps each stages to a list of nodes which are in that stage
StageNodes = {}

### This dictionary maps Nodes to a list of Chilren nodes. NodeChildren[n] gives list of nodes, which are the children
### of node n, NodeChildren[n] is an empty list if NodeStage[n] == T
NodeChildren = {}

### This dictionary maps Nodes to a list of numbers, representing the history of fractions lost up to that node in tree 
NodeFracHist = {}

### Initialize with root node
curnodeid = 1
Nodes.append(curnodeid)
NodeStage[curnodeid]=1
StageNodes[1] = [curnodeid]
NodeFracHist[curnodeid]=[0.0] ### I am using convention that nothing is "lost" at stage t=1 since we already know how
                              ###much is available
NodeChildren[1] = []

### Build the tree moving forward in stages
for t in range(2,T+1): ## 2...T
	StageNodes[t] = []
	for n in StageNodes[t-1]:
		for f in fracoutcomes:
			curnodeid+=1
			NodeFracHist[curnodeid]=NodeFracHist[n]+[f]
			NodeChildren[n] += [curnodeid]
			NodeChildren[curnodeid]=[]
			NodeStage[curnodeid]=t
			StageNodes[t]+=[curnodeid]
			Nodes+=[curnodeid]

### The following recursively defined function gets all leaf nodes (which correspond to scenarios) that originate from a node n
def getLeafs(n):
	if NodeStage[n] == T:
		return [n]
	else:
		tmp = []
		for c in NodeChildren[n]:
			tmp += getLeafs(c)	
		return tmp

### The full set of scenarios corresponds to the set of node indices that are Leafs in the tree from node 1
Scens = getLeafs(1)


### Let's see what we have! 

### here is the full observation in each scenario (path through tree)
for n in Scens:
	print(NodeFracHist[n])

Nscen = len(Scens)
print("number of scenarios: %d" % Nscen)

### print the node indices in each stage
for t in range(1,T+1):
	print(StageNodes[t])

## as an example, let's pick out node 3 and look at key information:
print("Node 3 time stage: %d" % NodeStage[3])
print("Node 3 scenarios that share this same history:")
print(getLeafs(3))



#### Here is some other data defined as parameters
prices = { 1 : 100,
			  2: 70,
			  3: 50 }
plevs = prices.keys() 
limits = { 1 : 98,
			  2 : 7,
			  3 : 100000}
initacres = { 1 : 8,
				  2 : 10,
				  3 : 20,
				  4 : 60 }
classes = initacres.keys()
yields = { 2 : 0.25,
			  3 : 0.51,
			  4 : 0.71 }


#### Now begin the implementation
m = Model("ext")

acres = {}
for t in range(1,T+1):
	for s in Scens:
		for i in [0,1,2,3,4]: ### includes 0 for keeping track of what is harvested now 
			acres[i,t,s] = m.addVar(name='acre%d_%d_%d' % (t,s,i))

harvest = {}
for t in range(1,T+1):
	for s in Scens:
		for i in [2,3,4]:  ## only the harvestable classes
			harvest[i,t,s] = m.addVar(name='harvest%d_%d_%d'% (t,s,i))

sellAtPrice = {}
for t in range(1,T+1):
	for s in Scens:
		for k in plevs:
			sellAtPrice[k,t,s] = m.addVar(name='sellatprice%d_%d_%d' % (k,t,s), obj=prices[k]/float(len(Scens)))

cumSoldAtPrice = {}
for t in range(1,T+1):
	for s in Scens:
		for k in [1,2]:
			cumSoldAtPrice[k,t,s] = m.addVar(name='cumsoldatprice%d_%d_%d' % (k,t,s),ub=limits[k])
		cumSoldAtPrice[3,t,s] = m.addVar(name='cumsoldatprice%d_%d_%d' % (3,t,s))


m.modelSense = GRB.MAXIMIZE
m.update()

### now the constraints per scenario

for s in Scens:
	for t in range(2,T+1):
		for i in [2,3]:
			m.addConstr(acres[i,t,s] == acres[i-1,t-1,s]*(1.0-NodeFracHist[s][t-1]) - harvest[i,t,s]) ##NodeFracHist has index t-1 due to zero
																										## indexing python convention
		m.addConstr(acres[4,t,s] == (acres[3,t-1,s]+acres[4,t-1,s])*(1.0-NodeFracHist[s][t-1]) - harvest[4,t,s])
		m.addConstr(acres[1,t,s] == acres[0,t-1,s]+sum(acres[i,t-1,s]*NodeFracHist[s][t-1] for i in [1,2,3,4]))
		m.addConstr(acres[0,t,s] == sum(harvest[i,t,s] for i in [2,3,4]))
	### do t=1 separately
	for i in [2,3,4]:
		m.addConstr(acres[i,1,s]==initacres[i] - harvest[i,1,s])
	m.addConstr(acres[1,1,s]==initacres[1])
	m.addConstr(acres[0,1,s]==sum(harvest[i,1,s] for i in [2,3,4]))

	for t in range(1,T+1):
		m.addConstr(sum(yields[i]*harvest[i,t,s] for i in [2,3,4]) == sum(sellAtPrice[k,t,s] for k in plevs))
		for k in plevs:
			if t==1:
				m.addConstr(cumSoldAtPrice[k,t,s] == sellAtPrice[k,t,s]) 
			else:
				m.addConstr(cumSoldAtPrice[k,t,s] == cumSoldAtPrice[k,t-1,s]+sellAtPrice[k,t,s])


### Finally add the nonancipativity constraints			
for n in Nodes:  
	scenSet = getLeafs(n)
	N = float(len(scenSet))
	t=NodeStage[n]
	for i in [0,1,2,3,4]:
		curRhs = sum(acres[i,t,s] for s in scenSet)
		for s in scenSet:	
			m.addConstr(N*acres[i,t,s] == curRhs)
	for k in plevs:
		curRhs = sum(cumSoldAtPrice[k,t,s] for s in scenSet)
		for s in scenSet:
			m.addConstr(N*cumSoldAtPrice[k,t,s] == curRhs)

   ### these NA constraints are not needed
	for k in plevs:
		curRhs = sum(sellAtPrice[k,t,s] for s in scenSet)
		for s in scenSet:
			m.addConstr(N*sellAtPrice[k,t,s] == curRhs)
	for i in [2,3,4]:
		curRhs = sum(harvest[i,t,s] for s in scenSet)
		for s in scenSet:
			m.addConstr(N*harvest[i,t,s] == curRhs)
			

m.update()			
m.optimize()

for s in Scens:
	print("scenario %d" % s)
	for t in range(1,T+1):
		print("stage %d" % t)
		print("harvest: [%f, %f, %f]" % (harvest[2,t,s].x,harvest[3,t,s].x,harvest[4,t,s].x))
		print("sellatprice: [%f, %f, %f]" % (sellAtPrice[1,t,s].x,sellAtPrice[2,t,s].x,sellAtPrice[3,t,s].x))



