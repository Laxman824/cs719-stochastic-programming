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

