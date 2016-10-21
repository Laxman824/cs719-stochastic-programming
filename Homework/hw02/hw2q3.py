import cPickle
from gurobipy import *

### This file was generated from a separate python file, using the `cPickle' module
### This is just for convenience -- data can be read in many ways in Python

dfile = open('forest4-10-31.pdat','r')

B = cPickle.load(dfile)  # set of base locations (list of strings) 
F = cPickle.load(dfile)  # set of fire districts (list of strings)
S = cPickle.load(dfile)  # set of scenarios (list of strings0
c = cPickle.load(dfile)  # dictionary of cost parameters c[i,j] is unit cost to relocate from i in B to j in B 
h = cPickle.load(dfile)  # dictionary or purchase params, h[i] is unit cost to add new resource at i in B
init = cPickle.load(dfile)  # dictionary of initial resource placements, init[i] is amt of resource at i in B
closesets = cPickle.load(dfile)  # dictionary of "close enough sets". closesets[f] is a list of base locations (subset of B)
demscens = cPickle.load(dfile)  # dictionary of demand scnearios. demscens[s,f] is demand of resource in scenario s at district f 
costscens = cPickle.load(dfile) # dictionary of cost scenarios. costscnes[s,f] is cost of resource shortage in scenario s at district f
budget = 500.0   ### hard coded here!!!
dfile.close() 

### NoTE: demscens and costscens are very sparse. In particular, for a given s, there are only one or two districts f
### that have nonzero demscens vals (and only these have costscens vals, since without demand there is no need for a cost param)
### Here I define the set of keys which exist in these sets

SFkeys = demscens.keys()
SFkeys = tuplelist(SFkeys)  ## make tuplelist for easy selection

### it may also be useful to have a "reverse" version of the dictionary closesets, which provides for each facility i in
### B, the set of districts it can serve. This is constructed here
closedists = {}
for i in B:
    closedists[i] = []

for f in F:
    for i in closesets[f]:
        closedists[i].append(f)

### This is just a check of the data. Probably you want to comment/delete these lines once you see the structure 

print B 
print F 
print S 
print c 
print h
print init
print closesets
print demscens
print costscens
print budget
