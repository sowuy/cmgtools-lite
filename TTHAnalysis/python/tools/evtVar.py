from CMGTools.TTHAnalysis.treeReAnalyzer import *

class EvtVars: 
    def __init__(self,vars):
        self.vars = vars
    
    def listBranches(self):
        return [(k,'F') for k in self.vars]
    
    def __call__(self,event):
        ret = {} 
        for k,v in self.vars.iteritems():
            ret[k] = v(event)
        return ret
