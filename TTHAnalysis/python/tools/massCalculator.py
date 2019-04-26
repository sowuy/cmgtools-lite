from CMGTools.TTHAnalysis.treeReAnalyzer2 import *
class MassCalculator:
    def __init__(self, lepFOColl, lepLooseColl):
        self.lepFOColl    = lepFOColl
        self.lepLooseColl = lepLooseColl
        
    def listBranches(self):
        return [ 'mZ1','m4l','minMllAFAS','nSFOSpairs' ] 

    def __call__(self,event):
        ret = {} 
        ret['mZ1']        = -1
        ret['m4l']        = -1
        ret['minMllAFAS'] = -1
        ret['nSFOSpairs'] = 0

        foLeps = [ l for l in Collection(event, self.lepFOColl)]
        loLeps = [ l for l in Collection(event, self.lepLooseColl)] 

        nossfpairs = [] 
        for l1 in foLeps+loLeps:
            for l2 in foLeps+loLeps:
                if l1==l2: continue
                if ret['minMllAFAS'] < 0 or ret['minMllAFAS'] > (l1.p4()+l2.p4()).M():
                    ret['minMllAFAS'] = (l1.p4()+l2.p4()).M()
                if l1.pdgId*l2.pdgId == -121 or l1.pdgId*l2.pdgId == -169:
                    if (ret['mZ1'] < 0 or abs( ret['mZ1']-91) > abs( (l1.p4()+l2.p4()).M() -91)):
                        ret['mZ1'] = (l1.p4()+l2.p4()).M()
                        if l1 in nossfpairs: continue
                        if l2 in nossfpairs: continue
                        nossfpairs.extend( [l1, l2]) 
                        ret['nSFOSpairs']  = ret['nSFOSpairs'] + 1 

                        if len(nossfpairs) == 4: 
                            ret['m4l'] = (nossfpairs[0].p4()+nossfpairs[1].p4()+nossfpairs[2].p4()+nossfpairs[3].p4()).M()
                        
        return ret

