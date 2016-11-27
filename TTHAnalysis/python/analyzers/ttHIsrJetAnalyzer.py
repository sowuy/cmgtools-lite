import math
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.HeppyCore.utils.deltar import deltaR

class ttHIsrJetAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(ttHIsrJetAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName) 
        self.jetPt = cfg_ana.jetPt

    def declareHandles(self):
        super(ttHIsrJetAnalyzer, self).declareHandles()
    
    def beginLoop(self, setup):
        super(ttHIsrJetAnalyzer,self).beginLoop(setup)
        
    def nIsrMatch(self, event):
        
        event.nisrMatch = 0
        jets = [j for j in event.cleanJets if j.pt()>self.jetPt]
        for jet in jets:
            matched = False
            for mc in event.genParticles: #WithMotherId:
                if matched: break
                if (mc.status()!=23 or abs(mc.pdgId())>5): continue
                momid = abs(mc.mother().pdgId())
                if not (momid==6 or momid==23 or momid==24 or momid==25 or momid>1e6): continue  # MT
                #check against daughter in case of hard initial splitting
                for idau in range(mc.numberOfDaughters()):
                    dR = deltaR(jet.eta(), jet.phi(), mc.daughter(idau).p4().eta(), mc.daughter(idau).p4().phi())
                    if dR<0.3:
                        matched = True
                        break
            if not matched:
                event.nisrMatch+=1
                
            

    def process(self, event):
        self.readCollections( event.input )

        event.nisrMatch = 0
    
        if self.cfg_comp.isMC:
            self.nIsrMatch(event)

        return True
