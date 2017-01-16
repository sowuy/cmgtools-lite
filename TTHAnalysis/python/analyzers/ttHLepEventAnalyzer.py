from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.HeppyCore.utils.deltar import deltaR


        
class ttHLepEventAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(ttHLepEventAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

    def declareHandles(self):
        super(ttHLepEventAnalyzer, self).declareHandles()

    def beginLoop(self, setup):
        super(ttHLepEventAnalyzer,self).beginLoop(setup)
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('all events')
        count.register('accepted events')


    def makeHadTopDecays(self, event):
        event.lightJets = [ j for j in event.cleanJets if not j.btagWP("CSVv2IVFM") ]
        event.minMWjj   = 999
        event.minMWjjPt = 0
        event.bestMWjj   = 0
        event.bestMWjjPt = 0
        event.bestMTopHad   = 0
        event.bestMTopHadPt = 0
        for i1,j1 in enumerate(event.lightJets):
            for i2 in xrange(i1+1,len(event.lightJets)):
                j2 = event.lightJets[i2]
                jjp4 = j1.p4() + j2.p4()
                mjj  = jjp4.M()
                if mjj > 30 and mjj < event.minMWjj:
                    event.minMWjj = mjj
                    event.minMWjjPt = jjp4.Pt()
                if abs(mjj-80.4) < abs(event.bestMWjj-80.4):
                    event.bestMWjj = mjj
                    event.bestMWjjPt = jjp4.Pt()
                    for bj in event.bjetsLoose:
                        if deltaR(bj.eta(),bj.phi(),j1.eta(),j1.phi()) < 0.1 or deltaR(bj.eta(),bj.phi(),j2.eta(),j2.phi()) < 0.1: continue
                        tp4 = jjp4 + bj.p4()
                        mtop = tp4.M()
                        if abs(mtop-172) < abs(event.bestMTopHad - 172):
                            event.bestMTopHad = mtop
                            event.bestMTopHadPt = tp4.Pt()

    def process(self, event):
        self.readCollections( event.input )
        self.counters.counter('events').inc('all events')

        if hasattr(self.cfg_ana, 'minJets25'):
            n25 = len([ j for j in event.cleanJets if j.pt() > 25 ])
            if n25 < self.cfg_ana.minJets25:
                return False
        ###Added by Pablo for the Edge
        event.nPFLep5 = 0
        event.nPFHad10 = 0
        if hasattr(event, 'selectedIsoCleanTrack'):
            if len(event.selectedIsoTrack)>0:
                for myTrack in event.selectedIsoTrack:
                    if abs(myTrack.pdgId()) == 11 or abs(myTrack.pdgId()) == 13:
                        if myTrack.pt()>5 and myTrack.absIso/myTrack.pt()<0.2:
                            event.nPFLep5 += 1
                    if abs(myTrack.pdgId()) == 211:
                        if myTrack.pt()>10 and myTrack.absIso/myTrack.pt()<0.1:
                            event.nPFHad10 += 1
        #################################
        self.makeHadTopDecays(event)

        self.counters.counter('events').inc('accepted events')
        return True
