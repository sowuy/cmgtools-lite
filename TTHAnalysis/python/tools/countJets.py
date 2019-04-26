from CMGTools.TTHAnalysis.treeReAnalyzer import *



class CountJets:
    def __init__(self, jetPts, jetSel, btagL_thr, btagM_thr, btag, doJetSums):
        self.jetPts    = jetPts    
        self.jetSel    = jetSel    
        self.btagL_thr = btagL_thr 
        self.btagM_thr = btagM_thr 
        self.btag      = btag
        self.systsJEC = {0:"_nom"  , 1:"_jesTotalUp", -1:"_jesTotalDown"}
        self.doJetSums  = doJetSums


    def listBranches(self):
        self.branches = [] 
        for pt in self.jetPts:
            for var in self.systsJEC:
                for sel in self.jetSel:
                    self.branches.append( 'n%s%d_Recl'%(sel, int(pt))           + self.systsJEC[var])
                    self.branches.append( 'nB%sMedium%d_Recl'%(sel, int(pt)) + self.systsJEC[var])
                    self.branches.append( 'nB%sLoose%d_Recl'%(sel, int(pt))  + self.systsJEC[var])
                if self.doJetSums:
                    self.branches.append( 'ht%d_Recl'%(int(pt))  + self.systsJEC[var])
                    self.branches.append( 'mht%d_Recl'%(int(pt))  + self.systsJEC[var])

        print self.branches[:]
        return self.branches[:]
    
    def __call__(self,event):
        allret = {} 
        year = '%d'%event.year
        jets = [ j for j in Collection(event,"JetSel","nJetSel")]
        
        for br in self.branches:
            allret[br] = 0

        if self.doJetSums:
            _mht  = {}; _ht = {} 
            leps = [l for l in Collection(event, "LepFO","nLepFO")]
            taus = [t for t in Collection(event, "TausFO","nTausFO")]
            for pt in self.jetPts:
                for var in self.systsJEC:
                    _mht['%s%d'%(var,int(pt))] = ROOT.TLorentzVector()
                    _ht ['%s%d'%(var,int(pt))] = 0
                    for l in leps+taus: 
                        vl = ROOT.TLorentzVector(); vl.SetPtEtaPhiM( l.pt, 0, l.phi, 0)
                        _mht['%s%d'%(var,int(pt))] = _mht['%s%d'%(var,int(pt))] - vl
                        _ht ['%s%d'%(var,int(pt))] = _ht ['%s%d'%(var,int(pt))] + l.pt

        for pt in self.jetPts:
            for var in self.systsJEC:
                for jet in jets:
                    for sel in self.jetSel:
                        if not self.jetSel[sel](jet): continue
                        if hasattr(jet, 'pt%s'%self.systsJEC[var]) and getattr(jet,'pt%s'%self.systsJEC[var]) > pt or (not var and not hasattr(jet, 'pt%s'%self.systsJEC[var]) and jet.pt > pt): 
                            namJ = 'n%s%d_Recl'%(sel, int(pt)) + self.systsJEC[var]
                            namb = 'nB%sLoose%d_Recl'%(sel, int(pt)) + self.systsJEC[var]
                            namB = 'nB%sMedium%d_Recl'%(sel, int(pt)) + self.systsJEC[var]
                            allret[namJ] = allret[namJ] + 1
                            print getattr(jet, self.btag)
                            if getattr(jet, self.btag) > self.btagL_thr[year]: allret[namb] = allret[namb] + 1
                            if getattr(jet, self.btag) > self.btagM_thr[year]: allret[namB] = allret[namB] + 1 
                    if abs(jet.eta) > 2.4: continue 
                    if hasattr(jet, 'pt%s'%self.systsJEC[var]) and getattr(jet,'pt%s'%self.systsJEC[var]) > pt or (not var and not hasattr(jet, 'pt%s'%self.systsJEC[var]) and jet.pt > pt): 
                        vj = ROOT.TLorentzVector(); vj.SetPtEtaPhiM( jet.pt, 0, jet.phi, 0)
                        _mht['%s%d'%(var,int(pt))] = _mht['%s%d'%(var,int(pt))] - vj
                        _ht ['%s%d'%(var,int(pt))] = _ht ['%s%d'%(var,int(pt))] + jet.pt

                allret['ht%d_Recl'%(int(pt))  + self.systsJEC[var]]  = _ht ['%s%d'%(var,int(pt))]
                allret['mht%d_Recl'%(int(pt))  + self.systsJEC[var]] = _mht['%s%d'%(var,int(pt))].Pt()
        return allret 
                    
