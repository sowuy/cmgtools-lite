from CMGTools.TTHAnalysis.treeReAnalyzer import *
from PhysicsTools.HeppyCore.utils.deltar import matchObjectCollection3
import ROOT

class MyVarProxy:
    def __init__(self,lep):
        self._ob = lep
    def __getitem__(self,name):
        return self.__getattr__(name)
    def __getattr__(self,name):
        if name in self.__dict__: return self.__dict__[name]
        else: return getattr(self._ob,name)
    def eta(self): return self._ob.eta
    def phi(self): return self._ob.phi
    def pt(self): return self._ob.pt
    def pdgId(self): return self._ob.pdgId

class TTVVariables:

    def __init__(self,label,leptonSel,jetSel):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.lepSel = leptonSel
        self.jetSel = jetSel
        self.systsJEC = {0:""}#, 1:"_jecUp", -1:"_jecDown"}
        self.debugprinted = False



    def listBranches(self):
        label = self.label
        biglist = [
            ("mZ"+label,"F"),
            ("pTZ"+label,"F"),
            ]
        for  var in self.systsJEC:
            biglist.extend(
                [('mTWLep'+label+self.systsJEC[var],"F"),
                 ('pTWLep'+label+self.systsJEC[var],"F"),
                 ('minDRWLepBJet'+label+self.systsJEC[var],"F"),
                 
                 ('mWHad1'+self.systsJEC[var]+label,"F"),
                 ('pTWHad1'+self.systsJEC[var]+label,"F"),
                 ('minDrWHadBJet1'+self.systsJEC[var]+label,"F"),
                 ('mWHad2'+self.systsJEC[var]+label,"F"),
                 ('pTWHad2'+self.systsJEC[var]+label,"F"),
                 ('minDrWHadBJet2'+self.systsJEC[var]+label,"F"),

                 #tops
                 ('mTHad'+self.systsJEC[var]+label,"F"),
                 ('pTTHad'+self.systsJEC[var]+label,"F"),

                 ('mTSLep'+self.systsJEC[var]+label,"F"),
                 ('pTSLep'+self.systsJEC[var]+label,"F"),

                 ])

        return biglist

    def buildBestZCandidate(self,leps,cut=lambda lep:True):
        pairs = []
        for i1,l1 in enumerate(leps):
            if not cut(l1): continue
            for i2,l2 in enumerate(leps):
                if i2<=i1: continue
                if not cut(l2): continue
                if l1.pdgId == -l2.pdgId:
                    zCand=(l1.p4() + l2.p4())
                    mz = zCand.M()
                    diff = abs(mz-91)
                    pairs.append( (diff, zCand, (l1,l2)) )
        if len(pairs):
            pairs.sort()
            #print "===<<>>",pairs
            return pairs[0][1], pairs[0][2]
        return None,[]


    def buildBestLepWCandidate(self, leps, met, metphi, jets, vetoLeps=[], cut=lambda lep:True, jetcut=lambda:True, mW=80.4):
        ws = []
        #print "-->> ",len(vetoLeps),len(leps)
        for i,l in enumerate(leps):
            ok=True
            for vl in vetoLeps:
                if l==vl:
                    ok=False
                    break
            if not ok: continue
            #if l==vetoLeps[0] or l==vetoLeps[1]: continue
            dPhi=deltaPhi(metphi,l.phi)
            mass=sqrt(2*l.pt*met*(1-cos(dPhi)))
            pt=sqrt( (l.pt+met*cos(dPhi))**2 + (met*sin(dPhi))**2 )
            phi=acos( (l.pt+met*cos(dPhi))/pt )+l.phi
            W=ROOT.TLorentzVector()          
            W.SetPtEtaPhiM(pt,0,phi,mass)
            diff=abs(mass-mW)
            ws.append( (diff,W) )
        if len(ws):
            ws.sort()
            #print "===<<<>>>",ws[0][1].M()
            tmpTops=[]
            for j in jets:
                #print " ======>>>  blou ", jets
                if not jetcut(j): continue
                #if j==w[2][1] or j==w[2][2]: continue
                topCand=(ws[0][1]+j.p4())
                diff=abs(topCand.M()-172)
                tmpTops.append( (diff,topCand) )
            if len(tmpTops):
                tmpTops.sort()
            else:
                return ws[0][1], None
                
            return ws[0][1], tmpTops[0][1]



        return None, None

    def buildBestHadWCandidate(self, jets, cut=lambda jet:True, mW=80.4):
        pairs = []
        for i1,j1 in enumerate(jets):
            #print "aqui"
            if not cut(j1): continue
            #print "aqua"
            for i2,j2 in enumerate(jets):
                if i2<=i1: continue
                if not cut(j2): continue
                wCand=(j1.p4() + j2.p4())
                mw = wCand.M()
                diff = abs(mw-mW)
                pairs.append( (diff, wCand, (j1,j2)) )
        #print "grou "
        #print "---------------> ", len(jets), len(pairs)
        if len(pairs):
            pairs.sort()
            #print pairs
            bestdr=[]
            besttops=[]
            #print " coin "
            #jets.sort(key=lambda j:1-j.btagCSV)
            for w in pairs:
                bestdr.append( min( [(deltaR(j.eta,w[1].Eta(),j.phi,w[1].Phi()) if (j!=w[2][0] and j!=w[2][1]) else 100) for j in jets[:2]]) if len(jets)>=2 else 100)
                tmpTops=[]
                #print " ======  blou "
                for j in jets:
                    if not cut(j): continue
                    #print "aqui?", w, "--->> ", len(w[2])
                    #print "///",w[2][0]
                    #print "..",w[2][1]
                    #print "...", j
                    if j==w[2][0] or j==w[2][1]: continue
                    #print "you"
                    #print "bjbjb "
                    topCand=(w[1]+j.p4())
                    diff=abs(topCand.M()-172)
                    tmpTops.append( (diff,topCand) )
                #print "yoyoyoy"
                if len(tmpTops):
                    tmpTops.sort()
                    #print "<><>"
                    besttops.append(tmpTops[0][1])
                    
            #print "===<<-->>",pairs
            theWs=[]
            theTops=[]
            for i,k in enumerate(pairs):
                #print "<><><><>",k
                theWs.append(k[1])
                if(len(besttops)>=i+1):
                    theTops.append(besttops[i])
                #print theWs[:1], theWs[:1][0]
                #print "-->> ",theWs[-1].M(), bestdr
            return theWs, bestdr, theTops
        return [],[],[]


    def __call__(self,event):
        self.ev = event
        ret = {};
        fullret = {}
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets={}
        for var in self.systsJEC:
            jets[var] =[j for j in Collection(event,"JetSel_Recl"+self.systsJEC[var],"nJetSel_Recl"+self.systsJEC[var])]
            jets[var].sort(key=lambda j:1-j.btagCSV)            
        met={}
        metphi={}
        met[0]=event.met_pt
        met[1]=getattr(event,"met_jecUp_pt",event.met_pt)
        met[-1]=getattr(event,"met_jecDown_pt",event.met_pt)
        metphi[0]= event.met_phi
        metphi[1]= getattr(event,"met_jecUp_phi",event.met_phi)
        metphi[-1]= getattr(event,"met_jecDown_phi",event.met_phi)
        #print " test 1 ============"
        bestZCand, lepsZ = self.buildBestZCandidate(leps, self.lepSel)
        bestWLepCands={}
        bestTSLepCands={}
        #print " test 2 ============"
        for var in self.systsJEC:
            #print " yolo ",var
            bestWLepCands[var], bestTSLepCands[var] = self.buildBestLepWCandidate(leps, met[var], metphi[var], jets[var], lepsZ ,self.lepSel,self.jetSel)
        hadWCands={}
        hadTopCands={}
        drBJets={}
        #print " test 3 ============"
        for var in self.systsJEC:
            hadWCands[var], drBJets[var], hadTopCands[var]=self.buildBestHadWCandidate(jets[var],self.jetSel)
            #print " <>>>> ", hadWCands[var], drBJets[var]

        #print " test 4 ============"
        ret['mZ']=bestZCand.M() if bestZCand!=None else -1
        ret['pTZ']=bestZCand.Pt() if bestZCand!=None else -1
        #print " test 5 ============"
        for var in self.systsJEC:
            ret['mTWLep'+self.systsJEC[var]]=bestWLepCands[var].M() if bestWLepCands[var]!=None else -1
            ret['pTWLep'+self.systsJEC[var]]=bestWLepCands[var].Pt() if bestWLepCands[var]!=None else -1
            #print "cui ", len(jets[var])
            ret['minDRWLepBJet'+self.systsJEC[var]]=min([deltaR(j.eta,bestWLepCands[var].Eta(),j.phi,bestWLepCands[var].Phi()) for j in jets[var][:2]]) if (bestWLepCands[var]!=None and len(jets[var])!=0) else 100
            #print "cui1 "
            #print " test 6 ============"
            ret['mWHad1'+self.systsJEC[var]]=-1
            ret['pTWHad1'+self.systsJEC[var]]=-1
            ret['minDrWHadBJet1'+self.systsJEC[var]]=100
            ret['mWHad2'+self.systsJEC[var]]=-1
            ret['pTWHad2'+self.systsJEC[var]]=-1
            ret['minDrWHadBJet2'+self.systsJEC[var]]=100
            for i,w in enumerate(hadWCands[var]):
                if i>=2: break
                #print " test 7 ============",i, w
                ret['mWHad'+str(i+1)+self.systsJEC[var]]=w.M()
                ret['pTWHad'+str(i+1)+self.systsJEC[var]]=w.Pt()
                ret['minDrWHadBJet'+str(i+1)+self.systsJEC[var]]=drBJets[var][i]

            ret['mTHad'+self.systsJEC[var]]=hadTopCands[var][0].M() if len(hadTopCands[var]) else -100
            ret['pTTHad'+self.systsJEC[var]]=hadTopCands[var][0].Pt() if len(hadTopCands[var]) else -100

            ret['mTSLep'+self.systsJEC[var]]=bestTSLepCands[var].M() if bestTSLepCands[var]!=None else -100
            ret['pTSLep'+self.systsJEC[var]]=bestTSLepCands[var].Pt() if bestTSLepCands[var]!=None else -100

        for k,v in ret.iteritems(): 
            fullret[k+self.label] = v
        return fullret



## =========================== End class ====================================
