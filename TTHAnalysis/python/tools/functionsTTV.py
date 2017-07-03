
from CMGTools.TTHAnalysis.tools.leptonJetReCleaner import LeptonJetReCleaner
from CMGTools.TTHAnalysis.tools.conept import conept_TTH

MODULES=[]

from CMGTools.TTHAnalysis.tools.combinedObjectTaggerForCleaning import *
from CMGTools.TTHAnalysis.tools.fastCombinedObjectRecleaner import *

from CMGTools.TTHAnalysis.tools.functionsTTH import _ttH_idEmu_cuts_E2, clean_and_FO_selection_TTH
from CMGTools.TTHAnalysis.tools.TTVVariables import TTVVariables

MODULES.append( ('leptonJetReCleanerTTV', lambda : LeptonJetReCleaner("Recl",
                                                                      looseLeptonSel = lambda lep : lep.miniRelIso < 0.4 and lep.sip3d < 8,
                                                                      cleaningLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep),
                                                                      FOLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep),
                                                                      tightLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep) and (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90,
                                                                      cleanJet = lambda lep,jet,dr : dr<0.4,
                                                                      selectJet = lambda jet: abs(jet.eta)<5 and jet.pt>25,
# and (abs(jet.eta)<3 or jet.puId)
                                                                      cleanTau = lambda tau: tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      looseTau = lambda tau: tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      tightTau = lambda tau: tau.idMVAdR03 >= 3,
                                                                      cleanJetsWithTaus = False,
                                                                      cleanTausWithLoose = False,
                                                                      doVetoZ = False,
                                                                      doVetoLMf = False,
                                                                      doVetoLMt = False,
                                                                      jetPt = 25,
                                                                      bJetPt = 25,
                                                                      coneptdef = lambda lep: conept_TTH(lep),
                                                                      storeJetVariables = True) ))



MODULES.append( ('ttZVariables', lambda: TTVVariables(label="Extra",
                                                      leptonSel=lambda lep: lep.isTight_Recl,
                                                      jetSel=lambda jet:True
                                                      ) ))
