def _ttH_idEmu_cuts_E2(lep):
    if (abs(lep.pdgId)!=11): return True
    if (lep.hadronicOverEm>=(0.10-0.03*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dEtaScTrkIn)>=(0.01-0.002*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dPhiScTrkIn)>=(0.04+0.03*(abs(lep.etaSc)>1.479))): return False
    if (lep.eInvMinusPInv<=-0.05): return False
    if (lep.eInvMinusPInv>=(0.01-0.005*(abs(lep.etaSc)>1.479))): return False
    if (lep.sigmaIEtaIEta>=(0.011+0.019*(abs(lep.etaSc)>1.479))): return False
    return True
def _ttH_idEmu_cuts_E2_obj(lep):
    if (abs(lep.pdgId())!=11): return True
    etasc = lep.superCluster().eta()
    if (lep.hadronicOverEm()>=(0.10-0.03*(abs(etasc)>1.479))): return False
    if (abs(lep.deltaEtaSuperClusterTrackAtVtx())>=(0.01-0.002*(abs(etasc)>1.479))): return False
    if (abs(lep.deltaPhiSuperClusterTrackAtVtx())>=(0.04+0.03*(abs(etasc)>1.479))): return False
    eInvMinusPInv = (1.0/lep.ecalEnergy() - lep.eSuperClusterOverP()/lep.ecalEnergy()) if lep.ecalEnergy()>0. else 9e9
    if (eInvMinusPInv<=-0.05): return False
    if (eInvMinusPInv>=(0.01-0.005*(abs(etasc)>1.479))): return False
    if (lep.full5x5_sigmaIetaIeta()>=(0.011+0.019*(abs(etasc)>1.479))): return False
    return True

def _ttH_idEmu_cuts_E3(lep):
    if (abs(lep.pdgId)!=11): return True
    if (lep.hoe>=(0.10-0.00*(abs(lep.deltaEtaSC+lep.eta)>1.479))): return False
    if (lep.eInvMinusPInv<=-0.04): return False
    if (lep.sieie>=(0.011+0.019*(abs(lep.deltaEtaSC+lep.eta)>1.479))): return False
    return True
# def _ttH_idEmu_cuts_E3_obj(lep):
#     if (abs(lep.pdgId())!=11): return True
#     etasc = lep.superCluster().eta()
#     if (lep.hadronicOverEm()>=(0.10-0.00*(abs(etasc)>1.479))): return False
#     eInvMinusPInv = (1.0/lep.ecalEnergy() - lep.eSuperClusterOverP()/lep.ecalEnergy()) if lep.ecalEnergy()>0. else 9e9
#     if (eInvMinusPInv<=-0.04): return False
#     if (lep.full5x5_sigmaIetaIeta()>=(0.011+0.019*(abs(etasc)>1.479))): return False
#     return True

def _soft_MuonId_2016ICHEP(lep):
    if (abs(lep.pdgId())!=13): return False
    if not lep.muonID("TMOneStationTight"): return False #TMOneStationTightMuonId
    if not lep.track().hitPattern().trackerLayersWithMeasurement() > 5: return False
    if not lep.track().hitPattern().pixelLayersWithMeasurement() > 0: return False
    if not (abs(lep.dxy())<0.3 and abs(lep.dz())<20): return False
    return True

def _medium_MuonId_2016ICHEP(lep):
    if (abs(lep.pdgId())!=13): return False
    if not (lep.physObj.isGlobalMuon() or lep.physObj.isTrackerMuon()): return False
    if not (lep.innerTrack().validFraction()>0.49): return False
    if lep.segmentCompatibility()>0.451: return True
    else:
        if not lep.globalTrack().isNonnull(): return False
        if not lep.isGlobalMuon: return False
        if not lep.globalTrack().normalizedChi2()<3: return False
        if not lep.combinedQuality().chi2LocalPosition<12: return False
        if not lep.combinedQuality().trkKink<20: return False 
        if not lep.segmentCompatibility()>0.303: return False

    return True


from CMGTools.TTHAnalysis.tools.leptonJetReCleaner import LeptonJetReCleaner
from CMGTools.TTHAnalysis.tools.conept import conept_TTH

MODULES=[]

from CMGTools.TTHAnalysis.tools.combinedObjectTaggerForCleaning import *
from CMGTools.TTHAnalysis.tools.fastCombinedObjectRecleaner import *
from CMGTools.TTHAnalysis.tools.objFloatCalc import ObjFloatCalc




def preselectMuon(lep):
    return lep.pt > 5 and abs(lep.eta) < 2.4 and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1 and lep.miniPFRelIso_all < 0.4 and lep.sip3d < 8

def preselectElectron(lep):
    return lep.pt > 7 and abs(lep.eta) < 2.5 and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.1 and lep.miniPFRelIso_all < 0.4  and lep.sip3d < 8 and lep.mvaFall17V1noIso_WPL and lep.lostHits <=1 

def preselectLepton(lep):
    return preselectElectron(lep) if abs(lep.pdgId) == 11 else preselectMuon(lep)


# def jetLepAwareJEC(lep,jet): 

#     l = ROOT.TLorentzVector(); l.SetPtEtaPhiM(lep.pt, lep.eta, lep.phi, lep.mass)
#     j = ROOT.TLorentzVector(); j.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass) 

#     corrFactor = 1 / (1-jet.rawFactor)

#     if ((j*(1-jet.rawFactor)-l).Rho()<1e-4): 
#         return l # matched to jet containing only the lepton

    
#     j = (j*(1-jet.rawFactor)-l*(1.0/jet.l1corrFactor))*corrFactor+l

#     return j
# #    return j.Pt(), 1-jet.rawFactor, jet.l1corrFactor, corrFactor



# def _pratio(lep,jet):
#     return lep.pt / jetLepAwareJEC(lep, jet).Pt() 

# def _ptrel(lep,jet):
#     m = jetLepAwareJEC(lep,jet)
#     l = ROOT.TLorentzVector(); l.SetPtEtaPhiM(lep.pt, lep.eta, lep.phi, lep.mass)
#     if ((m-l).Rho() < 1e-4): return 0 
#     return l.Perp((m-l).Vect())

 






from CMGTools.TTHAnalysis.tools.lepJetVars import LepJetVars
MODULES.append( ('lepJetVars', lambda : LepJetVars(None, 
                                                   dict( jetBTagDeepCSV = (lambda lep, jet: jet.btagDeepB, lambda lep : 0 )))))




from CMGTools.TTHAnalysis.tools.genLepVars import GenLepVars

MODULES.append( ('genLepVars', lambda : GenLepVars(None, 
                                                    dict( isMatchRightCharge = (lambda lep, gen: (lep.genPartFlav == 1 or lep.genPartFlav == 15) and (gen.pdgId == lep.pdgId), lambda lep : 0),
                                                          )
)))






def _bitFromInt(num, idx):
    # returns the bit in the idx's position of num 
    bitMap = "{0:b}".format(num)
    if idx > len(bitMap): return False
    return bool(int(bitMap[-idx]))


def clean_and_FO_selection_TTH(lep):
    return lep.conept>10 and lep.jetBTagDeepCSV<0.4941 and (abs(lep.pdgId)!=11 or (_ttH_idEmu_cuts_E3(lep) and lep.convVeto and lep.lostHits==0) \
    ) \
        and (lep.mvaTTH>0.90 or \
                 (abs(lep.pdgId)==13 and lep.jetBTagDeepCSV<0.07 and lep.segmentComp>0.3 and 1/(1+lep.jetRelIso)>0.60) or \
                 (abs(lep.pdgId)==11 and lep.jetBTagDeepCSV<0.07 and lep.mvaFall17V1noIso>0.5 and 1/(1+lep.jetRelIso)>0.60) \
                 )

def _FOTauSel(tau):
    return tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and _bitFromInt(tau.idMVAoldDMdR032017v2,2) and tau.idDecayMode

# MODULES.append( ('leptonJetFastReCleanerTTH_step1', lambda : CombinedObjectTaggerForCleaning("InternalRecl",
#                                                                                        looseLeptonSel = lambda lep : preselectLepton(lep),
#                                                                                        cleaningLeptonSel = clean_and_FO_selection_TTH,
#                                                                                        FOLeptonSel = clean_and_FO_selection_TTH,
#                                                                                        tightLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep) and (abs(lep.pdgId)!=13 or lep.mediumId>0) and lep.mvaTTH > 0.90,
#                                                                                        FOTauSel = lambda tau: _FOTauSel(tau), #  
#                                                                                        tightTauSel = lambda tau: _bitFromInt(tau.idMVAoldDMdR032017v2,3), 
#                                                                                        selectJet = lambda jet: abs(jet.eta)<2.4 and _bitFromInt(jet.jetId,2) and (jet.pt) > 15,
#                                                                                        coneptdef = lambda lep: conept_TTH(lep) ) ))
# MODULES.append( ('leptonJetFastReCleanerTTH_step2_mc',lambda : fastCombinedObjectRecleaner(label="Recl",
#                                                                                            inlabel="_InternalRecl",
#                                                                                            cleanTausWithLooseLeptons=True,
#                                                                                            cleanJetsWithFOTaus=True,
#                                                                                            doVetoZ=False,
#                                                                                            doVetoLMf=False,
#                                                                                            doVetoLMt=False,
#                                                                                            jetPts=[25,40],
#                                                                                            btagL_thr=0.1522,
#                                                                                            btagM_thr=0.4941,
#                                                                                            isMC = True) ))
# MODULES.append( ('leptonJetFastReCleanerTTH_step2_data',lambda : fastCombinedObjectRecleaner(label="Recl",
#                                                                                              inlabel="_InternalRecl",
#                                                                                              cleanTausWithLooseLeptons=True,
#                                                                                              cleanJetsWithFOTaus=True,
#                                                                                              doVetoZ=False,
#                                                                                              doVetoLMf=False,
#                                                                                              doVetoLMt=False,
#                                                                                              jetPts=[25,40],
#                                                                                              btagL_thr=0.1522,
#                                                                                              btagM_thr=0.4941,
#                                                                                              isMC = False) ))







# 1_recleaner: 
# run mc: lepJetVars genLepVars leptonJetFastReCleanerTTH_step1 leptonJetFastReCleanerTTH_step2_mc
# run data:  lepJetVars  leptonJetFastReCleanerTTH_step1 leptonJetFastReCleanerTTH_step2_data

# 2_TauTightFlag:
# TauTightFlag


from CMGTools.TTHAnalysis.tools.countJets import CountJets
MODULES.append( ('countJets',lambda : CountJets(jetPts=[25,40],
                                                jetSel={ "JetCentral" : lambda x : abs(x.eta) < 2.4,
                                                         "JetForward" : lambda x : abs(x.eta) > 2.4 } ,
                                                btagL_thr=0.1522,
                                                btagM_thr=0.4941,
                                                doJetSums=True,
                                                )))

# 4_eventVars: eventVars


from CMGTools.TTHAnalysis.tools.eventVars_2lss import EventVars2LSS
MODULES.append( ('eventVars', lambda : EventVars2LSS('','Recl')) )



from CMGTools.TTHAnalysis.tools.kinMVA_2D_2lss_3l import KinMVA_2D_2lss_3l
MODULES.append( ('kinMVA_2D_2lss_3l', lambda : KinMVA_2D_2lss_3l(os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/data/kinMVA/tth/%s_BDTG.weights.xml", useTT_2lss='v8,rTT,httTT', useMEM_3l = False)) )
MODULES.append( ('noTTMVA_2D_2lss_3l', lambda : KinMVA_2D_2lss_3l(os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/data/kinMVA/tth/%s_BDTG.weights.xml", useTT_2lss='', useMEM_3l = False)) )
MODULES.append( ('kinMEMMVA_2D_2lss_3l', lambda : KinMVA_2D_2lss_3l(os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/data/kinMVA/tth/%s_BDTG.weights.xml", useTT_2lss='v8,rTT,httTT', useMEM_3l = True)) )

from CMGTools.TTHAnalysis.tools.BDTv8_eventReco_cpp import BDTv8_eventReco
MODULES.append( ('oldcode_BDTv8_Hj', lambda : BDTv8_eventReco(os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_bloose_BDTG.weights.xml',
                                                      os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_btight_BDTG.weights.xml',
                                                      os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hj_csv_BDTG.weights.xml',
                                                      os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hjj_csv_BDTG.weights.xml',
                                                      selection = [
                lambda leps,jets,event : len(leps)>=2 and len(jets)>=3,
                lambda leps,jets,event : leps[0].conePt>20 and leps[1].conePt>10,
                ]
                                                      )) )

from CMGTools.TTHAnalysis.tools.BDT_eventReco_cpp import BDT_eventReco
MODULES.append( ('BDTv8_Hj', lambda : BDT_eventReco(os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_bloose_BDTG.weights.xml',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_btight_BDTG.weights.xml',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hj_2017_configA_dcsv_BDTG.weights.xml',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hjj_csv_BDTG.weights.xml',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/resTop_xgb_csv_order_deepCTag.xml.gz',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml.gz',
                                                    os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TF_jets_kinfit_httTT.root',
                                                    algostring = 'k_BDTv8_Hj',
                                                    csv_looseWP = 0.5426,
                                                    csv_mediumWP = 0.8484,
                                                    selection = [
                lambda leps,jets,event : len(leps)>=2 and len(jets)>=3,
                lambda leps,jets,event : leps[0].conePt>20 and leps[1].conePt>10,
                ]
                                                            )) )
MODULES.append( ('BDTrTT_Hj', lambda : BDT_eventReco(os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_bloose_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_btight_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hj_2017_configA_dcsv_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hjj_csv_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/resTop_xgb_csv_order_deepCTag.xml.gz',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml.gz',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TF_jets_kinfit_httTT.root',
                                                     algostring = 'k_rTT_Hj',
                                                     csv_looseWP = 0.5426,
                                                     csv_mediumWP = 0.8484,
                                                      selection = [
                lambda leps,jets,event : len(leps)>=2 and len(jets)>=3,
                lambda leps,jets,event : leps[0].conePt>20 and leps[1].conePt>10,
                ]
                                                     )) )
MODULES.append( ('BDThttTT_Hj', lambda : BDT_eventReco(os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_bloose_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TMVAClassification_btight_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hj_2017_configA_dcsv_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/Hjj_csv_BDTG.weights.xml',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/resTop_xgb_csv_order_deepCTag.xml.gz',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml.gz',
                                                     os.environ["CMSSW_BASE"]+'/src/CMGTools/TTHAnalysis/data/kinMVA/tth/TF_jets_kinfit_httTT.root',
                                                     algostring = 'k_httTT_Hj',
                                                     csv_looseWP = 0.5426,
                                                     csv_mediumWP = 0.8484,
                                                      selection = [
                lambda leps,jets,event : len(leps)>=2 and len(jets)>=3,
                lambda leps,jets,event : leps[0].conePt>20 and leps[1].conePt>10,
                ]
                                                     )) )

from CMGTools.TTHAnalysis.tools.evtTagger import EvtTagger

MODULES.append( ('isData', lambda : EvtTagger("isData",[
                lambda ev : not hasattr(ev,'xsec')
                    ])))


MODULES.append( ('Trigger_1e', lambda : EvtTagger("Trigger_1e",[
                lambda ev : (ev.HLT_Ele32_WPTight_Gsf if hasattr(ev, 'HLT_Ele32_WPTight_Gsf')  else False) or ev.HLT_Ele35_WPTight_Gsf
                    ])))
MODULES.append( ('Trigger_1m', lambda : EvtTagger("Trigger_1m",[
                lambda ev : ev.HLT_IsoMu24 or ev.HLT_IsoMu27
                    ])))
MODULES.append( ('Trigger_2e', lambda : EvtTagger("Trigger_2e",[
                lambda ev : ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or ev.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
                    ])))
MODULES.append( ('Trigger_2m', lambda : EvtTagger("Trigger_2m",[
                lambda ev : ev.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or ( ev.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 if hasattr(ev, 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8') else False)
                    ])))
MODULES.append( ('Trigger_em', lambda : EvtTagger("Trigger_em",[
                lambda ev : (ev.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL if hasattr(ev,'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL') else False) or \
                    ev.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or \
                    ev.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
                    ])))
MODULES.append( ('Trigger_3e', lambda : EvtTagger("Trigger_3e",[
                lambda ev : ev.HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL
                    ])))
MODULES.append( ('Trigger_3m', lambda : EvtTagger("Trigger_3m",[
                lambda ev : ev.HLT_TripleMu_12_10_5
                    ])))
MODULES.append( ('Trigger_mee', lambda : EvtTagger("Trigger_mee",[
                lambda ev : ev.HLT_Mu8_DiEle12_CaloIdL_TrackIdL
                    ])))
MODULES.append( ('Trigger_mme', lambda : EvtTagger("Trigger_mme",[
                lambda ev : ev.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ
                    ])))
MODULES.append( ('Trigger_2lss', lambda : EvtTagger("Trigger_2lss",[
                lambda ev : ev.Trigger_1e or ev.Trigger_1m or ev.Trigger_2e or ev.Trigger_2m or ev.Trigger_em ])))
MODULES.append( ('Trigger_3l', lambda : EvtTagger("Trigger_3l",[
                lambda ev : ev.Trigger_2lss or ev.Trigger_3e or ev.Trigger_3m or ev.Trigger_mee or ev.Trigger_mme ])))

# -m 'Trigger_1e' -m 'Trigger_1m' -m 'Trigger_2e' -m 'Trigger_2m' -m 'Trigger_em' -m 'Trigger_3e' -m 'Trigger_3m' -m 'Trigger_mee' -m 'Trigger_mme' -m 'Trigger_2lss' -m 'Trigger_3l'
    

from CMGTools.TTHAnalysis.tools.objTagger import ObjTagger
MODULES.append( ('TauTightFlag', lambda : ObjTagger("isTauTight","TausFO",
                                                    [lambda tau : _bitFromInt(tau.idMVAoldDMdR032017v2,3)]
                                                )))

MODULES.append( ('LepTightFlag', lambda : ObjTagger("isTight","LepFO",
                                                    [lambda lep : (abs(lep.pdgId)!=13 or lep.mediumId>0) and lep.mvaTTH > 0.90]
                                                    )))

MODULES.append( ('LepConePt', lambda : ObjFloatCalc(None,"LepFO",
                                                    dict(conePt = lambda lep: conept_TTH(lep))
                                                    )))
from CMGTools.TTHAnalysis.tools.massCalculator import MassCalculator
MODULES.append( ('MassCalculator', lambda : MassCalculator("LepFO","LepLoose")))


# -m 'LepTightFlag' -m 'LepConePt' -m 'MassCalculator'
                                                    

from CMGTools.TTHAnalysis.tools.bTagEventWeightsCSVFullShape import BTagEventWeightFriend
MODULES.append( ('eventBTagWeight', lambda : BTagEventWeightFriend(csvfile=os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/data/btag/DeepCSV_94XSF_V4_B_F.csv",
                                                                   discrname="btagDeepB")))

from CMGTools.TTHAnalysis.tools.higgsRecoTTH import HiggsRecoTTH
MODULES.append( ('higgsRecoTTH', lambda : HiggsRecoTTH(label="_Recl",
                                                       cut_BDT_rTT_score = 0.0,
                                                       cuts_mW_had = (60.,100.),
                                                       cuts_mH_vis = (80.,140.),
                                                       btagDeepCSVveto = 0.1522) ))

from CMGTools.TTHAnalysis.tools.ttHMCEventReco import TTHMCEventReco
MODULES.append( ('genLevelChain', lambda : TTHMCEventReco()) )

from CMGTools.TTHAnalysis.tools.matchRecoToPartonsTTH import MatchRecoToPartonsTTH
MODULES.append( ('matchPartons', lambda : MatchRecoToPartonsTTH(label="_Recl")) )

from CMGTools.TTHAnalysis.tools.vertexWeightFriend import VertexWeightFriend
# run on a big number of events, ideally one job per file using -N big_number, and not on skimmed trees when auto-reweighthing the pileup to avoid loss of statistical power!
MODULES.append( ('vtxWeight', lambda : VertexWeightFriend(myfile=None, targetfile=os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/data/pileup/puWeights_2017_41p4fb_rereco_69p2mb.root",
                                                          myhist=None,targethist="pileup",name="vtxWeight2017",
                                                          verbose=False,vtx_coll_to_reweight="nTrueInt",autoPU=True)) )

from CMGTools.TTHAnalysis.tools.bestHmmFriend import BestHmm
MODULES.append( ('bestHmm', lambda : BestHmm(label="_Recl")) )

from CMGTools.TTHAnalysis.tools.synchNtuples import SynchNtuples
MODULES.append( ('synchNtuples', lambda : SynchNtuples()))





