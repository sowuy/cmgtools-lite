import re, os, sys
from CMGTools.RootTools.samples.configTools import printSummary, mergeExtensions, doTestN, configureSplittingFromTime, cropToLumi
from CMGTools.RootTools.samples.autoAAAconfig import autoAAA
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()
def byCompName(components, regexps):
    return [ c for c in components if any(re.match(r, c.name) for r in regexps) ]

year = int(getHeppyOption("year", "2018"))
analysis = getHeppyOption("analysis", "main")
preprocessor = getHeppyOption("nanoPreProcessor")

if preprocessor:
    if year == 2018:
        from CMGTools.RootTools.samples.samples_13TeV_RunIIAutumn18MiniAOD import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2018_MiniAOD import samples as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2018 import all_triggers as triggers
    elif year == 2017:
        from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17MiniAOD import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2017 import dataSamples_31Mar2018 as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import all_triggers as triggers
    elif year == 2016:
        from CMGTools.RootTools.samples.samples_13TeV_RunIISummer16MiniAODv3 import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2016 import dataSamples_17Jul2018 as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2016 import all_triggers as triggers
else:
    if year == 2018:
        from CMGTools.RootTools.samples.samples_13TeV_RunIIAutumn18NanoAODv4 import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2018_NanoAOD import samples as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2018 import all_triggers as triggers
    elif year == 2017:
        from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17NanoAODv4 import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2017_NanoAOD import samples as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import all_triggers as triggers
    elif year == 2016:
        from CMGTools.RootTools.samples.samples_13TeV_RunIISummer16NanoAODv4 import samples as mcSamples_
        from CMGTools.RootTools.samples.samples_13TeV_DATA2016_NanoAOD import samples as allData
        from CMGTools.RootTools.samples.triggers_13TeV_DATA2016 import all_triggers as triggers

flavor = getHeppyOption("flavor")
DatasetsAndTriggers = []
if year == 2018:
    mcSamples = byCompName(mcSamples_, [
        "DYJetsToLL_M50$",
    ])
    if flavor == "Muon":
        DatasetsAndTriggers.append( ("SingleMuon", triggers["1mu_iso"]) )
    elif flavor == "Electron":
        DatasetsAndTriggers.append( ("EGamma", triggers["1e_iso"]) )

elif year == 2017:
    mcSamples = byCompName(mcSamples_, [
        "DYJetsToLL_M50$",
    ])
    if flavor == "Muon":
        DatasetsAndTriggers.append( ("SingleMuon", triggers["1mu_iso"]) )
    elif flavor == "Electron":
        DatasetsAndTriggers.append( ("SingleElectron", triggers["1e_iso"]) )


elif year == 2016:
    mcSamples = byCompName(mcSamples_, [
        "DYJetsToLL_M50$",
    ])
    if flavor == "Muon":
        DatasetsAndTriggers.append( ("SingleMuon", triggers["1mu_iso"]) )
    elif flavor == "Electron":
        DatasetsAndTriggers.append( ("SingleElectron", triggers["1e_iso"]) )


# make MC
mcTriggers = sum((trigs for (pd,trigs) in DatasetsAndTriggers), [])
for comp in mcSamples:
    comp.triggers = mcTriggers

# make data
dataSamples = []; vetoTriggers = []
for pd, triggers in DatasetsAndTriggers:
    for comp in byCompName(allData, [pd]):
        comp.triggers = triggers[:]
        comp.vetoTriggers = vetoTriggers[:]
        dataSamples.append(comp)
    vetoTriggers += triggers[:]

selectedComponents = mcSamples + dataSamples
if getHeppyOption('selectComponents'):
    selectedComponents = byCompName(selectedComponents, getHeppyOption('selectComponents').split(","))
autoAAA(selectedComponents, quiet=not(getHeppyOption("verboseAAA",False)))
configureSplittingFromTime(selectedComponents,100 if preprocessor else 10,4)
selectedComponents, _ = mergeExtensions(selectedComponents)

# create and set preprocessor if requested
if preprocessor:
    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][0]),cmsswArea=preproc_cmsswArea,keepOutput=True)
    if year==2018:
        preproc_data_ABC = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][1]),cmsswArea=preproc_cmsswArea,keepOutput=True)
        preproc_data_D = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][2]),cmsswArea=preproc_cmsswArea,keepOutput=True)
        for comp in selectedComponents:
            if comp.isData:
                comp.preprocessor = preproc_data_D if '2018D' in comp.name else preproc_data_ABC
            else:
                comp.preprocessor = preproc_mc
    else:
        preproc_data = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][1]),cmsswArea=preproc_cmsswArea,keepOutput=True)
        for comp in selectedComponents:
            comp.preprocessor = preproc_data if comp.isData else preproc_mc

# print summary of components to process
if getHeppyOption("justSummary"): 
    printSummary(selectedComponents)
    sys.exit(0)

from CMGTools.TTHAnalysis.tools.nanoAOD.addTnpTree import addTnpTree

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor


# in the cut string, keep only the main cuts to have it simpler
from CMGTools.TTHAnalysis.tools.nanoAOD.autoPuWeight import autoPuWeight
modules = [ autoPuWeight, addTnpTree(year,flavor) ] 
cut = 'n%s > 1'%flavor

branchsel_in = os.environ['CMSSW_BASE']+"/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/branchsel_in.txt"
branchsel_out = os.environ['CMSSW_BASE']+"/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/branchsel_out.txt"
compression = "ZLIB:3" #"LZ4:4" #"LZMA:9"

POSTPROCESSOR = PostProcessor(None, [], modules = modules,
        cut = cut, prefetch = True, longTermCache = True,
        branchsel = branchsel_in, outputbranchsel = branchsel_out, compression = compression)

test = getHeppyOption("test")
if test == "synch-2016":
    TTLep_pow = kreator.makeMCComponent("TTHnobb_pow", "/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 0.5085*(1-0.577) )
    TTLep_pow.files = ["/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/CMSSW_10_4_0/src/CMGTools/TTHAnalysis/cfg/F24F2D5E-DDEC-E811-AF50-90B11C08AD7D.root"]

    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][0]),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc
    selectedComponents = [TTLep_pow]

if test == "synch-2016-data":
    TTLep_pow = kreator.makeDataComponent("SingleMuon_2016H", "/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD", "CMS", ".*root", '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' )
    TTLep_pow.files = ["root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/00000/68B4A70D-998C-E811-816E-AC1F6B23C82E.root"]

    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/data94X2016_NANO.py'%(preproc_cmsswArea),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc
    selectedComponents = [TTLep_pow]

if test == "synch-2017-data":
    TTLep_pow = kreator.makeDataComponent("Electron_2017B", "/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD", "CMS", ".*root", '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt' )
    TTLep_pow.files = ["root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/40BCA89D-2C38-E811-9682-008CFAC93CF8.root"]

    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/data94Xv2_NANO.py'%(preproc_cmsswArea),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc
    selectedComponents = [TTLep_pow]


if test == "synch-2016-nano":
    TTLep_pow = kreator.makeMCComponent("TTHnobb_pow", "/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", "CMS", ".*root", 0.5085*(1-0.577) )
    TTLep_pow.files = ["/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/CMSSW_10_4_0/src/CMGTools/TTHAnalysis/cfg/tnp_test/TTHnobb_pow_Chunk0/cmsswPreProcessing.root"]
    selectedComponents = [TTLep_pow]
elif test == "synch-2017":
    TTLep_pow = kreator.makeMCComponent("TTHnobb_pow", "/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM", "CMS", ".*root", 0.5085*(1-0.577) )
    TTLep_pow.files = ["/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/CMSSW_10_4_0/src/CMGTools/TTHAnalysis/cfg/7C60AC2B-E76F-E811-9D60-0025905B860C.root"]

    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][0]),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc

    selectedComponents = [TTLep_pow]

elif test == "synch-2018":
    TTLep_pow = kreator.makeMCComponent("TTHnobb_pow", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", "CMS", ".*root", 0.5085*(1-0.577) )
    TTLep_pow.files = ["/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/CMSSW_10_4_0/src/CMGTools/TTHAnalysis/cfg/A912BBFA-D1A1-8544-A430-8C98C0767737.root"]


    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][0]),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc

    selectedComponents = [TTLep_pow]

elif test == "94X-MC-miniAOD":
    TTLep_pow = kreator.makeMCComponent("DYJets", "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM ", "CMS", ".*root", 1)
    TTLep_pow.files = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/120000/88F8E690-AAEA-E811-B764-0CC47A4C8E2E.root' ]
    localfile = os.path.expandvars("/tmp/$USER/%s" % os.path.basename(TTLep_pow.files[0]))
    if os.path.exists(localfile): TTLep_pow.files = [ localfile ] 
    from CMGTools.Production.nanoAODPreprocessor import nanoAODPreprocessor
    preproc_cfg = {2016: ("mc94X2016","data94X2016"),
                   2017: ("mc94Xv2","data94Xv2"),
                   2018: ("mc102X","data102X_ABC","data102X_D")}
    preproc_cmsswArea = "/afs/cern.ch/work/s/sesanche/private/TTH/NanoAOD/TnP/CMSSW_10_2_X_2019-05-23-1100/"
    preproc_mc = nanoAODPreprocessor(cfg='%s/src/PhysicsTools/NanoAOD/test/%s_NANO.py'%(preproc_cmsswArea,preproc_cfg[year][0]),cmsswArea=preproc_cmsswArea,keepOutput=True)
    TTLep_pow.preprocessor = preproc_mc
    selectedComponents = [TTLep_pow]
elif test == "102X-MC":
    TTLep_pow = kreator.makeMCComponent("TTLep_pow", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv4-Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/NANOAODSIM", "CMS", ".*root", 831.76*((3*0.108)**2), useAAA=True )
    TTLep_pow.files = TTLep_pow.files[:1]
    selectedComponents = [TTLep_pow]
elif test in ('2','3','3s'):
    doTestN(test, selectedComponents)
