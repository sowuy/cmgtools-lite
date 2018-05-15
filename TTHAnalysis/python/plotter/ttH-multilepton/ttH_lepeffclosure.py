#!/usr/bin/env python
import sys
import re
import os

ODIR=sys.argv[1]


TREES = "--Fs {P}/1_recleaner_230418_v2 --Fs {P}/5_triggerDecision_230418_v1 --Fs {P}/7_tauTightSel_v2 --FMCs {P}/8_vtxWeight2017_v1"
TREESONLYSKIM = "-P /pool/ciencias/HeppyTrees/ttH/LaEstructuraICHEP18/"
TREESONLYFULL = "-P /pool/ciencias/HeppyTrees/ttH/LaEstructuraICHEP18/"

def base():

    CORE=' '.join([TREES,TREESONLYSKIM])

    CORE+=" -f -j 40 -l 41.4 --s2v -L ttH-multilepton/functionsTTH.cc --tree treeProducerSusyMultilepton --mcc ttH-multilepton/mcc-lepeff-closure.txt --mcc ttH-multilepton/lepchoice-ttH-FO.txt --split-factor=-1 --WA prescaleFromSkim --perBin "# --neg"
    #CORE+=' '.join(["--plotgroup data_fakes%s+='.*_promptsub%s'"%(x,x) for x in ['','_FRe_norm_Up','_FRe_norm_Dn','_FRe_pt_Up','_FRe_pt_Dn','_FRe_be_Up','_FRe_be_Dn','_FRm_norm_Up','_FRm_norm_Dn','_FRm_pt_Up','_FRm_pt_Dn','_FRm_be_Up','_FRm_be_Dn']])+" --neglist '.*_promptsub.*' "
    RATIO= " --maxRatioRange 0.0  1.99 --ratioYNDiv 505 "
    RATIO2=" --showRatio --attachRatioPanel --fixRatioRange "
    LEGEND=" --legendColumns 2 --legendWidth 0.25 "
    LEGEND2=" --legendFontSize 0.042 "
    SPAM=" --noCms --topSpamSize 1.1 --lspam '#scale[1.1]{#bf{CMS}} #scale[0.9]{#it{Preliminary}}' "

    GO="%s ttH-multilepton/mca-2lss-mcdata-ttbar.txt ttH-multilepton/cuts_closure.txt "%CORE
    GO="%s -W 'vtxWeight2017'  "%GO#*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pt[iLepFO_Recl[0]],LepGood_eta[iLepFO_Recl[0]],2)*leptonSF_ttH(LepGood_pdgId[iLepFO_Recl[1]],LepGood_pt[iLepFO_Recl[1]],LepGood_eta[iLepFO_Recl[1]],2)*triggerSF_ttH(LepGood_pdgId[iLepFO_Recl[0]],LepGood_pdgId[iLepFO_Recl[1]],2)*eventBTagSF'"%GO
    GO+=" ttH-multilepton/closure_plots.txt {pseudo} ".format(pseudo='--pseudoData all' if False else '')
    GO=GO.replace(LEGEND, " --legendColumns 3 --legendWidth 0.46 ")
    GO=GO.replace(RATIO,  " --maxRatioRange 0.6  1.99 --ratioYNDiv 210 ")

    return GO

def runIt(GO,name,plots=[],noplots=[]):
    print 'python mcPlots.py',"--pdir %s/%s"%(ODIR,name),GO,' '.join(['--sP %s'%p for p in plots]),' '.join(['--xP %s'%p for p in noplots]),' '.join(sys.argv[3:])
    

if __name__ == '__main__':

    torun = sys.argv[2]

    x = base()
    #if '_data' in torun: x = x.replace('mca-2lss-mc.txt','mca-2lss-mcdata.txt')

    for tag, probe in [ ('Mu','El'), ('El','Mu') ]:
        for ty in ['pass','total']:
            #for gr in ['group']:
            for gr in ['nom', 'group']:
                x = x + " --xP '.*_pdgID' --xP '.*_pdgId' "
                x = x + " --xP '.*_isTight' "
                
                if ty == 'pass' and not '_nosf' in torun:
                    y = x.replace("-W 'vtxWeight2017'", " -W 'vtxWeight2017*_get_looseToTight_leptonSF_ttH({prob}_pdgId, {prob}_pt, {prob}_eta, 2)' ".format(prob='Electron' if probe == 'El' else 'Muon'))
                else: 
                    y = x 
                runIt(y + '-E {tag}Tight -E {tag}Tag_pt2510 {passing} {gr}'.format(tag=tag,
                                                                                   passing='-E {probe}Tight'.format(probe=probe) if ty=='pass' else '',
                                                                                   gr = ' --plotgroup TT_fake+=.*_fake --plotgroup TT_prompt+=.*_prompt ' if gr == 'group' else ''),
                      '%s_Tag%sProbe%s_%s_%s'%(torun,tag,probe,ty,gr))
                runIt(y + '-E {tag}Tight -E {tag}Tag_pt2510 {passing} -I opposite-sign {gr}'.format(tag=tag,
                                                                                                    passing='-E {probe}Tight'.format(probe=probe) if ty=='pass' else '',
                                                                                                    gr = ' --plotgroup TT_fake+=.*_fake --plotgroup TT_prompt+=.*_prompt ' if gr == 'group' else ''),
                      '%s_Tag%sProbe%s_%s_%s_SS'%(torun,tag,probe,ty,gr))
            

                        
