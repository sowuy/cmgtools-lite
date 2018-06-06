#include "TFile.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#include <iostream>
#include <algorithm>
#include "TRandom3.h"
#include <numeric>

float ttH_MVAto1D_6_2lss_Marco (float kinMVA_2lss_ttbar, float kinMVA_2lss_ttV){

  return 2*((kinMVA_2lss_ttbar>=-0.2)+(kinMVA_2lss_ttbar>=0.3))+(kinMVA_2lss_ttV>=-0.1)+1;

}
float ttH_MVAto1D_3_3l_Marco (float kinMVA_3l_ttbar, float kinMVA_3l_ttV){

  if (kinMVA_3l_ttbar<0.3 && kinMVA_3l_ttV<-0.1) return 1;
  else if (kinMVA_3l_ttbar>=0.3 && kinMVA_3l_ttV>=-0.1) return 3;
  else return 2;

}

#include "binning_2d_thresholds.h"
float ttH_MVAto1D_7_2lss_Marco (float kinMVA_2lss_ttbar, float kinMVA_2lss_ttV){

//________________
//|   |   |   | 7 |
//|   |   | 4 |___|
//| 1 | 2 |___| 6 |
//|   |   |   |___|
//|   |   | 3 | 5 |
//|___|___|___|___|
//

  if (kinMVA_2lss_ttbar<cuts_2lss_ttbar0) return 1;
  else if (kinMVA_2lss_ttbar<cuts_2lss_ttbar1) return 2;
  else if (kinMVA_2lss_ttbar<cuts_2lss_ttbar2) return 3+(kinMVA_2lss_ttV>=cuts_2lss_ttV0);
  else return 5+(kinMVA_2lss_ttV>=cuts_2lss_ttV1)+(kinMVA_2lss_ttV>=cuts_2lss_ttV2);

}
float ttH_MVAto1D_5_3l_Marco (float kinMVA_3l_ttbar, float kinMVA_3l_ttV){

  int reg = 2*((kinMVA_3l_ttbar>=cuts_3l_ttbar1)+(kinMVA_3l_ttbar>=cuts_3l_ttbar2))+(kinMVA_3l_ttV>=cuts_3l_ttV1)+1;
  if (reg==2) reg=1;
  if (reg>2) reg = reg-1;
  return reg;

}


float newBinning(float x, float y){
  float r =  4*((y>-0.16)+(y>0.28))+(x>-0.22)+(x>0.09)+(x>0.42)+1;
  if (r==9) r-=4;
  if (r>9) r-=1;
  return r;
}

#include "GetBinning.C"


float ttH_MVAto1D_6_flex (float kinMVA_2lss_ttbar, float kinMVA_2lss_ttV, int pdg1, int pdg2, float ttVcut, float ttcut1, float ttcut2){

  return 2*((kinMVA_2lss_ttbar>=ttcut1)+(kinMVA_2lss_ttbar>=ttcut2)) + (kinMVA_2lss_ttV>=ttVcut)+1;

}

float returnInputX(float x, float y) {return x;}

int ttH_catIndex_2lss(int LepGood1_pdgId, int LepGood2_pdgId, int LepGood1_charge, int nBJetMedium25){

//2lss_ee_neg
//2lss_ee_pos
//2lss_em_bl_neg
//2lss_em_bl_pos
//2lss_em_bt_neg
//2lss_em_bt_pos
//2lss_mm_bl_neg
//2lss_mm_bl_pos
//2lss_mm_bt_neg
//2lss_mm_bt_pos
   
  if (abs(LepGood1_pdgId)==11 && abs(LepGood2_pdgId)==11 && LepGood1_charge<0) return 2-1;
  if (abs(LepGood1_pdgId)==11 && abs(LepGood2_pdgId)==11 && LepGood1_charge>0) return 3-1;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && LepGood1_charge<0 && nBJetMedium25 < 2) return 4-1;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && LepGood1_charge>0 && nBJetMedium25 < 2) return 5-1;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && LepGood1_charge<0 && nBJetMedium25 >= 2) return 6-1;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && LepGood1_charge>0 && nBJetMedium25 >= 2) return 7-1;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && LepGood1_charge<0 && nBJetMedium25 < 2) return 8-1;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && LepGood1_charge>0 && nBJetMedium25 < 2) return 9-1;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && LepGood1_charge<0 && nBJetMedium25 >= 2) return 10-1;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && LepGood1_charge>0 && nBJetMedium25 >= 2) return 11-1;

 return -1;

}

int ttH_catIndex_2lss_nosign(int LepGood1_pdgId, int LepGood2_pdgId, int nBJetMedium25){

  if (abs(LepGood1_pdgId)==11 && abs(LepGood2_pdgId)==11) return 1;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && nBJetMedium25 < 2) return 2;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && nBJetMedium25 >= 2) return 3;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && nBJetMedium25 < 2) return 4;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && nBJetMedium25 >= 2) return 5;

 return -1;

}

int ttH_catIndex_3l(int LepGood1_charge, int LepGood2_charge, int LepGood3_charge, int nBJetMedium25){

//3l_bl_neg
//3l_bl_pos
//3l_bt_neg
//3l_bt_pos

  if ((LepGood1_charge+LepGood2_charge+LepGood3_charge)<0 && nBJetMedium25 < 2) return 11;
  if ((LepGood1_charge+LepGood2_charge+LepGood3_charge)>0 && nBJetMedium25 < 2) return 12;
  if ((LepGood1_charge+LepGood2_charge+LepGood3_charge)<0 && nBJetMedium25 >= 2) return 13;
  if ((LepGood1_charge+LepGood2_charge+LepGood3_charge)>0 && nBJetMedium25 >= 2) return 14;

 return -1;

}

TFile *_file_recoToLoose_leptonSF_mu1 = NULL;
TFile *_file_recoToLoose_leptonSF_mu2 = NULL;
TFile *_file_recoToLoose_leptonSF_mu3 = NULL;
TFile *_file_recoToLoose_leptonSF_mu4 = NULL;
TH2F *_histo_recoToLoose_leptonSF_mu1 = NULL;
TH2F *_histo_recoToLoose_leptonSF_mu2 = NULL;
TH2F *_histo_recoToLoose_leptonSF_mu3 = NULL;
TGraphAsymmErrors *_histo_recoToLoose_leptonSF_mu4 = NULL;
TFile *_file_recoToLoose_leptonSF_el = NULL;
TH2F *_histo_recoToLoose_leptonSF_el1 = NULL;
TH2F *_histo_recoToLoose_leptonSF_el2 = NULL;
TH2F *_histo_recoToLoose_leptonSF_el3 = NULL;
TFile *_file_recoToLoose_leptonSF_gsf = NULL;
TH2F *_histo_recoToLoose_leptonSF_gsf = NULL;

float _get_recoToLoose_leptonSF_ttH(int pdgid, float pt, float eta, int nlep, float var){

  // nlep is ignored for the loose selection

  if (!_histo_recoToLoose_leptonSF_mu1) {
    _file_recoToLoose_leptonSF_mu1 = new TFile("../../data/leptonSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root","read");
    _file_recoToLoose_leptonSF_mu2 = new TFile("../../data/leptonSF/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root","read");
    _file_recoToLoose_leptonSF_mu3 = new TFile("../../data/leptonSF/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","read");
    _file_recoToLoose_leptonSF_mu4 = new TFile("../../data/leptonSF/Tracking_EfficienciesAndSF_BCDEFGH.root","read");
    _histo_recoToLoose_leptonSF_mu1 = (TH2F*)(_file_recoToLoose_leptonSF_mu1->Get("SF"));
    _histo_recoToLoose_leptonSF_mu2 = (TH2F*)(_file_recoToLoose_leptonSF_mu2->Get("SF"));
    _histo_recoToLoose_leptonSF_mu3 = (TH2F*)(_file_recoToLoose_leptonSF_mu3->Get("SF"));
    _histo_recoToLoose_leptonSF_mu4 = (TGraphAsymmErrors*)(_file_recoToLoose_leptonSF_mu4->Get("ratio_eff_eta3_dr030e030_corr"));
  }
  if (!_histo_recoToLoose_leptonSF_el1) {
    _file_recoToLoose_leptonSF_el = new TFile("../../data/leptonSF/el_scaleFactors_Moriond17.root","read");
    _histo_recoToLoose_leptonSF_el1 = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("GsfElectronToMVAVLooseFOIDEmuTightIP2D"));
    _histo_recoToLoose_leptonSF_el2 = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToMini4"));
    _histo_recoToLoose_leptonSF_el3 = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToConvVetoIHit1"));
  }
  if (!_histo_recoToLoose_leptonSF_gsf) {
    _file_recoToLoose_leptonSF_gsf = new TFile("../../data/leptonSF/egammaEffi.txt_EGM2D.root","read");
    _histo_recoToLoose_leptonSF_gsf = (TH2F*)(_file_recoToLoose_leptonSF_gsf->Get("EGamma_SF2D"));
  }

  if (abs(pdgid)==13){

    // var is ignored for muons (handled in systsEnv.txt)

    float out = 1;

    TH2F *hist = _histo_recoToLoose_leptonSF_mu1;
    int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
    out *= hist->GetBinContent(ptbin,etabin);

    hist = _histo_recoToLoose_leptonSF_mu2;
    ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
    out *= hist->GetBinContent(ptbin,etabin);

    hist = _histo_recoToLoose_leptonSF_mu3;
    ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
    out *= hist->GetBinContent(ptbin,etabin);

    TGraphAsymmErrors *hist1 = _histo_recoToLoose_leptonSF_mu4;
    float eta1 = std::max(float(hist1->GetXaxis()->GetXmin()+1e-5), std::min(float(hist1->GetXaxis()->GetXmax()-1e-5), eta));
    out *= hist1->Eval(eta1);

    return out;

  }

  if (abs(pdgid)==11){
    TH2F *hist = _histo_recoToLoose_leptonSF_el1;
    int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(eta)));
    float out = hist->GetBinContent(ptbin,etabin)+var*hist->GetBinError(ptbin,etabin);
    hist = _histo_recoToLoose_leptonSF_el2;
    ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(eta)));
    out *= hist->GetBinContent(ptbin,etabin)+var*hist->GetBinError(ptbin,etabin);
    hist = _histo_recoToLoose_leptonSF_el3;
    ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(eta)));
    out *= hist->GetBinContent(ptbin,etabin)+var*hist->GetBinError(ptbin,etabin);

    hist = _histo_recoToLoose_leptonSF_gsf;
    etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(eta))); // careful, different convention
    ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pt)));
    out *= (hist->GetBinContent(etabin,ptbin)+var*(hist->GetBinError(ptbin,etabin) + 0.01*((pt<20) || (pt>80))));

    return out;
  }
  std::cout << "ERROR" << std::endl;
  std::abort();
  return -999;

}


std::vector<int> boundaries_runPeriod2017 = {297020,299337,302030,303435,304911};
std::vector<double> lumis_runPeriod2017 = {4.802,9.629,4.235,9.268,13.433};
bool cumul_lumis_runPeriod2017_isInit = false;
std::vector<float> cumul_lumis_runPeriod2017;

int runPeriod2017(int run){
  auto period = std::find_if(boundaries_runPeriod2017.begin(),boundaries_runPeriod2017.end(),[run](const int &y){return y>run;});
  return std::distance(boundaries_runPeriod2017.begin(),period)-1;
}



TRandom3 rand_generator_RunDependentMC(0);
int hashBasedRunDependentMC(int isData, int run, int lumi, int event){
  if (isData) return runPeriod2017(run);
  if (!cumul_lumis_runPeriod2017_isInit){
    cumul_lumis_runPeriod2017.push_back(0);
    float tot_lumi = std::accumulate(lumis_runPeriod2017.begin(),lumis_runPeriod2017.end(),float(0.0));
    for (uint i=0; i<lumis_runPeriod2017.size(); i++) cumul_lumis_runPeriod2017.push_back(cumul_lumis_runPeriod2017.back()+lumis_runPeriod2017[i]/tot_lumi);
    cumul_lumis_runPeriod2017_isInit = true;
  }
  Int_t x = 161248*run+2136324*lumi+12781432*event;
  unsigned int hash = TString::Hash(&x,sizeof(Int_t));
  rand_generator_RunDependentMC.SetSeed(hash);
  float val = rand_generator_RunDependentMC.Uniform();
  auto period = std::find_if(cumul_lumis_runPeriod2017.begin(),cumul_lumis_runPeriod2017.end(),[val](const float &y){return y>val;});
  return std::distance(cumul_lumis_runPeriod2017.begin(),period)-1;
}












TFile *_file_looseToTight_leptonSF_B_mu_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_B_mu_2lss = NULL;
TFile *_file_looseToTight_leptonSF_B_el_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_B_el_2lss = NULL;
TFile *_file_looseToTight_leptonSF_B_mu_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_B_mu_3l = NULL;
TFile *_file_looseToTight_leptonSF_B_el_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_B_el_3l = NULL;


TFile *_file_looseToTight_leptonSF_CD_mu_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_CD_mu_2lss = NULL;
TFile *_file_looseToTight_leptonSF_CD_el_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_CD_el_2lss = NULL;
TFile *_file_looseToTight_leptonSF_CD_mu_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_CD_mu_3l = NULL;
TFile *_file_looseToTight_leptonSF_CD_el_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_CD_el_3l = NULL;

TFile *_file_looseToTight_leptonSF_E_mu_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_E_mu_2lss = NULL;
TFile *_file_looseToTight_leptonSF_E_el_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_E_el_2lss = NULL;
TFile *_file_looseToTight_leptonSF_E_mu_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_E_mu_3l = NULL;
TFile *_file_looseToTight_leptonSF_E_el_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_E_el_3l = NULL;

TFile *_file_looseToTight_leptonSF_F_mu_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_F_mu_2lss = NULL;
TFile *_file_looseToTight_leptonSF_F_el_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_F_el_2lss = NULL;
TFile *_file_looseToTight_leptonSF_F_mu_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_F_mu_3l = NULL;
TFile *_file_looseToTight_leptonSF_F_el_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_F_el_3l = NULL;



float _get_looseToTight_leptonSF_ttH(int pdgid, float pt, float eta, int nlep, Int_t run, Int_t lumi, Int_t evt, float var=0){

  int era = hashBasedRunDependentMC(false, run, lumi, evt);

  // /nfs/fanae/user/sscruz/www/ttH/may31/SFs
  
  if (nlep == 2){
    if (!_histo_looseToTight_leptonSF_B_mu_2lss) {
      std::cout << "Openning files" << std::endl;
      _file_looseToTight_leptonSF_B_mu_2lss = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runB/lepMVAEffSF_m_2lss.root","read");
      _file_looseToTight_leptonSF_CD_mu_2lss = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runCD/lepMVAEffSF_m_2lss.root","read");
      _file_looseToTight_leptonSF_E_mu_2lss = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runE/lepMVAEffSF_m_2lss.root","read");
      _file_looseToTight_leptonSF_F_mu_2lss = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runF/lepMVAEffSF_m_2lss.root","read");

      _histo_looseToTight_leptonSF_B_mu_2lss = (TH2F*)(_file_looseToTight_leptonSF_B_mu_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_CD_mu_2lss = (TH2F*)(_file_looseToTight_leptonSF_CD_mu_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_E_mu_2lss = (TH2F*)(_file_looseToTight_leptonSF_E_mu_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_F_mu_2lss = (TH2F*)(_file_looseToTight_leptonSF_F_mu_2lss->Get("sf"));
    }
    if (!_histo_looseToTight_leptonSF_B_el_2lss) {
      std::cout << "Openning files2" << std::endl;
      _file_looseToTight_leptonSF_B_el_2lss  = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runB/lepMVAEffSF_e_2lss.root","read");
      _file_looseToTight_leptonSF_CD_el_2lss = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runCD/lepMVAEffSF_e_2lss.root","read");
      _file_looseToTight_leptonSF_E_el_2lss  = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runE/lepMVAEffSF_e_2lss.root","read");
      _file_looseToTight_leptonSF_F_el_2lss  = new TFile("/nfs/fanae/user/sscruz/www/ttH/may31/SFs/runF/lepMVAEffSF_e_2lss.root","read");

      _histo_looseToTight_leptonSF_B_el_2lss = (TH2F*)(_file_looseToTight_leptonSF_B_el_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_CD_el_2lss = (TH2F*)(_file_looseToTight_leptonSF_CD_el_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_E_el_2lss = (TH2F*)(_file_looseToTight_leptonSF_E_el_2lss->Get("sf"));
      _histo_looseToTight_leptonSF_F_el_2lss = (TH2F*)(_file_looseToTight_leptonSF_F_el_2lss->Get("sf"));


    }
  }
  else{
    std::cout << "3l selection is missing" << std::endl;
    assert(0);
  }

  TH2F *hist = 0;
  if (era == 0){
    if (abs(pdgid)==13) hist = (nlep>2) ? _histo_looseToTight_leptonSF_B_mu_3l : _histo_looseToTight_leptonSF_B_mu_2lss;
    else if (abs(pdgid)==11) hist = (nlep>2) ? _histo_looseToTight_leptonSF_B_el_3l : _histo_looseToTight_leptonSF_B_el_2lss;
  }
  else if (era ==1 || era ==2){
    if (abs(pdgid)==13) hist = (nlep>2) ? _histo_looseToTight_leptonSF_CD_mu_3l : _histo_looseToTight_leptonSF_CD_mu_2lss;
    else if (abs(pdgid)==11) hist = (nlep>2) ? _histo_looseToTight_leptonSF_CD_el_3l : _histo_looseToTight_leptonSF_CD_el_2lss;

  }

  else if (era == 3){
    if (abs(pdgid)==13) hist = (nlep>2) ? _histo_looseToTight_leptonSF_E_mu_3l : _histo_looseToTight_leptonSF_E_mu_2lss;
    else if (abs(pdgid)==11) hist = (nlep>2) ? _histo_looseToTight_leptonSF_E_el_3l : _histo_looseToTight_leptonSF_E_el_2lss;
  }
  else if (era ==4){
    if (abs(pdgid)==13) hist = (nlep>2) ? _histo_looseToTight_leptonSF_F_mu_3l : _histo_looseToTight_leptonSF_F_mu_2lss;
    else if (abs(pdgid)==11) hist = (nlep>2) ? _histo_looseToTight_leptonSF_F_el_3l : _histo_looseToTight_leptonSF_F_el_2lss;
  }
  else{
    std::cout << "Some run era is missing" << std::endl;
    assert(0);
  }

  // if (!hist) {std::cout << "ERROR" << std::endl; std::abort();}
  if (!hist){ return -999999; } // it shoudlnt matter for non emu events
  int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
  int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));




  return hist->GetBinContent(ptbin,etabin);

}

float leptonSF_ttH(int pdgid, float pt, float eta, int nlep, Int_t run, Int_t lumi, Int_t evt, float var=0){

  float recoToLoose = _get_recoToLoose_leptonSF_ttH(pdgid,pt,eta,nlep,var);
  float looseToTight = _get_looseToTight_leptonSF_ttH(pdgid,pt,eta,nlep,run, lumi, evt,var);
  float res = recoToLoose*looseToTight;
  if (!(res>0)) {std::cout << "ERROR" << std::endl; std::abort();}
  return res;

}

//TFile *file_triggerSF_ttH = NULL;
//TH2Poly* t2poly_triggerSF_ttH_mm = NULL;
//TH2Poly* t2poly_triggerSF_ttH_ee = NULL;
//TH2Poly* t2poly_triggerSF_ttH_em = NULL;
//TH2Poly* t2poly_triggerSF_ttH_3l = NULL;

float triggerSF_ttH(int pdgid1, int pdgid2, int nlep){
  
  if (nlep>=3) return 1;
  
  int comb = abs(pdgid1)+abs(pdgid2);
  if (comb==22) return 1.01; // ee
  else if (comb==24) return 1.01; // em
  else if (comb==26) return 1; // mm

  return 1;

}







//float triggerSF_ttH(int pdgid1, float pt1, int pdgid2, float pt2, int nlep, float var=0){
//
//  if (!file_triggerSF_ttH) {
//    file_triggerSF_ttH = new TFile("../../data/triggerSF/trig_eff_map_v4.root");
//    t2poly_triggerSF_ttH_mm = (TH2Poly*)(file_triggerSF_ttH->Get("SSuu2DPt_effic"));
//    t2poly_triggerSF_ttH_ee = (TH2Poly*)(file_triggerSF_ttH->Get("SSee2DPt__effic"));
//    t2poly_triggerSF_ttH_em = (TH2Poly*)(file_triggerSF_ttH->Get("SSeu2DPt_effic"));
//    t2poly_triggerSF_ttH_3l = (TH2Poly*)(file_triggerSF_ttH->Get("__3l2DPt_effic"));
//    if (!(t2poly_triggerSF_ttH_mm && t2poly_triggerSF_ttH_ee && t2poly_triggerSF_ttH_em && t2poly_triggerSF_ttH_3l)) {
//	std::cout << "Impossible to load trigger scale factors!" << std::endl;
//	file_triggerSF_ttH->ls();
//	file_triggerSF_ttH = NULL;
//      }
//  }
//  TH2Poly* hist = NULL;
//  if (nlep==2){
//    if (abs(pdgid1)==13 && abs(pdgid2)==13) hist = t2poly_triggerSF_ttH_mm;
//    else if (abs(pdgid1)==11 && abs(pdgid2)==11) hist = t2poly_triggerSF_ttH_ee;
//    else hist = t2poly_triggerSF_ttH_em;
//  }
//  else if (nlep==3) hist = t2poly_triggerSF_ttH_3l;
//  else std::cout << "Wrong options to trigger scale factors" << std::endl;
//  pt1 = std::max(float(hist->GetXaxis()->GetXmin()+1e-5), std::min(float(hist->GetXaxis()->GetXmax()-1e-5), pt1));
//  pt2 = std::max(float(hist->GetYaxis()->GetXmin()+1e-5), std::min(float(hist->GetYaxis()->GetXmax()-1e-5), pt2));
//  int bin = hist->FindBin(pt1,pt2);
//  float eff = hist->GetBinContent(bin) + var * hist->GetBinError(bin);
//
//  if (nlep>2) return eff;
//  int cat = (abs(pdgid1)==11) + (abs(pdgid2)==11);
//  if (cat==2) return eff*1.02;
//  else if (cat==1) return eff*1.02;
//  else return eff*1.01;
//
//
//}


float ttH_2lss_ifflav(int LepGood1_pdgId, int LepGood2_pdgId, float ret_ee, float ret_em, float ret_mm){
  if (abs(LepGood1_pdgId)==11 && abs(LepGood2_pdgId)==11) return ret_ee;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)))       return ret_em;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13) return ret_mm;
  std::cerr << "ERROR: invalid input " << abs(LepGood1_pdgId) << ", " << abs(LepGood1_pdgId) << std::endl;
  assert(0);
}
float ttH_2lss_ifflavnb(int LepGood1_pdgId, int LepGood2_pdgId, int nBJetMedium25, float ret_ee, float ret_em_bl, float ret_em_bt, float ret_mm_bl, float ret_mm_bt){
  if (abs(LepGood1_pdgId)==11 && abs(LepGood2_pdgId)==11) return ret_ee;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && nBJetMedium25 < 2) return ret_em_bl;
  if ((abs(LepGood1_pdgId) != abs(LepGood2_pdgId)) && nBJetMedium25 >= 2) return ret_em_bt;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && nBJetMedium25 < 2) return ret_mm_bl;
  if (abs(LepGood1_pdgId)==13 && abs(LepGood2_pdgId)==13 && nBJetMedium25 >= 2) return ret_mm_bt;
  std::cerr << "ERROR: invalid input " << abs(LepGood1_pdgId) << ", " << abs(LepGood1_pdgId) <<  ", " << nBJetMedium25 << std::endl;
  assert(0);
}


TH1D* _histo_eff_el_2lss_pass  = 0;
TH1D* _histo_eff_el_2lss_total = 0;
TH1D* _histo_eff_mu_2lss_pass  = 0;
TH1D* _histo_eff_mu_2lss_total = 0;
TEfficiency* eff_el_2lss = 0;
TEfficiency* eff_mu_2lss = 0;

float getEff( float LepGood1_pt, int LepGood1_pdgId, int nlep  ){

  if (!eff_el_2lss) {
    TFile* _file_eff_el_2lss_pass   = new TFile("/nfs/fanae/user/sscruz/TTH/may03/CMSSW_9_4_4/src/CMGTools/TTHAnalysis/python/plotter/closure_withData/closure_withData_TagMuProbeEl_pass_group/closure_plots.root","read");
    _histo_eff_el_2lss_pass  = (TH1D*) ((TH1D*)_file_eff_el_2lss_pass->Get("Electron_pt_TT_prompt"))->Clone("Electron_pt_TT_prompt_pass");
    TFile* _file_eff_el_2lss_total  = new TFile("/nfs/fanae/user/sscruz/TTH/may03/CMSSW_9_4_4/src/CMGTools/TTHAnalysis/python/plotter/closure_withData/closure_withData_TagMuProbeEl_total_group/closure_plots.root","read");
    _histo_eff_el_2lss_total = (TH1D*) ((TH1D*) _file_eff_el_2lss_total->Get("Electron_pt_TT_prompt"))->Clone("Electron_pt_TT_prompt_total");
    _histo_eff_el_2lss_pass->Draw();
    eff_el_2lss = new TEfficiency(*_histo_eff_el_2lss_pass, *_histo_eff_el_2lss_total);
  }

  if (!eff_mu_2lss) {
    TFile* _file_eff_mu_2lss_pass   = new TFile("/nfs/fanae/user/sscruz/TTH/may03/CMSSW_9_4_4/src/CMGTools/TTHAnalysis/python/plotter/closure_withData/closure_withData_TagElProbeMu_pass_group/closure_plots.root","read");
    _histo_eff_mu_2lss_pass  =  (TH1D*) ((TH1D*) _file_eff_mu_2lss_pass->Get("Muon_pt_TT_prompt"))->Clone("Muon_pt_TT_prompt_pass");
    TFile* _file_eff_mu_2lss_total  = new TFile("/nfs/fanae/user/sscruz/TTH/may03/CMSSW_9_4_4/src/CMGTools/TTHAnalysis/python/plotter/closure_withData/closure_withData_TagElProbeMu_total_group/closure_plots.root","read");
    _histo_eff_mu_2lss_total =  (TH1D*) ((TH1D*) _file_eff_mu_2lss_total->Get("Muon_pt_TT_prompt"))->Clone("Muon_pt_TT_prompt_total");

    eff_mu_2lss = new TEfficiency(*_histo_eff_mu_2lss_pass, *_histo_eff_mu_2lss_total);
  }
  
  if (nlep != 2) {
    std::cerr << "ERROR: invalid input. nlep != 2 not supported yet" << std::endl;
    assert(0);
  }

  TEfficiency *eff = 0;
  if (abs(LepGood1_pdgId)==13) eff = eff_mu_2lss;
  else                eff = eff_el_2lss;

  return eff->GetEfficiency( eff->FindFixBin( LepGood1_pt ) );

}



float getAntiSF(int pdgid,float pt,float eta, int nlep, Int_t run, Int_t lumi, Int_t evt,int var=0){
  float SF  = _get_looseToTight_leptonSF_ttH(pdgid,pt,eta,nlep,run, lumi, evt,var);
  float eff =  getEff( pt, pdgid, nlep  );
  //std::cout << (1-SF*eff)/(1-eff) << " " << eff << std::endl;
  return (1-SF*eff)/(1-eff);

}


