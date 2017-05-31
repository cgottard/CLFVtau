//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 31 15:43:01 2017 by ROOT version 6.04/14
// from TTree nominal/tree
// found on file: Signal_v9ta.root
//////////////////////////////////////////////////////////

#ifndef nominal_h
#define nominal_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>

// Header file for the classes stored in the TTree if any.



class nominal {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<float>   *PDFinfo_X1;
   std::vector<float>   *PDFinfo_X2;
   std::vector<int>     *PDFinfo_PDGID1;
   std::vector<int>     *PDFinfo_PDGID2;
   std::vector<float>   *PDFinfo_Q;
   std::vector<float>   *PDFinfo_XF1;
   std::vector<float>   *PDFinfo_XF2;
   Float_t         weight_mc;
   Float_t         weight_pileup;
   Float_t         weight_leptonSF;
   Float_t         weight_tauSF;
   Float_t         weight_bTagSF_77;
   Float_t         weight_jvt;
   Float_t         weight_pileup_UP;
   Float_t         weight_pileup_DOWN;
   Float_t         weight_leptonSF_EL_SF_Trigger_UP;
   Float_t         weight_leptonSF_EL_SF_Trigger_DOWN;
   Float_t         weight_leptonSF_EL_SF_Reco_UP;
   Float_t         weight_leptonSF_EL_SF_Reco_DOWN;
   Float_t         weight_leptonSF_EL_SF_ID_UP;
   Float_t         weight_leptonSF_EL_SF_ID_DOWN;
   Float_t         weight_leptonSF_EL_SF_Isol_UP;
   Float_t         weight_leptonSF_EL_SF_Isol_DOWN;
   Float_t         weight_leptonSF_MU_SF_Trigger_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_Trigger_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Trigger_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_Trigger_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_DOWN;
   Float_t         weight_indiv_SF_EL_Trigger;
   Float_t         weight_indiv_SF_EL_Trigger_UP;
   Float_t         weight_indiv_SF_EL_Trigger_DOWN;
   Float_t         weight_indiv_SF_EL_Reco;
   Float_t         weight_indiv_SF_EL_Reco_UP;
   Float_t         weight_indiv_SF_EL_Reco_DOWN;
   Float_t         weight_indiv_SF_EL_ID;
   Float_t         weight_indiv_SF_EL_ID_UP;
   Float_t         weight_indiv_SF_EL_ID_DOWN;
   Float_t         weight_indiv_SF_EL_Isol;
   Float_t         weight_indiv_SF_EL_Isol_UP;
   Float_t         weight_indiv_SF_EL_Isol_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeID;
   Float_t         weight_indiv_SF_EL_ChargeID_UP;
   Float_t         weight_indiv_SF_EL_ChargeID_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeMisID;
   Float_t         weight_indiv_SF_EL_ChargeMisID_STAT_UP;
   Float_t         weight_indiv_SF_EL_ChargeMisID_STAT_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeMisID_SYST_UP;
   Float_t         weight_indiv_SF_EL_ChargeMisID_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_Trigger;
   Float_t         weight_indiv_SF_MU_Trigger_STAT_UP;
   Float_t         weight_indiv_SF_MU_Trigger_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_Trigger_SYST_UP;
   Float_t         weight_indiv_SF_MU_Trigger_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_ID;
   Float_t         weight_indiv_SF_MU_ID_STAT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_Isol;
   Float_t         weight_indiv_SF_MU_Isol_STAT_UP;
   Float_t         weight_indiv_SF_MU_Isol_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_Isol_SYST_UP;
   Float_t         weight_indiv_SF_MU_Isol_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_TTVA;
   Float_t         weight_indiv_SF_MU_TTVA_STAT_UP;
   Float_t         weight_indiv_SF_MU_TTVA_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_TTVA_SYST_UP;
   Float_t         weight_indiv_SF_MU_TTVA_SYST_DOWN;
   Float_t         weight_tauSF_ELEOLR_UP;
   Float_t         weight_tauSF_ELEOLR_DOWN;
   Float_t         weight_tauSF_JETID_UP;
   Float_t         weight_tauSF_JETID_DOWN;
   Float_t         weight_tauSF_RECO_UP;
   Float_t         weight_tauSF_RECO_DOWN;
   Float_t         weight_jvt_UP;
   Float_t         weight_jvt_DOWN;
   std::vector<float>   *weight_bTagSF_77_eigenvars_B_up;
   std::vector<float>   *weight_bTagSF_77_eigenvars_C_up;
   std::vector<float>   *weight_bTagSF_77_eigenvars_Light_up;
   std::vector<float>   *weight_bTagSF_77_eigenvars_B_down;
   std::vector<float>   *weight_bTagSF_77_eigenvars_C_down;
   std::vector<float>   *weight_bTagSF_77_eigenvars_Light_down;
   Float_t         weight_bTagSF_77_extrapolation_up;
   Float_t         weight_bTagSF_77_extrapolation_down;
   Float_t         weight_bTagSF_77_extrapolation_from_charm_up;
   Float_t         weight_bTagSF_77_extrapolation_from_charm_down;
   ULong64_t       eventNumber;
   UInt_t          runNumber;
   UInt_t          randomRunNumber;
   UInt_t          mcChannelNumber;
   Float_t         mu;
   UInt_t          backgroundFlags;
   UInt_t          hasBadMuon;
   std::vector<float>   *el_pt;
   std::vector<float>   *el_eta;
   std::vector<float>   *el_cl_eta;
   std::vector<float>   *el_phi;
   std::vector<float>   *el_e;
   std::vector<float>   *el_charge;
   std::vector<float>   *el_topoetcone20;
   std::vector<float>   *el_ptvarcone20;
   std::vector<char>    *el_CF;
   std::vector<float>   *el_d0sig;
   std::vector<float>   *el_delta_z0_sintheta;
   std::vector<int>     *el_true_type;
   std::vector<int>     *el_true_origin;
   std::vector<int>     *el_true_typebkg;
   std::vector<int>     *el_true_originbkg;
   std::vector<float>   *mu_pt;
   std::vector<float>   *mu_eta;
   std::vector<float>   *mu_phi;
   std::vector<float>   *mu_e;
   std::vector<float>   *mu_charge;
   std::vector<float>   *mu_topoetcone20;
   std::vector<float>   *mu_ptvarcone30;
   std::vector<float>   *mu_d0sig;
   std::vector<float>   *mu_delta_z0_sintheta;
   std::vector<int>     *mu_true_type;
   std::vector<int>     *mu_true_origin;
   std::vector<float>   *tau_pt;
   std::vector<float>   *tau_eta;
   std::vector<float>   *tau_phi;
   std::vector<float>   *tau_charge;
   std::vector<float>   *jet_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_e;
   std::vector<float>   *jet_mv2c00;
   std::vector<float>   *jet_mv2c10;
   std::vector<float>   *jet_mv2c20;
   std::vector<float>   *jet_ip3dsv1;
   std::vector<float>   *jet_jvt;
   std::vector<char>    *jet_passfjvt;
   std::vector<int>     *jet_truthflav;
   std::vector<int>     *jet_truthPartonLabel;
   std::vector<char>    *jet_isTrueHS;
   std::vector<char>    *jet_isbtagged_77;
   Float_t         met_met;
   Float_t         met_phi;
   Int_t           Selection2015;
   Int_t           Selection2016;
   Char_t          HLT_e60_lhmedium_nod0;
   Char_t          HLT_mu26_ivarmedium;
   Char_t          HLT_e26_lhtight_nod0_ivarloose;
   Char_t          HLT_e140_lhloose_nod0;
   Char_t          HLT_mu20_iloose_L1MU15;
   Char_t          HLT_mu50;
   Char_t          HLT_e60_lhmedium;
   Char_t          HLT_e24_lhmedium_L1EM20VH;
   Char_t          HLT_e120_lhloose;
   std::vector<char>    *el_trigMatch_HLT_e60_lhmedium_nod0;
   std::vector<char>    *el_trigMatch_HLT_e26_lhtight_nod0_ivarloose;
   std::vector<char>    *el_trigMatch_HLT_e140_lhloose_nod0;
   std::vector<char>    *el_trigMatch_HLT_e60_lhmedium;
   std::vector<char>    *el_trigMatch_HLT_e24_lhmedium_L1EM20VH;
   std::vector<char>    *el_trigMatch_HLT_e120_lhloose;
   std::vector<char>    *mu_trigMatch_HLT_mu26_ivarmedium;
   std::vector<char>    *mu_trigMatch_HLT_mu50;
   std::vector<char>    *mu_trigMatch_HLT_mu20_iloose_L1MU15;
   Float_t         nEl;
   Float_t         nMu;
   Float_t         nJets;
   Float_t         nBJets70;
   Float_t         nBJets77;
   Float_t         nBJets85;
   Float_t         pT_1lep;
   Float_t         pT_2lep;
   Float_t         pT_3lep;
   Float_t         eta_1lep;
   Float_t         eta_2lep;
   Float_t         eta_3lep;
   Float_t         phi_1lep;
   Float_t         phi_2lep;
   Float_t         phi_3lep;
   Float_t         pdgid_1lep;
   Float_t         pdgid_2lep;
   Float_t         pdgid_3lep;
   Float_t         charge_1lep;
   Float_t         charge_2lep;
   Float_t         charge_3lep;
   Float_t         pT_1jet;
   Float_t         pT_2jet;
   Float_t         pT_3jet;
   Float_t         ht;
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           lep1_origin;
   Int_t           lep2_origin;
   Int_t           lep3_origin;
   Int_t           lep1_pdgorigin;
   Int_t           lep2_pdgorigin;
   Int_t           lep3_pdgorigin;
   Int_t           lep1_istight;
   Int_t           lep2_istight;
   Int_t           lep3_istight;
   Int_t           lep1_type;
   Int_t           lep2_type;
   Int_t           lep3_type;
   Float_t         lep1_d0sig;
   Float_t         lep2_d0sig;
   Float_t         lep3_d0sig;
   Float_t         lep1_z0sintheta;
   Float_t         lep2_z0sintheta;
   Float_t         lep3_z0sintheta;
   Float_t         lep1_topoetcone;
   Float_t         lep2_topoetcone;
   Float_t         lep3_topoetcone;
   Float_t         lep1_topoetcone40;
   Float_t         lep2_topoetcone40;
   Float_t         lep3_topoetcone40;
   Float_t         lep1_ptvarcone;
   Float_t         lep2_ptvarcone;
   Float_t         lep3_ptvarcone;
   Float_t         weight0tag;
   Float_t         weight1tag;
   Float_t         weight2tag;
   Bool_t          lep1_Isol_LooseTrackOnly;
   Bool_t          lep1_Isol_Loose;
   Bool_t          lep1_Isol_Gradient;
   Bool_t          lep1_Isol_GradientLoose;
   Bool_t          lep1_Isol_FixedCutTightTrackOnly;
   Bool_t          lep1_Isol_FixedCutLoose;
   Bool_t          lep1_Isol_FixedCutTight;
   Bool_t          lep2_Isol_LooseTrackOnly;
   Bool_t          lep2_Isol_Loose;
   Bool_t          lep2_Isol_Gradient;
   Bool_t          lep2_Isol_GradientLoose;
   Bool_t          lep2_Isol_FixedCutTightTrackOnly;
   Bool_t          lep2_Isol_FixedCutLoose;
   Bool_t          lep2_Isol_FixedCutTight;
   Bool_t          lep3_Isol_LooseTrackOnly;
   Bool_t          lep3_Isol_Loose;
   Bool_t          lep3_Isol_Gradient;
   Bool_t          lep3_Isol_GradientLoose;
   Bool_t          lep3_Isol_FixedCutTightTrackOnly;
   Bool_t          lep3_Isol_FixedCutLoose;
   Bool_t          lep3_Isol_FixedCutTight;
   Bool_t          hasTrigMatchedEl;
   Bool_t          hasTrigMatchedMu;
   Int_t           TopHeavyFlavorFilterFlag;
   Bool_t          MEgamma;
   Int_t           numberTruePhotons;
   std::vector<float>   *MTW;
   std::vector<float>   *el_ECID;
   std::vector<float>   *OSSF_M;
   Int_t           nOSSF;
   Int_t           nOS;
   Float_t         mll;
   Bool_t          prompt;
   Float_t         xs_weight;

   // List of branches
   TBranch        *b_PDFinfo_X1;   //!
   TBranch        *b_PDFinfo_X2;   //!
   TBranch        *b_PDFinfo_PDGID1;   //!
   TBranch        *b_PDFinfo_PDGID2;   //!
   TBranch        *b_PDFinfo_Q;   //!
   TBranch        *b_PDFinfo_XF1;   //!
   TBranch        *b_PDFinfo_XF2;   //!
   TBranch        *b_weight_mc;   //!
   TBranch        *b_weight_pileup;   //!
   TBranch        *b_weight_leptonSF;   //!
   TBranch        *b_weight_tauSF;   //!
   TBranch        *b_weight_bTagSF_77;   //!
   TBranch        *b_weight_jvt;   //!
   TBranch        *b_weight_pileup_UP;   //!
   TBranch        *b_weight_pileup_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Trigger_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Trigger_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Reco_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Reco_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_ID_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_ID_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Isol_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Isol_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ID;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Trigger;   //!
   TBranch        *b_weight_indiv_SF_MU_Trigger_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Trigger_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Trigger_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Trigger_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_LOWPT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_LOWPT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_SYST_DOWN;   //!
   TBranch        *b_weight_tauSF_ELEOLR_UP;   //!
   TBranch        *b_weight_tauSF_ELEOLR_DOWN;   //!
   TBranch        *b_weight_tauSF_JETID_UP;   //!
   TBranch        *b_weight_tauSF_JETID_DOWN;   //!
   TBranch        *b_weight_tauSF_RECO_UP;   //!
   TBranch        *b_weight_tauSF_RECO_DOWN;   //!
   TBranch        *b_weight_jvt_UP;   //!
   TBranch        *b_weight_jvt_DOWN;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_B_up;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_C_up;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_Light_up;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_B_down;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_C_down;   //!
   TBranch        *b_weight_bTagSF_77_eigenvars_Light_down;   //!
   TBranch        *b_weight_bTagSF_77_extrapolation_up;   //!
   TBranch        *b_weight_bTagSF_77_extrapolation_down;   //!
   TBranch        *b_weight_bTagSF_77_extrapolation_from_charm_up;   //!
   TBranch        *b_weight_bTagSF_77_extrapolation_from_charm_down;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_randomRunNumber;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_backgroundFlags;   //!
   TBranch        *b_hasBadMuon;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_cl_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_topoetcone20;   //!
   TBranch        *b_el_ptvarcone20;   //!
   TBranch        *b_el_CF;   //!
   TBranch        *b_el_d0sig;   //!
   TBranch        *b_el_delta_z0_sintheta;   //!
   TBranch        *b_el_true_type;   //!
   TBranch        *b_el_true_origin;   //!
   TBranch        *b_el_true_typebkg;   //!
   TBranch        *b_el_true_originbkg;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_topoetcone20;   //!
   TBranch        *b_mu_ptvarcone30;   //!
   TBranch        *b_mu_d0sig;   //!
   TBranch        *b_mu_delta_z0_sintheta;   //!
   TBranch        *b_mu_true_type;   //!
   TBranch        *b_mu_true_origin;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_mv2c00;   //!
   TBranch        *b_jet_mv2c10;   //!
   TBranch        *b_jet_mv2c20;   //!
   TBranch        *b_jet_ip3dsv1;   //!
   TBranch        *b_jet_jvt;   //!
   TBranch        *b_jet_passfjvt;   //!
   TBranch        *b_jet_truthflav;   //!
   TBranch        *b_jet_truthPartonLabel;   //!
   TBranch        *b_jet_isTrueHS;   //!
   TBranch        *b_jet_isbtagged_77;   //!
   TBranch        *b_met_met;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_Selection2015;   //!
   TBranch        *b_Selection2016;   //!
   TBranch        *b_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_HLT_mu26_ivarmedium;   //!
   TBranch        *b_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_HLT_mu50;   //!
   TBranch        *b_HLT_e60_lhmedium;   //!
   TBranch        *b_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_HLT_e120_lhloose;   //!
   TBranch        *b_el_trigMatch_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_el_trigMatch_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_el_trigMatch_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_el_trigMatch_HLT_e60_lhmedium;   //!
   TBranch        *b_el_trigMatch_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_el_trigMatch_HLT_e120_lhloose;   //!
   TBranch        *b_mu_trigMatch_HLT_mu26_ivarmedium;   //!
   TBranch        *b_mu_trigMatch_HLT_mu50;   //!
   TBranch        *b_mu_trigMatch_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_nEl;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_nBJets70;   //!
   TBranch        *b_nBJets77;   //!
   TBranch        *b_nBJets85;   //!
   TBranch        *b_pT_1lep;   //!
   TBranch        *b_pT_2lep;   //!
   TBranch        *b_pT_3lep;   //!
   TBranch        *b_eta_1lep;   //!
   TBranch        *b_eta_2lep;   //!
   TBranch        *b_eta_3lep;   //!
   TBranch        *b_phi_1lep;   //!
   TBranch        *b_phi_2lep;   //!
   TBranch        *b_phi_3lep;   //!
   TBranch        *b_pdgid_1lep;   //!
   TBranch        *b_pdgid_2lep;   //!
   TBranch        *b_pdgid_3lep;   //!
   TBranch        *b_charge_1lep;   //!
   TBranch        *b_charge_2lep;   //!
   TBranch        *b_charge_3lep;   //!
   TBranch        *b_pT_1jet;   //!
   TBranch        *b_pT_2jet;   //!
   TBranch        *b_pT_3jet;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lep1_origin;   //!
   TBranch        *b_lep2_origin;   //!
   TBranch        *b_lep3_origin;   //!
   TBranch        *b_lep1_pdgorigin;   //!
   TBranch        *b_lep2_pdgorigin;   //!
   TBranch        *b_lep3_pdgorigin;   //!
   TBranch        *b_lep1_istight;   //!
   TBranch        *b_lep2_istight;   //!
   TBranch        *b_lep3_istight;   //!
   TBranch        *b_lep1_type;   //!
   TBranch        *b_lep2_type;   //!
   TBranch        *b_lep3_type;   //!
   TBranch        *b_lep1_d0sig;   //!
   TBranch        *b_lep2_d0sig;   //!
   TBranch        *b_lep3_d0sig;   //!
   TBranch        *b_lep1_z0sintheta;   //!
   TBranch        *b_lep2_z0sintheta;   //!
   TBranch        *b_lep3_z0sintheta;   //!
   TBranch        *b_lep1_topoetcone;   //!
   TBranch        *b_lep2_topoetcone;   //!
   TBranch        *b_lep3_topoetcone;   //!
   TBranch        *b_lep1_topoetcone40;   //!
   TBranch        *b_lep2_topoetcone40;   //!
   TBranch        *b_lep3_topoetcone40;   //!
   TBranch        *b_lep1_ptvarcone;   //!
   TBranch        *b_lep2_ptvarcone;   //!
   TBranch        *b_lep3_ptvarcone;   //!
   TBranch        *b_weight0tag;   //!
   TBranch        *b_weight1tag;   //!
   TBranch        *b_weight2tag;   //!
   TBranch        *b_lep1_Isol_LooseTrackOnly;   //!
   TBranch        *b_lep1_Isol_Loose;   //!
   TBranch        *b_lep1_Isol_Gradient;   //!
   TBranch        *b_lep1_Isol_GradientLoose;   //!
   TBranch        *b_lep1_Isol_FixedCutTightTrackOnly;   //!
   TBranch        *b_lep1_Isol_FixedCutLoose;   //!
   TBranch        *b_lep1_Isol_FixedCutTight;   //!
   TBranch        *b_lep2_Isol_LooseTrackOnly;   //!
   TBranch        *b_lep2_Isol_Loose;   //!
   TBranch        *b_lep2_Isol_Gradient;   //!
   TBranch        *b_lep2_Isol_GradientLoose;   //!
   TBranch        *b_lep2_Isol_FixedCutTightTrackOnly;   //!
   TBranch        *b_lep2_Isol_FixedCutLoose;   //!
   TBranch        *b_lep2_Isol_FixedCutTight;   //!
   TBranch        *b_lep3_Isol_LooseTrackOnly;   //!
   TBranch        *b_lep3_Isol_Loose;   //!
   TBranch        *b_lep3_Isol_Gradient;   //!
   TBranch        *b_lep3_Isol_GradientLoose;   //!
   TBranch        *b_lep3_Isol_FixedCutTightTrackOnly;   //!
   TBranch        *b_lep3_Isol_FixedCutLoose;   //!
   TBranch        *b_lep3_Isol_FixedCutTight;   //!
   TBranch        *b_hasTrigMatchedEl;   //!
   TBranch        *b_hasTrigMatchedMu;   //!
   TBranch        *b_TopHeavyFlavorFilterFlag;   //!
   TBranch        *b_MEgamma;   //!
   TBranch        *b_numberTruePhotons;   //!
   TBranch        *b_MTW;   //!
   TBranch        *b_el_ECID;   //!
   TBranch        *b_OSSF_M;   //!
   TBranch        *b_nOSSF;   //!
   TBranch        *b_nOS;   //!
   TBranch        *b_mll;   //!
   TBranch        *b_prompt;   //!
   TBranch        *b_xs_weight;   //!

   nominal(TTree *tree=0);
   virtual ~nominal();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif