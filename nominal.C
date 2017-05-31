#define nominal_cxx
#include "nominal.h"
#include <TStyle.h>
#include <TCanvas.h>

void nominal::Loop()
{
//   In a ROOT session, you can do:
//      root> .L nominal.C
//      root> nominal t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

nominal::nominal(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Signal_v9ta.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Signal_v9ta.root");
      }
      f->GetObject("nominal",tree);

   }
   Init(tree);
}

nominal::~nominal()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t nominal::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t nominal::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void nominal::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PDFinfo_X1 = 0;
   PDFinfo_X2 = 0;
   PDFinfo_PDGID1 = 0;
   PDFinfo_PDGID2 = 0;
   PDFinfo_Q = 0;
   PDFinfo_XF1 = 0;
   PDFinfo_XF2 = 0;
   weight_bTagSF_77_eigenvars_B_up = 0;
   weight_bTagSF_77_eigenvars_C_up = 0;
   weight_bTagSF_77_eigenvars_Light_up = 0;
   weight_bTagSF_77_eigenvars_B_down = 0;
   weight_bTagSF_77_eigenvars_C_down = 0;
   weight_bTagSF_77_eigenvars_Light_down = 0;
   el_pt = 0;
   el_eta = 0;
   el_cl_eta = 0;
   el_phi = 0;
   el_e = 0;
   el_charge = 0;
   el_topoetcone20 = 0;
   el_ptvarcone20 = 0;
   el_CF = 0;
   el_d0sig = 0;
   el_delta_z0_sintheta = 0;
   el_true_type = 0;
   el_true_origin = 0;
   el_true_typebkg = 0;
   el_true_originbkg = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_e = 0;
   mu_charge = 0;
   mu_topoetcone20 = 0;
   mu_ptvarcone30 = 0;
   mu_d0sig = 0;
   mu_delta_z0_sintheta = 0;
   mu_true_type = 0;
   mu_true_origin = 0;
   tau_pt = 0;
   tau_eta = 0;
   tau_phi = 0;
   tau_charge = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_e = 0;
   jet_mv2c00 = 0;
   jet_mv2c10 = 0;
   jet_mv2c20 = 0;
   jet_ip3dsv1 = 0;
   jet_jvt = 0;
   jet_passfjvt = 0;
   jet_truthflav = 0;
   jet_truthPartonLabel = 0;
   jet_isTrueHS = 0;
   jet_isbtagged_77 = 0;
   el_trigMatch_HLT_e60_lhmedium_nod0 = 0;
   el_trigMatch_HLT_e26_lhtight_nod0_ivarloose = 0;
   el_trigMatch_HLT_e140_lhloose_nod0 = 0;
   el_trigMatch_HLT_e60_lhmedium = 0;
   el_trigMatch_HLT_e24_lhmedium_L1EM20VH = 0;
   el_trigMatch_HLT_e120_lhloose = 0;
   mu_trigMatch_HLT_mu26_ivarmedium = 0;
   mu_trigMatch_HLT_mu50 = 0;
   mu_trigMatch_HLT_mu20_iloose_L1MU15 = 0;
   MTW = 0;
   el_ECID = 0;
   OSSF_M = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PDFinfo_X1", &PDFinfo_X1, &b_PDFinfo_X1);
   fChain->SetBranchAddress("PDFinfo_X2", &PDFinfo_X2, &b_PDFinfo_X2);
   fChain->SetBranchAddress("PDFinfo_PDGID1", &PDFinfo_PDGID1, &b_PDFinfo_PDGID1);
   fChain->SetBranchAddress("PDFinfo_PDGID2", &PDFinfo_PDGID2, &b_PDFinfo_PDGID2);
   fChain->SetBranchAddress("PDFinfo_Q", &PDFinfo_Q, &b_PDFinfo_Q);
   fChain->SetBranchAddress("PDFinfo_XF1", &PDFinfo_XF1, &b_PDFinfo_XF1);
   fChain->SetBranchAddress("PDFinfo_XF2", &PDFinfo_XF2, &b_PDFinfo_XF2);
   fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   fChain->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
   fChain->SetBranchAddress("weight_leptonSF", &weight_leptonSF, &b_weight_leptonSF);
   fChain->SetBranchAddress("weight_tauSF", &weight_tauSF, &b_weight_tauSF);
   fChain->SetBranchAddress("weight_bTagSF_77", &weight_bTagSF_77, &b_weight_bTagSF_77);
   fChain->SetBranchAddress("weight_jvt", &weight_jvt, &b_weight_jvt);
   fChain->SetBranchAddress("weight_pileup_UP", &weight_pileup_UP, &b_weight_pileup_UP);
   fChain->SetBranchAddress("weight_pileup_DOWN", &weight_pileup_DOWN, &b_weight_pileup_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Trigger_UP", &weight_leptonSF_EL_SF_Trigger_UP, &b_weight_leptonSF_EL_SF_Trigger_UP);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Trigger_DOWN", &weight_leptonSF_EL_SF_Trigger_DOWN, &b_weight_leptonSF_EL_SF_Trigger_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Reco_UP", &weight_leptonSF_EL_SF_Reco_UP, &b_weight_leptonSF_EL_SF_Reco_UP);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Reco_DOWN", &weight_leptonSF_EL_SF_Reco_DOWN, &b_weight_leptonSF_EL_SF_Reco_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_ID_UP", &weight_leptonSF_EL_SF_ID_UP, &b_weight_leptonSF_EL_SF_ID_UP);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_ID_DOWN", &weight_leptonSF_EL_SF_ID_DOWN, &b_weight_leptonSF_EL_SF_ID_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Isol_UP", &weight_leptonSF_EL_SF_Isol_UP, &b_weight_leptonSF_EL_SF_Isol_UP);
   fChain->SetBranchAddress("weight_leptonSF_EL_SF_Isol_DOWN", &weight_leptonSF_EL_SF_Isol_DOWN, &b_weight_leptonSF_EL_SF_Isol_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_STAT_UP", &weight_leptonSF_MU_SF_Trigger_STAT_UP, &b_weight_leptonSF_MU_SF_Trigger_STAT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_STAT_DOWN", &weight_leptonSF_MU_SF_Trigger_STAT_DOWN, &b_weight_leptonSF_MU_SF_Trigger_STAT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_SYST_UP", &weight_leptonSF_MU_SF_Trigger_SYST_UP, &b_weight_leptonSF_MU_SF_Trigger_SYST_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_SYST_DOWN", &weight_leptonSF_MU_SF_Trigger_SYST_DOWN, &b_weight_leptonSF_MU_SF_Trigger_SYST_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_UP", &weight_leptonSF_MU_SF_ID_STAT_UP, &b_weight_leptonSF_MU_SF_ID_STAT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_DOWN", &weight_leptonSF_MU_SF_ID_STAT_DOWN, &b_weight_leptonSF_MU_SF_ID_STAT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_UP", &weight_leptonSF_MU_SF_ID_SYST_UP, &b_weight_leptonSF_MU_SF_ID_SYST_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_DOWN", &weight_leptonSF_MU_SF_ID_SYST_DOWN, &b_weight_leptonSF_MU_SF_ID_SYST_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP", &weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP, &b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN", &weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN, &b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP", &weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP, &b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN", &weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN, &b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_STAT_UP", &weight_leptonSF_MU_SF_Isol_STAT_UP, &b_weight_leptonSF_MU_SF_Isol_STAT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_STAT_DOWN", &weight_leptonSF_MU_SF_Isol_STAT_DOWN, &b_weight_leptonSF_MU_SF_Isol_STAT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_SYST_UP", &weight_leptonSF_MU_SF_Isol_SYST_UP, &b_weight_leptonSF_MU_SF_Isol_SYST_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_SYST_DOWN", &weight_leptonSF_MU_SF_Isol_SYST_DOWN, &b_weight_leptonSF_MU_SF_Isol_SYST_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_STAT_UP", &weight_leptonSF_MU_SF_TTVA_STAT_UP, &b_weight_leptonSF_MU_SF_TTVA_STAT_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_STAT_DOWN", &weight_leptonSF_MU_SF_TTVA_STAT_DOWN, &b_weight_leptonSF_MU_SF_TTVA_STAT_DOWN);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_SYST_UP", &weight_leptonSF_MU_SF_TTVA_SYST_UP, &b_weight_leptonSF_MU_SF_TTVA_SYST_UP);
   fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_SYST_DOWN", &weight_leptonSF_MU_SF_TTVA_SYST_DOWN, &b_weight_leptonSF_MU_SF_TTVA_SYST_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger", &weight_indiv_SF_EL_Trigger, &b_weight_indiv_SF_EL_Trigger);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger_UP", &weight_indiv_SF_EL_Trigger_UP, &b_weight_indiv_SF_EL_Trigger_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger_DOWN", &weight_indiv_SF_EL_Trigger_DOWN, &b_weight_indiv_SF_EL_Trigger_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco", &weight_indiv_SF_EL_Reco, &b_weight_indiv_SF_EL_Reco);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_UP", &weight_indiv_SF_EL_Reco_UP, &b_weight_indiv_SF_EL_Reco_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_DOWN", &weight_indiv_SF_EL_Reco_DOWN, &b_weight_indiv_SF_EL_Reco_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID", &weight_indiv_SF_EL_ID, &b_weight_indiv_SF_EL_ID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID_UP", &weight_indiv_SF_EL_ID_UP, &b_weight_indiv_SF_EL_ID_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID_DOWN", &weight_indiv_SF_EL_ID_DOWN, &b_weight_indiv_SF_EL_ID_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol", &weight_indiv_SF_EL_Isol, &b_weight_indiv_SF_EL_Isol);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_UP", &weight_indiv_SF_EL_Isol_UP, &b_weight_indiv_SF_EL_Isol_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_DOWN", &weight_indiv_SF_EL_Isol_DOWN, &b_weight_indiv_SF_EL_Isol_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID", &weight_indiv_SF_EL_ChargeID, &b_weight_indiv_SF_EL_ChargeID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID_UP", &weight_indiv_SF_EL_ChargeID_UP, &b_weight_indiv_SF_EL_ChargeID_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID_DOWN", &weight_indiv_SF_EL_ChargeID_DOWN, &b_weight_indiv_SF_EL_ChargeID_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID", &weight_indiv_SF_EL_ChargeMisID, &b_weight_indiv_SF_EL_ChargeMisID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_STAT_UP", &weight_indiv_SF_EL_ChargeMisID_STAT_UP, &b_weight_indiv_SF_EL_ChargeMisID_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_STAT_DOWN", &weight_indiv_SF_EL_ChargeMisID_STAT_DOWN, &b_weight_indiv_SF_EL_ChargeMisID_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_SYST_UP", &weight_indiv_SF_EL_ChargeMisID_SYST_UP, &b_weight_indiv_SF_EL_ChargeMisID_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_SYST_DOWN", &weight_indiv_SF_EL_ChargeMisID_SYST_DOWN, &b_weight_indiv_SF_EL_ChargeMisID_SYST_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Trigger", &weight_indiv_SF_MU_Trigger, &b_weight_indiv_SF_MU_Trigger);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Trigger_STAT_UP", &weight_indiv_SF_MU_Trigger_STAT_UP, &b_weight_indiv_SF_MU_Trigger_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Trigger_STAT_DOWN", &weight_indiv_SF_MU_Trigger_STAT_DOWN, &b_weight_indiv_SF_MU_Trigger_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Trigger_SYST_UP", &weight_indiv_SF_MU_Trigger_SYST_UP, &b_weight_indiv_SF_MU_Trigger_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Trigger_SYST_DOWN", &weight_indiv_SF_MU_Trigger_SYST_DOWN, &b_weight_indiv_SF_MU_Trigger_SYST_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID", &weight_indiv_SF_MU_ID, &b_weight_indiv_SF_MU_ID);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_UP", &weight_indiv_SF_MU_ID_STAT_UP, &b_weight_indiv_SF_MU_ID_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_DOWN", &weight_indiv_SF_MU_ID_STAT_DOWN, &b_weight_indiv_SF_MU_ID_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_UP", &weight_indiv_SF_MU_ID_SYST_UP, &b_weight_indiv_SF_MU_ID_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_DOWN", &weight_indiv_SF_MU_ID_SYST_DOWN, &b_weight_indiv_SF_MU_ID_SYST_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_LOWPT_UP", &weight_indiv_SF_MU_ID_STAT_LOWPT_UP, &b_weight_indiv_SF_MU_ID_STAT_LOWPT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN", &weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN, &b_weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_LOWPT_UP", &weight_indiv_SF_MU_ID_SYST_LOWPT_UP, &b_weight_indiv_SF_MU_ID_SYST_LOWPT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN", &weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN, &b_weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Isol", &weight_indiv_SF_MU_Isol, &b_weight_indiv_SF_MU_Isol);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_STAT_UP", &weight_indiv_SF_MU_Isol_STAT_UP, &b_weight_indiv_SF_MU_Isol_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_STAT_DOWN", &weight_indiv_SF_MU_Isol_STAT_DOWN, &b_weight_indiv_SF_MU_Isol_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_SYST_UP", &weight_indiv_SF_MU_Isol_SYST_UP, &b_weight_indiv_SF_MU_Isol_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_SYST_DOWN", &weight_indiv_SF_MU_Isol_SYST_DOWN, &b_weight_indiv_SF_MU_Isol_SYST_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA", &weight_indiv_SF_MU_TTVA, &b_weight_indiv_SF_MU_TTVA);
   fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_STAT_UP", &weight_indiv_SF_MU_TTVA_STAT_UP, &b_weight_indiv_SF_MU_TTVA_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_STAT_DOWN", &weight_indiv_SF_MU_TTVA_STAT_DOWN, &b_weight_indiv_SF_MU_TTVA_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_SYST_UP", &weight_indiv_SF_MU_TTVA_SYST_UP, &b_weight_indiv_SF_MU_TTVA_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_SYST_DOWN", &weight_indiv_SF_MU_TTVA_SYST_DOWN, &b_weight_indiv_SF_MU_TTVA_SYST_DOWN);
   fChain->SetBranchAddress("weight_tauSF_ELEOLR_UP", &weight_tauSF_ELEOLR_UP, &b_weight_tauSF_ELEOLR_UP);
   fChain->SetBranchAddress("weight_tauSF_ELEOLR_DOWN", &weight_tauSF_ELEOLR_DOWN, &b_weight_tauSF_ELEOLR_DOWN);
   fChain->SetBranchAddress("weight_tauSF_JETID_UP", &weight_tauSF_JETID_UP, &b_weight_tauSF_JETID_UP);
   fChain->SetBranchAddress("weight_tauSF_JETID_DOWN", &weight_tauSF_JETID_DOWN, &b_weight_tauSF_JETID_DOWN);
   fChain->SetBranchAddress("weight_tauSF_RECO_UP", &weight_tauSF_RECO_UP, &b_weight_tauSF_RECO_UP);
   fChain->SetBranchAddress("weight_tauSF_RECO_DOWN", &weight_tauSF_RECO_DOWN, &b_weight_tauSF_RECO_DOWN);
   fChain->SetBranchAddress("weight_jvt_UP", &weight_jvt_UP, &b_weight_jvt_UP);
   fChain->SetBranchAddress("weight_jvt_DOWN", &weight_jvt_DOWN, &b_weight_jvt_DOWN);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_B_up", &weight_bTagSF_77_eigenvars_B_up, &b_weight_bTagSF_77_eigenvars_B_up);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_C_up", &weight_bTagSF_77_eigenvars_C_up, &b_weight_bTagSF_77_eigenvars_C_up);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_Light_up", &weight_bTagSF_77_eigenvars_Light_up, &b_weight_bTagSF_77_eigenvars_Light_up);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_B_down", &weight_bTagSF_77_eigenvars_B_down, &b_weight_bTagSF_77_eigenvars_B_down);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_C_down", &weight_bTagSF_77_eigenvars_C_down, &b_weight_bTagSF_77_eigenvars_C_down);
   fChain->SetBranchAddress("weight_bTagSF_77_eigenvars_Light_down", &weight_bTagSF_77_eigenvars_Light_down, &b_weight_bTagSF_77_eigenvars_Light_down);
   fChain->SetBranchAddress("weight_bTagSF_77_extrapolation_up", &weight_bTagSF_77_extrapolation_up, &b_weight_bTagSF_77_extrapolation_up);
   fChain->SetBranchAddress("weight_bTagSF_77_extrapolation_down", &weight_bTagSF_77_extrapolation_down, &b_weight_bTagSF_77_extrapolation_down);
   fChain->SetBranchAddress("weight_bTagSF_77_extrapolation_from_charm_up", &weight_bTagSF_77_extrapolation_from_charm_up, &b_weight_bTagSF_77_extrapolation_from_charm_up);
   fChain->SetBranchAddress("weight_bTagSF_77_extrapolation_from_charm_down", &weight_bTagSF_77_extrapolation_from_charm_down, &b_weight_bTagSF_77_extrapolation_from_charm_down);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("randomRunNumber", &randomRunNumber, &b_randomRunNumber);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("backgroundFlags", &backgroundFlags, &b_backgroundFlags);
   fChain->SetBranchAddress("hasBadMuon", &hasBadMuon, &b_hasBadMuon);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_cl_eta", &el_cl_eta, &b_el_cl_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_e", &el_e, &b_el_e);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_topoetcone20", &el_topoetcone20, &b_el_topoetcone20);
   fChain->SetBranchAddress("el_ptvarcone20", &el_ptvarcone20, &b_el_ptvarcone20);
   fChain->SetBranchAddress("el_CF", &el_CF, &b_el_CF);
   fChain->SetBranchAddress("el_d0sig", &el_d0sig, &b_el_d0sig);
   fChain->SetBranchAddress("el_delta_z0_sintheta", &el_delta_z0_sintheta, &b_el_delta_z0_sintheta);
   fChain->SetBranchAddress("el_true_type", &el_true_type, &b_el_true_type);
   fChain->SetBranchAddress("el_true_origin", &el_true_origin, &b_el_true_origin);
   fChain->SetBranchAddress("el_true_typebkg", &el_true_typebkg, &b_el_true_typebkg);
   fChain->SetBranchAddress("el_true_originbkg", &el_true_originbkg, &b_el_true_originbkg);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_topoetcone20", &mu_topoetcone20, &b_mu_topoetcone20);
   fChain->SetBranchAddress("mu_ptvarcone30", &mu_ptvarcone30, &b_mu_ptvarcone30);
   fChain->SetBranchAddress("mu_d0sig", &mu_d0sig, &b_mu_d0sig);
   fChain->SetBranchAddress("mu_delta_z0_sintheta", &mu_delta_z0_sintheta, &b_mu_delta_z0_sintheta);
   fChain->SetBranchAddress("mu_true_type", &mu_true_type, &b_mu_true_type);
   fChain->SetBranchAddress("mu_true_origin", &mu_true_origin, &b_mu_true_origin);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_mv2c00", &jet_mv2c00, &b_jet_mv2c00);
   fChain->SetBranchAddress("jet_mv2c10", &jet_mv2c10, &b_jet_mv2c10);
   fChain->SetBranchAddress("jet_mv2c20", &jet_mv2c20, &b_jet_mv2c20);
   fChain->SetBranchAddress("jet_ip3dsv1", &jet_ip3dsv1, &b_jet_ip3dsv1);
   fChain->SetBranchAddress("jet_jvt", &jet_jvt, &b_jet_jvt);
   fChain->SetBranchAddress("jet_passfjvt", &jet_passfjvt, &b_jet_passfjvt);
   fChain->SetBranchAddress("jet_truthflav", &jet_truthflav, &b_jet_truthflav);
   fChain->SetBranchAddress("jet_truthPartonLabel", &jet_truthPartonLabel, &b_jet_truthPartonLabel);
   fChain->SetBranchAddress("jet_isTrueHS", &jet_isTrueHS, &b_jet_isTrueHS);
   fChain->SetBranchAddress("jet_isbtagged_77", &jet_isbtagged_77, &b_jet_isbtagged_77);
   fChain->SetBranchAddress("met_met", &met_met, &b_met_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("Selection2015", &Selection2015, &b_Selection2015);
   fChain->SetBranchAddress("Selection2016", &Selection2016, &b_Selection2016);
   fChain->SetBranchAddress("HLT_e60_lhmedium_nod0", &HLT_e60_lhmedium_nod0, &b_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("HLT_mu26_ivarmedium", &HLT_mu26_ivarmedium, &b_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("HLT_e26_lhtight_nod0_ivarloose", &HLT_e26_lhtight_nod0_ivarloose, &b_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("HLT_e140_lhloose_nod0", &HLT_e140_lhloose_nod0, &b_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("HLT_mu20_iloose_L1MU15", &HLT_mu20_iloose_L1MU15, &b_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
   fChain->SetBranchAddress("HLT_e60_lhmedium", &HLT_e60_lhmedium, &b_HLT_e60_lhmedium);
   fChain->SetBranchAddress("HLT_e24_lhmedium_L1EM20VH", &HLT_e24_lhmedium_L1EM20VH, &b_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("HLT_e120_lhloose", &HLT_e120_lhloose, &b_HLT_e120_lhloose);
   fChain->SetBranchAddress("el_trigMatch_HLT_e60_lhmedium_nod0", &el_trigMatch_HLT_e60_lhmedium_nod0, &b_el_trigMatch_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("el_trigMatch_HLT_e26_lhtight_nod0_ivarloose", &el_trigMatch_HLT_e26_lhtight_nod0_ivarloose, &b_el_trigMatch_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("el_trigMatch_HLT_e140_lhloose_nod0", &el_trigMatch_HLT_e140_lhloose_nod0, &b_el_trigMatch_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("el_trigMatch_HLT_e60_lhmedium", &el_trigMatch_HLT_e60_lhmedium, &b_el_trigMatch_HLT_e60_lhmedium);
   fChain->SetBranchAddress("el_trigMatch_HLT_e24_lhmedium_L1EM20VH", &el_trigMatch_HLT_e24_lhmedium_L1EM20VH, &b_el_trigMatch_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("el_trigMatch_HLT_e120_lhloose", &el_trigMatch_HLT_e120_lhloose, &b_el_trigMatch_HLT_e120_lhloose);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu26_ivarmedium", &mu_trigMatch_HLT_mu26_ivarmedium, &b_mu_trigMatch_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu50", &mu_trigMatch_HLT_mu50, &b_mu_trigMatch_HLT_mu50);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu20_iloose_L1MU15", &mu_trigMatch_HLT_mu20_iloose_L1MU15, &b_mu_trigMatch_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("nEl", &nEl, &b_nEl);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nBJets70", &nBJets70, &b_nBJets70);
   fChain->SetBranchAddress("nBJets77", &nBJets77, &b_nBJets77);
   fChain->SetBranchAddress("nBJets85", &nBJets85, &b_nBJets85);
   fChain->SetBranchAddress("pT_1lep", &pT_1lep, &b_pT_1lep);
   fChain->SetBranchAddress("pT_2lep", &pT_2lep, &b_pT_2lep);
   fChain->SetBranchAddress("pT_3lep", &pT_3lep, &b_pT_3lep);
   fChain->SetBranchAddress("eta_1lep", &eta_1lep, &b_eta_1lep);
   fChain->SetBranchAddress("eta_2lep", &eta_2lep, &b_eta_2lep);
   fChain->SetBranchAddress("eta_3lep", &eta_3lep, &b_eta_3lep);
   fChain->SetBranchAddress("phi_1lep", &phi_1lep, &b_phi_1lep);
   fChain->SetBranchAddress("phi_2lep", &phi_2lep, &b_phi_2lep);
   fChain->SetBranchAddress("phi_3lep", &phi_3lep, &b_phi_3lep);
   fChain->SetBranchAddress("pdgid_1lep", &pdgid_1lep, &b_pdgid_1lep);
   fChain->SetBranchAddress("pdgid_2lep", &pdgid_2lep, &b_pdgid_2lep);
   fChain->SetBranchAddress("pdgid_3lep", &pdgid_3lep, &b_pdgid_3lep);
   fChain->SetBranchAddress("charge_1lep", &charge_1lep, &b_charge_1lep);
   fChain->SetBranchAddress("charge_2lep", &charge_2lep, &b_charge_2lep);
   fChain->SetBranchAddress("charge_3lep", &charge_3lep, &b_charge_3lep);
   fChain->SetBranchAddress("pT_1jet", &pT_1jet, &b_pT_1jet);
   fChain->SetBranchAddress("pT_2jet", &pT_2jet, &b_pT_2jet);
   fChain->SetBranchAddress("pT_3jet", &pT_3jet, &b_pT_3jet);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lep1_origin", &lep1_origin, &b_lep1_origin);
   fChain->SetBranchAddress("lep2_origin", &lep2_origin, &b_lep2_origin);
   fChain->SetBranchAddress("lep3_origin", &lep3_origin, &b_lep3_origin);
   fChain->SetBranchAddress("lep1_pdgorigin", &lep1_pdgorigin, &b_lep1_pdgorigin);
   fChain->SetBranchAddress("lep2_pdgorigin", &lep2_pdgorigin, &b_lep2_pdgorigin);
   fChain->SetBranchAddress("lep3_pdgorigin", &lep3_pdgorigin, &b_lep3_pdgorigin);
   fChain->SetBranchAddress("lep1_istight", &lep1_istight, &b_lep1_istight);
   fChain->SetBranchAddress("lep2_istight", &lep2_istight, &b_lep2_istight);
   fChain->SetBranchAddress("lep3_istight", &lep3_istight, &b_lep3_istight);
   fChain->SetBranchAddress("lep1_type", &lep1_type, &b_lep1_type);
   fChain->SetBranchAddress("lep2_type", &lep2_type, &b_lep2_type);
   fChain->SetBranchAddress("lep3_type", &lep3_type, &b_lep3_type);
   fChain->SetBranchAddress("lep1_d0sig", &lep1_d0sig, &b_lep1_d0sig);
   fChain->SetBranchAddress("lep2_d0sig", &lep2_d0sig, &b_lep2_d0sig);
   fChain->SetBranchAddress("lep3_d0sig", &lep3_d0sig, &b_lep3_d0sig);
   fChain->SetBranchAddress("lep1_z0sintheta", &lep1_z0sintheta, &b_lep1_z0sintheta);
   fChain->SetBranchAddress("lep2_z0sintheta", &lep2_z0sintheta, &b_lep2_z0sintheta);
   fChain->SetBranchAddress("lep3_z0sintheta", &lep3_z0sintheta, &b_lep3_z0sintheta);
   fChain->SetBranchAddress("lep1_topoetcone", &lep1_topoetcone, &b_lep1_topoetcone);
   fChain->SetBranchAddress("lep2_topoetcone", &lep2_topoetcone, &b_lep2_topoetcone);
   fChain->SetBranchAddress("lep3_topoetcone", &lep3_topoetcone, &b_lep3_topoetcone);
   fChain->SetBranchAddress("lep1_topoetcone40", &lep1_topoetcone40, &b_lep1_topoetcone40);
   fChain->SetBranchAddress("lep2_topoetcone40", &lep2_topoetcone40, &b_lep2_topoetcone40);
   fChain->SetBranchAddress("lep3_topoetcone40", &lep3_topoetcone40, &b_lep3_topoetcone40);
   fChain->SetBranchAddress("lep1_ptvarcone", &lep1_ptvarcone, &b_lep1_ptvarcone);
   fChain->SetBranchAddress("lep2_ptvarcone", &lep2_ptvarcone, &b_lep2_ptvarcone);
   fChain->SetBranchAddress("lep3_ptvarcone", &lep3_ptvarcone, &b_lep3_ptvarcone);
   fChain->SetBranchAddress("weight0tag", &weight0tag, &b_weight0tag);
   fChain->SetBranchAddress("weight1tag", &weight1tag, &b_weight1tag);
   fChain->SetBranchAddress("weight2tag", &weight2tag, &b_weight2tag);
   fChain->SetBranchAddress("lep1_Isol_LooseTrackOnly", &lep1_Isol_LooseTrackOnly, &b_lep1_Isol_LooseTrackOnly);
   fChain->SetBranchAddress("lep1_Isol_Loose", &lep1_Isol_Loose, &b_lep1_Isol_Loose);
   fChain->SetBranchAddress("lep1_Isol_Gradient", &lep1_Isol_Gradient, &b_lep1_Isol_Gradient);
   fChain->SetBranchAddress("lep1_Isol_GradientLoose", &lep1_Isol_GradientLoose, &b_lep1_Isol_GradientLoose);
   fChain->SetBranchAddress("lep1_Isol_FixedCutTightTrackOnly", &lep1_Isol_FixedCutTightTrackOnly, &b_lep1_Isol_FixedCutTightTrackOnly);
   fChain->SetBranchAddress("lep1_Isol_FixedCutLoose", &lep1_Isol_FixedCutLoose, &b_lep1_Isol_FixedCutLoose);
   fChain->SetBranchAddress("lep1_Isol_FixedCutTight", &lep1_Isol_FixedCutTight, &b_lep1_Isol_FixedCutTight);
   fChain->SetBranchAddress("lep2_Isol_LooseTrackOnly", &lep2_Isol_LooseTrackOnly, &b_lep2_Isol_LooseTrackOnly);
   fChain->SetBranchAddress("lep2_Isol_Loose", &lep2_Isol_Loose, &b_lep2_Isol_Loose);
   fChain->SetBranchAddress("lep2_Isol_Gradient", &lep2_Isol_Gradient, &b_lep2_Isol_Gradient);
   fChain->SetBranchAddress("lep2_Isol_GradientLoose", &lep2_Isol_GradientLoose, &b_lep2_Isol_GradientLoose);
   fChain->SetBranchAddress("lep2_Isol_FixedCutTightTrackOnly", &lep2_Isol_FixedCutTightTrackOnly, &b_lep2_Isol_FixedCutTightTrackOnly);
   fChain->SetBranchAddress("lep2_Isol_FixedCutLoose", &lep2_Isol_FixedCutLoose, &b_lep2_Isol_FixedCutLoose);
   fChain->SetBranchAddress("lep2_Isol_FixedCutTight", &lep2_Isol_FixedCutTight, &b_lep2_Isol_FixedCutTight);
   fChain->SetBranchAddress("lep3_Isol_LooseTrackOnly", &lep3_Isol_LooseTrackOnly, &b_lep3_Isol_LooseTrackOnly);
   fChain->SetBranchAddress("lep3_Isol_Loose", &lep3_Isol_Loose, &b_lep3_Isol_Loose);
   fChain->SetBranchAddress("lep3_Isol_Gradient", &lep3_Isol_Gradient, &b_lep3_Isol_Gradient);
   fChain->SetBranchAddress("lep3_Isol_GradientLoose", &lep3_Isol_GradientLoose, &b_lep3_Isol_GradientLoose);
   fChain->SetBranchAddress("lep3_Isol_FixedCutTightTrackOnly", &lep3_Isol_FixedCutTightTrackOnly, &b_lep3_Isol_FixedCutTightTrackOnly);
   fChain->SetBranchAddress("lep3_Isol_FixedCutLoose", &lep3_Isol_FixedCutLoose, &b_lep3_Isol_FixedCutLoose);
   fChain->SetBranchAddress("lep3_Isol_FixedCutTight", &lep3_Isol_FixedCutTight, &b_lep3_Isol_FixedCutTight);
   fChain->SetBranchAddress("hasTrigMatchedEl", &hasTrigMatchedEl, &b_hasTrigMatchedEl);
   fChain->SetBranchAddress("hasTrigMatchedMu", &hasTrigMatchedMu, &b_hasTrigMatchedMu);
   fChain->SetBranchAddress("TopHeavyFlavorFilterFlag", &TopHeavyFlavorFilterFlag, &b_TopHeavyFlavorFilterFlag);
   fChain->SetBranchAddress("MEgamma", &MEgamma, &b_MEgamma);
   fChain->SetBranchAddress("numberTruePhotons", &numberTruePhotons, &b_numberTruePhotons);
   fChain->SetBranchAddress("MTW", &MTW, &b_MTW);
   fChain->SetBranchAddress("el_ECID", &el_ECID, &b_el_ECID);
   fChain->SetBranchAddress("OSSF_M", &OSSF_M, &b_OSSF_M);
   fChain->SetBranchAddress("nOSSF", &nOSSF, &b_nOSSF);
   fChain->SetBranchAddress("nOS", &nOS, &b_nOS);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
   fChain->SetBranchAddress("prompt", &prompt, &b_prompt);
   fChain->SetBranchAddress("xs_weight", &xs_weight, &b_xs_weight);
   Notify();
}

Bool_t nominal::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void nominal::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t nominal::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}