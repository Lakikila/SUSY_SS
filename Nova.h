//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul 18 15:26:46 2021 by ROOT version 6.22/08
// from TTree trees_SRRPV_/trees_SRRPV_
// found on file: 403615_new.root
//////////////////////////////////////////////////////////

#ifndef Nova_h
#define Nova_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class Nova {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           NPV;
   vector<int>     *TriggerDecisions;
   Float_t         actualInteractionsPerCrossing;
   Float_t         averageInteractionsPerCrossing;
   UInt_t          bcid;
   ULong64_t       eventNumber;
   vector<float>   *fatjet_C2;
   vector<float>   *fatjet_D2;
   vector<float>   *fatjet_ECF1;
   vector<float>   *fatjet_ECF2;
   vector<float>   *fatjet_ECF3;
   vector<float>   *fatjet_NTrimSubjets;
   vector<float>   *fatjet_Split12;
   vector<float>   *fatjet_Split23;
   vector<float>   *fatjet_Split34;
   vector<double>  *fatjet_e;
   vector<double>  *fatjet_eta;
   vector<double>  *fatjet_phi;
   vector<double>  *fatjet_pt;
   vector<float>   *fatjet_tau1_wta;
   vector<float>   *fatjet_tau21_wta;
   vector<float>   *fatjet_tau2_wta;
   vector<float>   *fatjet_tau32_wta;
   vector<float>   *fatjet_tau3_wta;
   vector<float>   *fatjet_ungrtrk500;
   vector<int>     *jet_PartonTruthLabelID;
   vector<int>     *jet_bTag;
   vector<double>  *jet_deltaR0_20_matched_truth_particle_barcode;
   vector<double>  *jet_e;
   vector<double>  *jet_eta;
   vector<int>     *jet_isSig;
   vector<int>     *jet_passOR;
   vector<double>  *jet_phi;
   vector<double>  *jet_pt;
   UInt_t          lumiBlock;
   UInt_t          mcChannelNumber;
   Float_t         mcEventWeight;
   Float_t         minmax_jet_mass;
   Int_t           pass_HLT_ht1000_L1J100;
   Int_t           pass_HLT_ht700_L1J75;
   Long64_t        pileupReweightHash;
   Float_t         pileupWeight;
   UInt_t          runNumber;
   vector<int>     *truth_NeutralinoFromGluino_ParentBarcode;
   vector<int>     *truth_NeutralinoFromGluino_barcode;
   vector<double>  *truth_NeutralinoFromGluino_charge;
   vector<double>  *truth_NeutralinoFromGluino_e;
   vector<double>  *truth_NeutralinoFromGluino_eta;
   vector<int>     *truth_NeutralinoFromGluino_pdgID;
   vector<double>  *truth_NeutralinoFromGluino_phi;
   vector<double>  *truth_NeutralinoFromGluino_pt;
   vector<int>     *truth_QuarkFromGluino_ParentBarcode;
   vector<int>     *truth_QuarkFromGluino_barcode;
   vector<double>  *truth_QuarkFromGluino_charge;
   vector<double>  *truth_QuarkFromGluino_e;
   vector<double>  *truth_QuarkFromGluino_eta;
   vector<int>     *truth_QuarkFromGluino_pdgID;
   vector<double>  *truth_QuarkFromGluino_phi;
   vector<double>  *truth_QuarkFromGluino_pt;
   vector<int>     *truth_QuarkFromNeutralino_ParentBarcode;
   vector<int>     *truth_QuarkFromNeutralino_barcode;
   vector<double>  *truth_QuarkFromNeutralino_charge;
   vector<double>  *truth_QuarkFromNeutralino_e;
   vector<double>  *truth_QuarkFromNeutralino_eta;
   vector<int>     *truth_QuarkFromNeutralino_pdgID;
   vector<double>  *truth_QuarkFromNeutralino_phi;
   vector<double>  *truth_QuarkFromNeutralino_pt;
   vector<double>  *truth_jet_e;
   vector<double>  *truth_jet_eta;
   vector<double>  *truth_jet_phi;
   vector<double>  *truth_jet_pt;
   vector<double>  *truth_parent__charge;
   vector<int>     *truth_parent_barcode;
   vector<double>  *truth_parent_e;
   vector<double>  *truth_parent_eta;
   vector<double>  *truth_parent_m;
   vector<int>     *truth_parent_pdgId;
   vector<double>  *truth_parent_phi;
   vector<double>  *truth_parent_pt;
   Float_t         normweight;

   // List of branches
   TBranch        *b_NPV;   //!
   TBranch        *b_TriggerDecisions;   //!
   TBranch        *b_actualInteractionsPerCrossing;   //!
   TBranch        *b_averageInteractionsPerCrossing;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_fatjet_C2;   //!
   TBranch        *b_fatjet_D2;   //!
   TBranch        *b_fatjet_ECF1;   //!
   TBranch        *b_fatjet_ECF2;   //!
   TBranch        *b_fatjet_ECF3;   //!
   TBranch        *b_fatjet_NTrimSubjets;   //!
   TBranch        *b_fatjet_Split12;   //!
   TBranch        *b_fatjet_Split23;   //!
   TBranch        *b_fatjet_Split34;   //!
   TBranch        *b_fatjet_e;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_tau1_wta;   //!
   TBranch        *b_fatjet_tau21_wta;   //!
   TBranch        *b_fatjet_tau2_wta;   //!
   TBranch        *b_fatjet_tau32_wta;   //!
   TBranch        *b_fatjet_tau3_wta;   //!
   TBranch        *b_fatjet_ungrtrk500;   //!
   TBranch        *b_jet_PartonTruthLabelID;   //!
   TBranch        *b_jet_bTag;   //!
   TBranch        *b_jet_deltaR0_20_matched_truth_particle_barcode;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_isSig;   //!
   TBranch        *b_jet_passOR;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_mcEventWeight;   //!
   TBranch        *b_minmax_jet_mass;   //!
   TBranch        *b_pass_HLT_ht1000_L1J100;   //!
   TBranch        *b_pass_HLT_ht700_L1J75;   //!
   TBranch        *b_pileupReweightHash;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_truth_NeutralinoFromGluino_ParentBarcode;   //!
   TBranch        *b_truth_NeutralinoFromGluino_barcode;   //!
   TBranch        *b_truth_NeutralinoFromGluino_charge;   //!
   TBranch        *b_truth_NeutralinoFromGluino_e;   //!
   TBranch        *b_truth_NeutralinoFromGluino_eta;   //!
   TBranch        *b_truth_NeutralinoFromGluino_pdgID;   //!
   TBranch        *b_truth_NeutralinoFromGluino_phi;   //!
   TBranch        *b_truth_NeutralinoFromGluino_pt;   //!
   TBranch        *b_truth_QuarkFromGluino_ParentBarcode;   //!
   TBranch        *b_truth_QuarkFromGluino_barcode;   //!
   TBranch        *b_truth_QuarkFromGluino_charge;   //!
   TBranch        *b_truth_QuarkFromGluino_e;   //!
   TBranch        *b_truth_QuarkFromGluino_eta;   //!
   TBranch        *b_truth_QuarkFromGluino_pdgID;   //!
   TBranch        *b_truth_QuarkFromGluino_phi;   //!
   TBranch        *b_truth_QuarkFromGluino_pt;   //!
   TBranch        *b_truth_QuarkFromNeutralino_ParentBarcode;   //!
   TBranch        *b_truth_QuarkFromNeutralino_barcode;   //!
   TBranch        *b_truth_QuarkFromNeutralino_charge;   //!
   TBranch        *b_truth_QuarkFromNeutralino_e;   //!
   TBranch        *b_truth_QuarkFromNeutralino_eta;   //!
   TBranch        *b_truth_QuarkFromNeutralino_pdgID;   //!
   TBranch        *b_truth_QuarkFromNeutralino_phi;   //!
   TBranch        *b_truth_QuarkFromNeutralino_pt;   //!
   TBranch        *b_truth_jet_e;   //!
   TBranch        *b_truth_jet_eta;   //!
   TBranch        *b_truth_jet_phi;   //!
   TBranch        *b_truth_jet_pt;   //!
   TBranch        *b_truth_parent__charge;   //!
   TBranch        *b_truth_parent_barcode;   //!
   TBranch        *b_truth_parent_e;   //!
   TBranch        *b_truth_parent_eta;   //!
   TBranch        *b_truth_parent_m;   //!
   TBranch        *b_truth_parent_pdgId;   //!
   TBranch        *b_truth_parent_phi;   //!
   TBranch        *b_truth_parent_pt;   //!
   TBranch        *b_normweight;   //!

   Nova(TTree *tree=0);
   virtual ~Nova();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Nova_cxx
Nova::Nova(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("403615_new.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("403615_new.root");
      }
      f->GetObject("trees_SRRPV_",tree);

   }
   Init(tree);
}

Nova::~Nova()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Nova::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Nova::LoadTree(Long64_t entry)
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

void Nova::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TriggerDecisions = 0;
   fatjet_C2 = 0;
   fatjet_D2 = 0;
   fatjet_ECF1 = 0;
   fatjet_ECF2 = 0;
   fatjet_ECF3 = 0;
   fatjet_NTrimSubjets = 0;
   fatjet_Split12 = 0;
   fatjet_Split23 = 0;
   fatjet_Split34 = 0;
   fatjet_e = 0;
   fatjet_eta = 0;
   fatjet_phi = 0;
   fatjet_pt = 0;
   fatjet_tau1_wta = 0;
   fatjet_tau21_wta = 0;
   fatjet_tau2_wta = 0;
   fatjet_tau32_wta = 0;
   fatjet_tau3_wta = 0;
   fatjet_ungrtrk500 = 0;
   jet_PartonTruthLabelID = 0;
   jet_bTag = 0;
   jet_deltaR0_20_matched_truth_particle_barcode = 0;
   jet_e = 0;
   jet_eta = 0;
   jet_isSig = 0;
   jet_passOR = 0;
   jet_phi = 0;
   jet_pt = 0;
   truth_NeutralinoFromGluino_ParentBarcode = 0;
   truth_NeutralinoFromGluino_barcode = 0;
   truth_NeutralinoFromGluino_charge = 0;
   truth_NeutralinoFromGluino_e = 0;
   truth_NeutralinoFromGluino_eta = 0;
   truth_NeutralinoFromGluino_pdgID = 0;
   truth_NeutralinoFromGluino_phi = 0;
   truth_NeutralinoFromGluino_pt = 0;
   truth_QuarkFromGluino_ParentBarcode = 0;
   truth_QuarkFromGluino_barcode = 0;
   truth_QuarkFromGluino_charge = 0;
   truth_QuarkFromGluino_e = 0;
   truth_QuarkFromGluino_eta = 0;
   truth_QuarkFromGluino_pdgID = 0;
   truth_QuarkFromGluino_phi = 0;
   truth_QuarkFromGluino_pt = 0;
   truth_QuarkFromNeutralino_ParentBarcode = 0;
   truth_QuarkFromNeutralino_barcode = 0;
   truth_QuarkFromNeutralino_charge = 0;
   truth_QuarkFromNeutralino_e = 0;
   truth_QuarkFromNeutralino_eta = 0;
   truth_QuarkFromNeutralino_pdgID = 0;
   truth_QuarkFromNeutralino_phi = 0;
   truth_QuarkFromNeutralino_pt = 0;
   truth_jet_e = 0;
   truth_jet_eta = 0;
   truth_jet_phi = 0;
   truth_jet_pt = 0;
   truth_parent__charge = 0;
   truth_parent_barcode = 0;
   truth_parent_e = 0;
   truth_parent_eta = 0;
   truth_parent_m = 0;
   truth_parent_pdgId = 0;
   truth_parent_phi = 0;
   truth_parent_pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
   fChain->SetBranchAddress("TriggerDecisions", &TriggerDecisions, &b_TriggerDecisions);
   fChain->SetBranchAddress("actualInteractionsPerCrossing", &actualInteractionsPerCrossing, &b_actualInteractionsPerCrossing);
   fChain->SetBranchAddress("averageInteractionsPerCrossing", &averageInteractionsPerCrossing, &b_averageInteractionsPerCrossing);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("fatjet_C2", &fatjet_C2, &b_fatjet_C2);
   fChain->SetBranchAddress("fatjet_D2", &fatjet_D2, &b_fatjet_D2);
   fChain->SetBranchAddress("fatjet_ECF1", &fatjet_ECF1, &b_fatjet_ECF1);
   fChain->SetBranchAddress("fatjet_ECF2", &fatjet_ECF2, &b_fatjet_ECF2);
   fChain->SetBranchAddress("fatjet_ECF3", &fatjet_ECF3, &b_fatjet_ECF3);
   fChain->SetBranchAddress("fatjet_NTrimSubjets", &fatjet_NTrimSubjets, &b_fatjet_NTrimSubjets);
   fChain->SetBranchAddress("fatjet_Split12", &fatjet_Split12, &b_fatjet_Split12);
   fChain->SetBranchAddress("fatjet_Split23", &fatjet_Split23, &b_fatjet_Split23);
   fChain->SetBranchAddress("fatjet_Split34", &fatjet_Split34, &b_fatjet_Split34);
   fChain->SetBranchAddress("fatjet_e", &fatjet_e, &b_fatjet_e);
   fChain->SetBranchAddress("fatjet_eta", &fatjet_eta, &b_fatjet_eta);
   fChain->SetBranchAddress("fatjet_phi", &fatjet_phi, &b_fatjet_phi);
   fChain->SetBranchAddress("fatjet_pt", &fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_tau1_wta", &fatjet_tau1_wta, &b_fatjet_tau1_wta);
   fChain->SetBranchAddress("fatjet_tau21_wta", &fatjet_tau21_wta, &b_fatjet_tau21_wta);
   fChain->SetBranchAddress("fatjet_tau2_wta", &fatjet_tau2_wta, &b_fatjet_tau2_wta);
   fChain->SetBranchAddress("fatjet_tau32_wta", &fatjet_tau32_wta, &b_fatjet_tau32_wta);
   fChain->SetBranchAddress("fatjet_tau3_wta", &fatjet_tau3_wta, &b_fatjet_tau3_wta);
   fChain->SetBranchAddress("fatjet_ungrtrk500", &fatjet_ungrtrk500, &b_fatjet_ungrtrk500);
   fChain->SetBranchAddress("jet_PartonTruthLabelID", &jet_PartonTruthLabelID, &b_jet_PartonTruthLabelID);
   fChain->SetBranchAddress("jet_bTag", &jet_bTag, &b_jet_bTag);
   fChain->SetBranchAddress("jet_deltaR0.20_matched_truth_particle_barcode", &jet_deltaR0_20_matched_truth_particle_barcode, &b_jet_deltaR0_20_matched_truth_particle_barcode);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_isSig", &jet_isSig, &b_jet_isSig);
   fChain->SetBranchAddress("jet_passOR", &jet_passOR, &b_jet_passOR);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("mcEventWeight", &mcEventWeight, &b_mcEventWeight);
   fChain->SetBranchAddress("minmax_jet_mass", &minmax_jet_mass, &b_minmax_jet_mass);
   fChain->SetBranchAddress("pass_HLT_ht1000_L1J100", &pass_HLT_ht1000_L1J100, &b_pass_HLT_ht1000_L1J100);
   fChain->SetBranchAddress("pass_HLT_ht700_L1J75", &pass_HLT_ht700_L1J75, &b_pass_HLT_ht700_L1J75);
   fChain->SetBranchAddress("pileupReweightHash", &pileupReweightHash, &b_pileupReweightHash);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_ParentBarcode", &truth_NeutralinoFromGluino_ParentBarcode, &b_truth_NeutralinoFromGluino_ParentBarcode);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_barcode", &truth_NeutralinoFromGluino_barcode, &b_truth_NeutralinoFromGluino_barcode);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_charge", &truth_NeutralinoFromGluino_charge, &b_truth_NeutralinoFromGluino_charge);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_e", &truth_NeutralinoFromGluino_e, &b_truth_NeutralinoFromGluino_e);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_eta", &truth_NeutralinoFromGluino_eta, &b_truth_NeutralinoFromGluino_eta);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_pdgID", &truth_NeutralinoFromGluino_pdgID, &b_truth_NeutralinoFromGluino_pdgID);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_phi", &truth_NeutralinoFromGluino_phi, &b_truth_NeutralinoFromGluino_phi);
   fChain->SetBranchAddress("truth_NeutralinoFromGluino_pt", &truth_NeutralinoFromGluino_pt, &b_truth_NeutralinoFromGluino_pt);
   fChain->SetBranchAddress("truth_QuarkFromGluino_ParentBarcode", &truth_QuarkFromGluino_ParentBarcode, &b_truth_QuarkFromGluino_ParentBarcode);
   fChain->SetBranchAddress("truth_QuarkFromGluino_barcode", &truth_QuarkFromGluino_barcode, &b_truth_QuarkFromGluino_barcode);
   fChain->SetBranchAddress("truth_QuarkFromGluino_charge", &truth_QuarkFromGluino_charge, &b_truth_QuarkFromGluino_charge);
   fChain->SetBranchAddress("truth_QuarkFromGluino_e", &truth_QuarkFromGluino_e, &b_truth_QuarkFromGluino_e);
   fChain->SetBranchAddress("truth_QuarkFromGluino_eta", &truth_QuarkFromGluino_eta, &b_truth_QuarkFromGluino_eta);
   fChain->SetBranchAddress("truth_QuarkFromGluino_pdgID", &truth_QuarkFromGluino_pdgID, &b_truth_QuarkFromGluino_pdgID);
   fChain->SetBranchAddress("truth_QuarkFromGluino_phi", &truth_QuarkFromGluino_phi, &b_truth_QuarkFromGluino_phi);
   fChain->SetBranchAddress("truth_QuarkFromGluino_pt", &truth_QuarkFromGluino_pt, &b_truth_QuarkFromGluino_pt);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_ParentBarcode", &truth_QuarkFromNeutralino_ParentBarcode, &b_truth_QuarkFromNeutralino_ParentBarcode);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_barcode", &truth_QuarkFromNeutralino_barcode, &b_truth_QuarkFromNeutralino_barcode);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_charge", &truth_QuarkFromNeutralino_charge, &b_truth_QuarkFromNeutralino_charge);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_e", &truth_QuarkFromNeutralino_e, &b_truth_QuarkFromNeutralino_e);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_eta", &truth_QuarkFromNeutralino_eta, &b_truth_QuarkFromNeutralino_eta);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_pdgID", &truth_QuarkFromNeutralino_pdgID, &b_truth_QuarkFromNeutralino_pdgID);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_phi", &truth_QuarkFromNeutralino_phi, &b_truth_QuarkFromNeutralino_phi);
   fChain->SetBranchAddress("truth_QuarkFromNeutralino_pt", &truth_QuarkFromNeutralino_pt, &b_truth_QuarkFromNeutralino_pt);
   fChain->SetBranchAddress("truth_jet_e", &truth_jet_e, &b_truth_jet_e);
   fChain->SetBranchAddress("truth_jet_eta", &truth_jet_eta, &b_truth_jet_eta);
   fChain->SetBranchAddress("truth_jet_phi", &truth_jet_phi, &b_truth_jet_phi);
   fChain->SetBranchAddress("truth_jet_pt", &truth_jet_pt, &b_truth_jet_pt);
   fChain->SetBranchAddress("truth_parent__charge", &truth_parent__charge, &b_truth_parent__charge);
   fChain->SetBranchAddress("truth_parent_barcode", &truth_parent_barcode, &b_truth_parent_barcode);
   fChain->SetBranchAddress("truth_parent_e", &truth_parent_e, &b_truth_parent_e);
   fChain->SetBranchAddress("truth_parent_eta", &truth_parent_eta, &b_truth_parent_eta);
   fChain->SetBranchAddress("truth_parent_m", &truth_parent_m, &b_truth_parent_m);
   fChain->SetBranchAddress("truth_parent_pdgId", &truth_parent_pdgId, &b_truth_parent_pdgId);
   fChain->SetBranchAddress("truth_parent_phi", &truth_parent_phi, &b_truth_parent_phi);
   fChain->SetBranchAddress("truth_parent_pt", &truth_parent_pt, &b_truth_parent_pt);
   fChain->SetBranchAddress("normweight", &normweight, &b_normweight);
   Notify();
}

Bool_t Nova::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Nova::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Nova::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Nova_cxx
