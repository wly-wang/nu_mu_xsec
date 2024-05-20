//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 12 21:21:51 2024 by ROOT version 6.28/06
// from TTree T_PFeval/T_PFeval
// found on file: ./ROOT_files/checkout_run3_overlaycv.root
//////////////////////////////////////////////////////////

#ifndef overlay_pfeval_run3_h
#define overlay_pfeval_run3_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "TObjArray.h"

class overlay_pfeval_run3 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Int_t           neutrino_type;
   Float_t         reco_nuvtxX;
   Float_t         reco_nuvtxY;
   Float_t         reco_nuvtxZ;
   vector<float>   *reco_vec_showervtxX;
   vector<float>   *reco_vec_showervtxY;
   vector<float>   *reco_vec_showervtxZ;
   vector<float>   *reco_vec_showerKE;
   Float_t         reco_showervtxX;
   Float_t         reco_showervtxY;
   Float_t         reco_showervtxZ;
   Float_t         reco_showerKE;
   Float_t         reco_showerMomentum[4];
   Float_t         reco_muonvtxX;
   Float_t         reco_muonvtxY;
   Float_t         reco_muonvtxZ;
   Float_t         reco_muonMomentum[4];
   Int_t           reco_Nproton;
   Float_t         reco_protonvtxX;
   Float_t         reco_protonvtxY;
   Float_t         reco_protonvtxZ;
   Float_t         reco_protonMomentum[4];
   Float_t         nuvtx_diff;
   Float_t         showervtx_diff;
   Float_t         muonvtx_diff;
   Int_t           mcflux_run;
   Int_t           mcflux_evtno;
   Int_t           mcflux_ndecay;
   Int_t           mcflux_ntype;
   Int_t           mcflux_ptype;
   Int_t           mcflux_tptype;
   Float_t         mcflux_nuEnergy;
   Float_t         mcflux_vx;
   Float_t         mcflux_vy;
   Float_t         mcflux_vz;
   Float_t         mcflux_genx;
   Float_t         mcflux_geny;
   Float_t         mcflux_genz;
   Float_t         mcflux_dk2gen;
   Float_t         mcflux_gen2vtx;
   Float_t         truth_corr_nuvtxX;
   Float_t         truth_corr_nuvtxY;
   Float_t         truth_corr_nuvtxZ;
   Float_t         truth_corr_showervtxX;
   Float_t         truth_corr_showervtxY;
   Float_t         truth_corr_showervtxZ;
   Float_t         truth_showerKE;
   Float_t         truth_showerMomentum[4];
   Int_t           truth_showerPdg;
   Int_t           truth_showerMother;
   Float_t         truth_corr_muonvtxX;
   Float_t         truth_corr_muonvtxY;
   Float_t         truth_corr_muonvtxZ;
   Float_t         truth_muonvtxX;
   Float_t         truth_muonvtxY;
   Float_t         truth_muonvtxZ;
   Float_t         truth_muonendX;
   Float_t         truth_muonendY;
   Float_t         truth_muonendZ;
   Float_t         truth_muonMomentum[4];
   Float_t         truth_nuEnergy;
   Float_t         truth_energyInside;
   Float_t         truth_electronInside;
   Int_t           truth_nuPdg;
   Bool_t          truth_isCC;
   Float_t         truth_vtxX;
   Float_t         truth_vtxY;
   Float_t         truth_vtxZ;
   Float_t         truth_nuTime;
   Int_t           truth_nuIntType;
   Int_t           truth_nuScatType;
   Int_t           truth_Npi0;
   Int_t           truth_NprimPio;
   Float_t         truth_pio_energy_1;
   Float_t         truth_pio_energy_2;
   Float_t         truth_pio_angle;
   Int_t           truth_NCDelta;
   Int_t           truth_single_photon;
   Float_t         truth_photon_angle;
   Float_t         truth_photon_dis;
   Float_t         truth_nu_pos[4];
   Float_t         truth_nu_momentum[4];
   Int_t           truth_Ntrack;
   Int_t           truth_id[1708];   //[truth_Ntrack]
   Int_t           truth_pdg[1708];   //[truth_Ntrack]
   vector<string>  *truth_process;
   Int_t           truth_mother[1708];   //[truth_Ntrack]
   Float_t         truth_startXYZT[1708][4];   //[truth_Ntrack]
   Float_t         truth_endXYZT[1708][4];   //[truth_Ntrack]
   Float_t         truth_startMomentum[1708][4];   //[truth_Ntrack]
   Float_t         truth_endMomentum[1708][4];   //[truth_Ntrack]
   vector<vector<int> > *truth_daughters;
   TObjArray       *fMC_trackPosition;
   Int_t           mc_isnu;
   Int_t           mc_nGeniePrimaries;
   Int_t           mc_nu_pdg;
   Int_t           mc_nu_ccnc;
   Int_t           mc_nu_mode;
   Int_t           mc_nu_intType;
   Int_t           mc_nu_target;
   Int_t           mc_hitnuc;
   Int_t           mc_hitquark;
   Double_t        mc_nu_Q2;
   Double_t        mc_nu_W;
   Double_t        mc_nu_X;
   Double_t        mc_nu_Y;
   Double_t        mc_nu_Pt;
   Double_t        mc_nu_Theta;
   Float_t         mc_nu_pos[4];
   Float_t         mc_nu_mom[4];
   Int_t           reco_Ntrack;
   Int_t           reco_id[253];   //[reco_Ntrack]
   Int_t           reco_pdg[253];   //[reco_Ntrack]
   vector<string>  *reco_process;
   Int_t           reco_mother[253];   //[reco_Ntrack]
   Float_t         reco_startXYZT[253][4];   //[reco_Ntrack]
   Float_t         reco_endXYZT[253][4];   //[reco_Ntrack]
   Float_t         reco_startMomentum[253][4];   //[reco_Ntrack]
   Float_t         reco_endMomentum[253][4];   //[reco_Ntrack]
   vector<vector<int> > *reco_daughters;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_neutrino_type;   //!
   TBranch        *b_reco_nuvtxX;   //!
   TBranch        *b_reco_nuvtxY;   //!
   TBranch        *b_reco_nuvtxZ;   //!
   TBranch        *b_reco_vec_showervtxX;   //!
   TBranch        *b_reco_vec_showervtxY;   //!
   TBranch        *b_reco_vec_showervtxZ;   //!
   TBranch        *b_reco_vec_showerKE;   //!
   TBranch        *b_reco_showervtxX;   //!
   TBranch        *b_reco_showervtxY;   //!
   TBranch        *b_reco_showervtxZ;   //!
   TBranch        *b_reco_showerKE;   //!
   TBranch        *b_reco_showerMomentum;   //!
   TBranch        *b_reco_muonvtxX;   //!
   TBranch        *b_reco_muonvtxY;   //!
   TBranch        *b_reco_muonvtxZ;   //!
   TBranch        *b_reco_muonMomentum;   //!
   TBranch        *b_reco_Nproton;   //!
   TBranch        *b_reco_protonvtxX;   //!
   TBranch        *b_reco_protonvtxY;   //!
   TBranch        *b_reco_protonvtxZ;   //!
   TBranch        *b_reco_protonMomentum;   //!
   TBranch        *b_nuvtx_diff;   //!
   TBranch        *b_showervtx_diff;   //!
   TBranch        *b_muonvtx_diff;   //!
   TBranch        *b_mcflux_run;   //!
   TBranch        *b_mcflux_evtno;   //!
   TBranch        *b_mcflux_ndecay;   //!
   TBranch        *b_mcflux_ntype;   //!
   TBranch        *b_mcflux_ptype;   //!
   TBranch        *b_mcflux_tptype;   //!
   TBranch        *b_mcflux_nuEnergy;   //!
   TBranch        *b_mcflux_vx;   //!
   TBranch        *b_mcflux_vy;   //!
   TBranch        *b_mcflux_vz;   //!
   TBranch        *b_mcflux_genx;   //!
   TBranch        *b_mcflux_geny;   //!
   TBranch        *b_mcflux_genz;   //!
   TBranch        *b_mcflux_dk2gen;   //!
   TBranch        *b_mcflux_gen2vtx;   //!
   TBranch        *b_truth_corr_nuvtxX;   //!
   TBranch        *b_truth_corr_nuvtxY;   //!
   TBranch        *b_truth_corr_nuvtxZ;   //!
   TBranch        *b_truth_corr_showervtxX;   //!
   TBranch        *b_truth_corr_showervtxY;   //!
   TBranch        *b_truth_corr_showervtxZ;   //!
   TBranch        *b_truth_showerKE;   //!
   TBranch        *b_truth_showerMomentum;   //!
   TBranch        *b_truth_showerPdg;   //!
   TBranch        *b_truth_showerMother;   //!
   TBranch        *b_truth_corr_muonvtxX;   //!
   TBranch        *b_truth_corr_muonvtxY;   //!
   TBranch        *b_truth_corr_muonvtxZ;   //!
   TBranch        *b_truth_muonvtxX;   //!
   TBranch        *b_truth_muonvtxY;   //!
   TBranch        *b_truth_muonvtxZ;   //!
   TBranch        *b_truth_muonendX;   //!
   TBranch        *b_truth_muonendY;   //!
   TBranch        *b_truth_muonendZ;   //!
   TBranch        *b_truth_muonMomentum;   //!
   TBranch        *b_truth_nuEnergy;   //!
   TBranch        *b_truth_energyInside;   //!
   TBranch        *b_truth_electronInside;   //!
   TBranch        *b_truth_nuPdg;   //!
   TBranch        *b_truth_isCC;   //!
   TBranch        *b_truth_vtxX;   //!
   TBranch        *b_truth_vtxY;   //!
   TBranch        *b_truth_vtxZ;   //!
   TBranch        *b_truth_nuTime;   //!
   TBranch        *b_truth_nuIntType;   //!
   TBranch        *b_truth_nuScatType;   //!
   TBranch        *b_truth_Npi0;   //!
   TBranch        *b_truth_NprimPio;   //!
   TBranch        *b_truth_pio_energy_1;   //!
   TBranch        *b_truth_pio_energy_2;   //!
   TBranch        *b_truth_pio_angle;   //!
   TBranch        *b_truth_NCDelta;   //!
   TBranch        *b_truth_single_photon;   //!
   TBranch        *b_truth_photon_angle;   //!
   TBranch        *b_truth_photon_dis;   //!
   TBranch        *b_truth_nu_pos;   //!
   TBranch        *b_truth_nu_momentum;   //!
   TBranch        *b_truth_Ntrack;   //!
   TBranch        *b_truth_id;   //!
   TBranch        *b_truth_pdg;   //!
   TBranch        *b_truth_process;   //!
   TBranch        *b_truth_mother;   //!
   TBranch        *b_truth_startXYZT;   //!
   TBranch        *b_truth_endXYZT;   //!
   TBranch        *b_truth_startMomentum;   //!
   TBranch        *b_truth_endMomentum;   //!
   TBranch        *b_truth_daughters;   //!
   TBranch        *b_fMC_trackPosition;   //!
   TBranch        *b_mc_isnu;   //!
   TBranch        *b_mc_nGeniePrimaries;   //!
   TBranch        *b_mc_nu_pdg;   //!
   TBranch        *b_mc_nu_ccnc;   //!
   TBranch        *b_mc_nu_mode;   //!
   TBranch        *b_mc_nu_intType;   //!
   TBranch        *b_mc_nu_target;   //!
   TBranch        *b_mc_hitnuc;   //!
   TBranch        *b_mc_hitquark;   //!
   TBranch        *b_mc_nu_Q2;   //!
   TBranch        *b_mc_nu_W;   //!
   TBranch        *b_mc_nu_X;   //!
   TBranch        *b_mc_nu_Y;   //!
   TBranch        *b_mc_nu_Pt;   //!
   TBranch        *b_mc_nu_Theta;   //!
   TBranch        *b_mc_nu_pos;   //!
   TBranch        *b_mc_nu_mom;   //!
   TBranch        *b_reco_Ntrack;   //!
   TBranch        *b_reco_id;   //!
   TBranch        *b_reco_pdg;   //!
   TBranch        *b_reco_process;   //!
   TBranch        *b_reco_mother;   //!
   TBranch        *b_reco_startXYZT;   //!
   TBranch        *b_reco_endXYZT;   //!
   TBranch        *b_reco_startMomentum;   //!
   TBranch        *b_reco_endMomentum;   //!
   TBranch        *b_reco_daughters;   //!

   overlay_pfeval_run3(TTree *tree=0);
   virtual ~overlay_pfeval_run3();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef overlay_pfeval_run3_cxx
overlay_pfeval_run3::overlay_pfeval_run3(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./ROOT_files/checkout_run3_overlaycv.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./ROOT_files/checkout_run3_overlaycv.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("./ROOT_files/checkout_run3_overlaycv.root:/wcpselection");
      dir->GetObject("T_PFeval",tree);

   }
   Init(tree);
}

overlay_pfeval_run3::~overlay_pfeval_run3()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t overlay_pfeval_run3::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t overlay_pfeval_run3::LoadTree(Long64_t entry)
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

void overlay_pfeval_run3::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   reco_vec_showervtxX = 0;
   reco_vec_showervtxY = 0;
   reco_vec_showervtxZ = 0;
   reco_vec_showerKE = 0;
   truth_process = 0;
   truth_daughters = 0;
   fMC_trackPosition = 0;
   reco_process = 0;
   reco_daughters = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("neutrino_type", &neutrino_type, &b_neutrino_type);
   fChain->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX, &b_reco_nuvtxX);
   fChain->SetBranchAddress("reco_nuvtxY", &reco_nuvtxY, &b_reco_nuvtxY);
   fChain->SetBranchAddress("reco_nuvtxZ", &reco_nuvtxZ, &b_reco_nuvtxZ);
   fChain->SetBranchAddress("reco_vec_showervtxX", &reco_vec_showervtxX, &b_reco_vec_showervtxX);
   fChain->SetBranchAddress("reco_vec_showervtxY", &reco_vec_showervtxY, &b_reco_vec_showervtxY);
   fChain->SetBranchAddress("reco_vec_showervtxZ", &reco_vec_showervtxZ, &b_reco_vec_showervtxZ);
   fChain->SetBranchAddress("reco_vec_showerKE", &reco_vec_showerKE, &b_reco_vec_showerKE);
   fChain->SetBranchAddress("reco_showervtxX", &reco_showervtxX, &b_reco_showervtxX);
   fChain->SetBranchAddress("reco_showervtxY", &reco_showervtxY, &b_reco_showervtxY);
   fChain->SetBranchAddress("reco_showervtxZ", &reco_showervtxZ, &b_reco_showervtxZ);
   fChain->SetBranchAddress("reco_showerKE", &reco_showerKE, &b_reco_showerKE);
   fChain->SetBranchAddress("reco_showerMomentum", reco_showerMomentum, &b_reco_showerMomentum);
   fChain->SetBranchAddress("reco_muonvtxX", &reco_muonvtxX, &b_reco_muonvtxX);
   fChain->SetBranchAddress("reco_muonvtxY", &reco_muonvtxY, &b_reco_muonvtxY);
   fChain->SetBranchAddress("reco_muonvtxZ", &reco_muonvtxZ, &b_reco_muonvtxZ);
   fChain->SetBranchAddress("reco_muonMomentum", reco_muonMomentum, &b_reco_muonMomentum);
   fChain->SetBranchAddress("reco_Nproton", &reco_Nproton, &b_reco_Nproton);
   fChain->SetBranchAddress("reco_protonvtxX", &reco_protonvtxX, &b_reco_protonvtxX);
   fChain->SetBranchAddress("reco_protonvtxY", &reco_protonvtxY, &b_reco_protonvtxY);
   fChain->SetBranchAddress("reco_protonvtxZ", &reco_protonvtxZ, &b_reco_protonvtxZ);
   fChain->SetBranchAddress("reco_protonMomentum", reco_protonMomentum, &b_reco_protonMomentum);
   fChain->SetBranchAddress("nuvtx_diff", &nuvtx_diff, &b_nuvtx_diff);
   fChain->SetBranchAddress("showervtx_diff", &showervtx_diff, &b_showervtx_diff);
   fChain->SetBranchAddress("muonvtx_diff", &muonvtx_diff, &b_muonvtx_diff);
   fChain->SetBranchAddress("mcflux_run", &mcflux_run, &b_mcflux_run);
   fChain->SetBranchAddress("mcflux_evtno", &mcflux_evtno, &b_mcflux_evtno);
   fChain->SetBranchAddress("mcflux_ndecay", &mcflux_ndecay, &b_mcflux_ndecay);
   fChain->SetBranchAddress("mcflux_ntype", &mcflux_ntype, &b_mcflux_ntype);
   fChain->SetBranchAddress("mcflux_ptype", &mcflux_ptype, &b_mcflux_ptype);
   fChain->SetBranchAddress("mcflux_tptype", &mcflux_tptype, &b_mcflux_tptype);
   fChain->SetBranchAddress("mcflux_nuEnergy", &mcflux_nuEnergy, &b_mcflux_nuEnergy);
   fChain->SetBranchAddress("mcflux_vx", &mcflux_vx, &b_mcflux_vx);
   fChain->SetBranchAddress("mcflux_vy", &mcflux_vy, &b_mcflux_vy);
   fChain->SetBranchAddress("mcflux_vz", &mcflux_vz, &b_mcflux_vz);
   fChain->SetBranchAddress("mcflux_genx", &mcflux_genx, &b_mcflux_genx);
   fChain->SetBranchAddress("mcflux_geny", &mcflux_geny, &b_mcflux_geny);
   fChain->SetBranchAddress("mcflux_genz", &mcflux_genz, &b_mcflux_genz);
   fChain->SetBranchAddress("mcflux_dk2gen", &mcflux_dk2gen, &b_mcflux_dk2gen);
   fChain->SetBranchAddress("mcflux_gen2vtx", &mcflux_gen2vtx, &b_mcflux_gen2vtx);
   fChain->SetBranchAddress("truth_corr_nuvtxX", &truth_corr_nuvtxX, &b_truth_corr_nuvtxX);
   fChain->SetBranchAddress("truth_corr_nuvtxY", &truth_corr_nuvtxY, &b_truth_corr_nuvtxY);
   fChain->SetBranchAddress("truth_corr_nuvtxZ", &truth_corr_nuvtxZ, &b_truth_corr_nuvtxZ);
   fChain->SetBranchAddress("truth_corr_showervtxX", &truth_corr_showervtxX, &b_truth_corr_showervtxX);
   fChain->SetBranchAddress("truth_corr_showervtxY", &truth_corr_showervtxY, &b_truth_corr_showervtxY);
   fChain->SetBranchAddress("truth_corr_showervtxZ", &truth_corr_showervtxZ, &b_truth_corr_showervtxZ);
   fChain->SetBranchAddress("truth_showerKE", &truth_showerKE, &b_truth_showerKE);
   fChain->SetBranchAddress("truth_showerMomentum", truth_showerMomentum, &b_truth_showerMomentum);
   fChain->SetBranchAddress("truth_showerPdg", &truth_showerPdg, &b_truth_showerPdg);
   fChain->SetBranchAddress("truth_showerMother", &truth_showerMother, &b_truth_showerMother);
   fChain->SetBranchAddress("truth_corr_muonvtxX", &truth_corr_muonvtxX, &b_truth_corr_muonvtxX);
   fChain->SetBranchAddress("truth_corr_muonvtxY", &truth_corr_muonvtxY, &b_truth_corr_muonvtxY);
   fChain->SetBranchAddress("truth_corr_muonvtxZ", &truth_corr_muonvtxZ, &b_truth_corr_muonvtxZ);
   fChain->SetBranchAddress("truth_muonvtxX", &truth_muonvtxX, &b_truth_muonvtxX);
   fChain->SetBranchAddress("truth_muonvtxY", &truth_muonvtxY, &b_truth_muonvtxY);
   fChain->SetBranchAddress("truth_muonvtxZ", &truth_muonvtxZ, &b_truth_muonvtxZ);
   fChain->SetBranchAddress("truth_muonendX", &truth_muonendX, &b_truth_muonendX);
   fChain->SetBranchAddress("truth_muonendY", &truth_muonendY, &b_truth_muonendY);
   fChain->SetBranchAddress("truth_muonendZ", &truth_muonendZ, &b_truth_muonendZ);
   fChain->SetBranchAddress("truth_muonMomentum", truth_muonMomentum, &b_truth_muonMomentum);
   fChain->SetBranchAddress("truth_nuEnergy", &truth_nuEnergy, &b_truth_nuEnergy);
   fChain->SetBranchAddress("truth_energyInside", &truth_energyInside, &b_truth_energyInside);
   fChain->SetBranchAddress("truth_electronInside", &truth_electronInside, &b_truth_electronInside);
   fChain->SetBranchAddress("truth_nuPdg", &truth_nuPdg, &b_truth_nuPdg);
   fChain->SetBranchAddress("truth_isCC", &truth_isCC, &b_truth_isCC);
   fChain->SetBranchAddress("truth_vtxX", &truth_vtxX, &b_truth_vtxX);
   fChain->SetBranchAddress("truth_vtxY", &truth_vtxY, &b_truth_vtxY);
   fChain->SetBranchAddress("truth_vtxZ", &truth_vtxZ, &b_truth_vtxZ);
   fChain->SetBranchAddress("truth_nuTime", &truth_nuTime, &b_truth_nuTime);
   fChain->SetBranchAddress("truth_nuIntType", &truth_nuIntType, &b_truth_nuIntType);
   fChain->SetBranchAddress("truth_nuScatType", &truth_nuScatType, &b_truth_nuScatType);
   fChain->SetBranchAddress("truth_Npi0", &truth_Npi0, &b_truth_Npi0);
   fChain->SetBranchAddress("truth_NprimPio", &truth_NprimPio, &b_truth_NprimPio);
   fChain->SetBranchAddress("truth_pio_energy_1", &truth_pio_energy_1, &b_truth_pio_energy_1);
   fChain->SetBranchAddress("truth_pio_energy_2", &truth_pio_energy_2, &b_truth_pio_energy_2);
   fChain->SetBranchAddress("truth_pio_angle", &truth_pio_angle, &b_truth_pio_angle);
   fChain->SetBranchAddress("truth_NCDelta", &truth_NCDelta, &b_truth_NCDelta);
   fChain->SetBranchAddress("truth_single_photon", &truth_single_photon, &b_truth_single_photon);
   fChain->SetBranchAddress("truth_photon_angle", &truth_photon_angle, &b_truth_photon_angle);
   fChain->SetBranchAddress("truth_photon_dis", &truth_photon_dis, &b_truth_photon_dis);
   fChain->SetBranchAddress("truth_nu_pos", truth_nu_pos, &b_truth_nu_pos);
   fChain->SetBranchAddress("truth_nu_momentum", truth_nu_momentum, &b_truth_nu_momentum);
   fChain->SetBranchAddress("truth_Ntrack", &truth_Ntrack, &b_truth_Ntrack);
   fChain->SetBranchAddress("truth_id", truth_id, &b_truth_id);
   fChain->SetBranchAddress("truth_pdg", truth_pdg, &b_truth_pdg);
   fChain->SetBranchAddress("truth_process", &truth_process, &b_truth_process);
   fChain->SetBranchAddress("truth_mother", truth_mother, &b_truth_mother);
   fChain->SetBranchAddress("truth_startXYZT", truth_startXYZT, &b_truth_startXYZT);
   fChain->SetBranchAddress("truth_endXYZT", truth_endXYZT, &b_truth_endXYZT);
   fChain->SetBranchAddress("truth_startMomentum", truth_startMomentum, &b_truth_startMomentum);
   fChain->SetBranchAddress("truth_endMomentum", truth_endMomentum, &b_truth_endMomentum);
   fChain->SetBranchAddress("truth_daughters", &truth_daughters, &b_truth_daughters);
   fChain->SetBranchAddress("fMC_trackPosition", &fMC_trackPosition, &b_fMC_trackPosition);
   fChain->SetBranchAddress("mc_isnu", &mc_isnu, &b_mc_isnu);
   fChain->SetBranchAddress("mc_nGeniePrimaries", &mc_nGeniePrimaries, &b_mc_nGeniePrimaries);
   fChain->SetBranchAddress("mc_nu_pdg", &mc_nu_pdg, &b_mc_nu_pdg);
   fChain->SetBranchAddress("mc_nu_ccnc", &mc_nu_ccnc, &b_mc_nu_ccnc);
   fChain->SetBranchAddress("mc_nu_mode", &mc_nu_mode, &b_mc_nu_mode);
   fChain->SetBranchAddress("mc_nu_intType", &mc_nu_intType, &b_mc_nu_intType);
   fChain->SetBranchAddress("mc_nu_target", &mc_nu_target, &b_mc_nu_target);
   fChain->SetBranchAddress("mc_hitnuc", &mc_hitnuc, &b_mc_hitnuc);
   fChain->SetBranchAddress("mc_hitquark", &mc_hitquark, &b_mc_hitquark);
   fChain->SetBranchAddress("mc_nu_Q2", &mc_nu_Q2, &b_mc_nu_Q2);
   fChain->SetBranchAddress("mc_nu_W", &mc_nu_W, &b_mc_nu_W);
   fChain->SetBranchAddress("mc_nu_X", &mc_nu_X, &b_mc_nu_X);
   fChain->SetBranchAddress("mc_nu_Y", &mc_nu_Y, &b_mc_nu_Y);
   fChain->SetBranchAddress("mc_nu_Pt", &mc_nu_Pt, &b_mc_nu_Pt);
   fChain->SetBranchAddress("mc_nu_Theta", &mc_nu_Theta, &b_mc_nu_Theta);
   fChain->SetBranchAddress("mc_nu_pos", mc_nu_pos, &b_mc_nu_pos);
   fChain->SetBranchAddress("mc_nu_mom", mc_nu_mom, &b_mc_nu_mom);
   fChain->SetBranchAddress("reco_Ntrack", &reco_Ntrack, &b_reco_Ntrack);
   fChain->SetBranchAddress("reco_id", reco_id, &b_reco_id);
   fChain->SetBranchAddress("reco_pdg", reco_pdg, &b_reco_pdg);
   fChain->SetBranchAddress("reco_process", &reco_process, &b_reco_process);
   fChain->SetBranchAddress("reco_mother", reco_mother, &b_reco_mother);
   fChain->SetBranchAddress("reco_startXYZT", reco_startXYZT, &b_reco_startXYZT);
   fChain->SetBranchAddress("reco_endXYZT", reco_endXYZT, &b_reco_endXYZT);
   fChain->SetBranchAddress("reco_startMomentum", reco_startMomentum, &b_reco_startMomentum);
   fChain->SetBranchAddress("reco_endMomentum", reco_endMomentum, &b_reco_endMomentum);
   fChain->SetBranchAddress("reco_daughters", &reco_daughters, &b_reco_daughters);
   Notify();
}

Bool_t overlay_pfeval_run3::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void overlay_pfeval_run3::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t overlay_pfeval_run3::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef overlay_pfeval_run3_cxx
