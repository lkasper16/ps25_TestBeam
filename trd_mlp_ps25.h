//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 17 22:57:15 2023 by ROOT version 6.14/04
// from TTree gem_hits/GEM TTree with single track hit info
// found on file: trd_singleTrackHits_Run_003200.root
//////////////////////////////////////////////////////////


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"


// Declaration of leaf types
Int_t           event_num;
Int_t           tgem_nhit;
Int_t           tgem_nclu;
Int_t           tgem_ntracks;
Int_t           mmg1_nhit;
Int_t           mmg1_nclu;
Int_t           mmg1_ntracks;
//float           ecal_energy;
//float           presh_energy;
//float           mult_energy;

vector<bool>    *parID;
vector<int>     *xpos;
vector<float>   *zpos;
vector<float>   *dedx;
vector<float>   *zHist;
vector<float>   *xposc;
vector<float>   *zposc;
vector<float>   *dedxc;
vector<float>   *widthc;

float     xposc_max;
float     zposc_max;
float     dedxc_max;
float     widthc_max;
float     dedxc_tot;

// List of branches
TBranch        *b_event_num;   //!
//TBranch        *b_ecal_energy;   //!
//TBranch        *b_presh_energy;   //!
//TBranch        *b_mult_energy;   //!
TBranch        *b_tgem_nhit;   //!
TBranch        *b_tgem_nclu;
TBranch        *b_tgem_ntracks;
TBranch        *b_mmg1_nhit;   //!
TBranch        *b_mmg1_nclu;
TBranch        *b_mmg1_ntracks;
TBranch        *b_xpos;   //!
TBranch        *b_zpos;   //!
TBranch        *b_dedx;   //!
TBranch        *b_parID;   //!
TBranch        *b_zHist;   //!
TBranch        *b_xposc;   //!
TBranch        *b_zposc;   //!
TBranch        *b_dedxc;   //!
TBranch        *b_widthc;   //!
TBranch        *b_xposc_max;
TBranch        *b_zposc_max;
TBranch        *b_dedxc_max;
TBranch        *b_widthc_max;
TBranch        *b_dedxc_tot;
