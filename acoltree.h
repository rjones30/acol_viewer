//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  4 11:50:45 2014 by ROOT version 5.34/21
// from TTree N9:raw_XP/PXI records
// found on file: rawdata/ac_20141101_181632.root
//////////////////////////////////////////////////////////

#ifndef acoltree_h
#define acoltree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class acoltree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *fSummary;

   // Declaration of leaf types
   struct record_t {
      Long64_t tsec;
      Long64_t tnsec;
      Float_t  data[8192];
      Float_t  current;
   };
   struct record_t record;

   // List of branches
   TBranch        *b_record;

   acoltree(TTree *tree=0);
   virtual ~acoltree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef acoltree_cxx
acoltree::acoltree(TTree *tree) : fChain(0), fSummary(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TString rootfile("rawdata/ac_20141030_052216.root");
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(rootfile);
      if (!f || !f->IsOpen()) {
         f = new TFile(rootfile);
      }
      f->GetObject("N9:raw_XP", tree);
   }
   Init(tree);
}

acoltree::~acoltree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t acoltree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t acoltree::LoadTree(Long64_t entry)
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

void acoltree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("record", &record.tsec, &b_record);
   Notify();
}

Bool_t acoltree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void acoltree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t acoltree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return (entry >= 0)? 1 : 0;
}

#endif // #ifdef acoltree_cxx
