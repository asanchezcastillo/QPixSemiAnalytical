// -----------------------------------------------------------------------------
//  ROOTFileManager.cpp
//
//  Class definition of the ROOT file manager
//   * Author: Everybody is an author!
//   * Creation date: 31 August 2020
// -----------------------------------------------------------------------------

#include "RootFileManager.h"

// ROOT includes
#include "TObject.h"

// math includes
#include <math.h>

//--------------------------------------------------------------------------
ROOTFileManager::ROOTFileManager(char* inputFileName, char* outputFileName, std::string TreeName)
{
    
    Initialize(inputFileName, outputFileName, TreeName);
}

//--------------------------------------------------------------------------
ROOTFileManager::~ROOTFileManager()
{}

//--------------------------------------------------------------------------

void ROOTFileManager::Initialize(char* inputFileName, char* outputFileName, std::string TreeName)
{
    TreeName_=TreeName;
    const char* TreeName__ = TreeName_.c_str();
    //Input file:
    input_file = new TFile(inputFileName);
    input_tree = (TTree*)input_file->Get(TreeName__);
    //TFile *OutputFile = TFile::Open(outputFileName, "RECREATE");
    set_branch_addresses(input_tree);
}

//--------------------------------------------------------------------------
void ROOTFileManager::set_branch_addresses(TTree *input_tree)
  {
    if(TreeName_=="InitialParticle")
    {
    input_tree->SetBranchAddress("InitialParticlePDG", &InitialPDG_);
    input_tree->SetBranchAddress("InteractionTime", &InteractionTime_);
    input_tree->SetBranchAddress("InitialParticleEnergy", &InitialEnergy_);
    input_tree->SetBranchAddress("InteractionVertexX", &InteractionVertexX_);
    input_tree->SetBranchAddress("InteractionVertexY", &InteractionVertexY_);
    input_tree->SetBranchAddress("InteractionVertexZ", &InteractionVertexZ_);
    }
    else if(TreeName_=="Hits")
    {
    input_tree->SetBranchAddress("run",&run);
    input_tree->SetBranchAddress("hit_start_x",&hitX_start);
    input_tree->SetBranchAddress("hit_end_x",&hitX_end);
    input_tree->SetBranchAddress("hit_start_y",&hitY_start);
    input_tree->SetBranchAddress("hit_end_y",&hitY_end);
    input_tree->SetBranchAddress("hit_start_z",&hitZ_start);
    input_tree->SetBranchAddress("hit_end_z",&hitZ_end);
    input_tree->SetBranchAddress("hit_energy_deposit",&edep);
    input_tree->SetBranchAddress("hit_start_t",&time_start);
    input_tree->SetBranchAddress("hit_end_t",&time_end);
    input_tree->SetBranchAddress("hit_length",&length);
    input_tree->SetBranchAddress("hit_pdg",&pdg);
    }
    else if(TreeName_=="PrimaryParticle")
    {
    input_tree->SetBranchAddress("PrimaryParticlePDG", &PrimaryPDG_);
    input_tree->SetBranchAddress("PrimaryParticleEnergy", &PrimaryEnergy_);
    input_tree->SetBranchAddress("PrimaryParticlePx", &PrimaryPx_);
    input_tree->SetBranchAddress("PrimaryParticlePy", &PrimaryPy_);
    input_tree->SetBranchAddress("PrimaryParticlePz", &PrimaryPz_);
    }
    else if(TreeName_=="Background")
    {
    input_tree->SetBranchAddress("DecayTime", &BackgroundDecayTime_);
    input_tree->SetBranchAddress("AtomicNumber", &BackgroundAtomicNumber_);
    input_tree->SetBranchAddress("AtomicMass", &BackgroundAtomicMass_);
    input_tree->SetBranchAddress("BackgroundVertexX", &BackgroundVertexX_);
    input_tree->SetBranchAddress("BackgroundVertexY", &BackgroundVertexY_);
    input_tree->SetBranchAddress("BackgroundVertexZ", &BackgroundVertexZ_);

    }
}

void ROOTFileManager::EventFill()
{
  //  OutputTree->Branch("eventID", &runID);
   // OutputTree->Branch("SimPhotonsperOpChVUV",&SavePhotons);
  //  OutputTree->Branch("GeneratedPhotons",&generated_counter);
  //  OutputTree->Branch("DetectedPhotons",&nPhotons);
  //  OutputTree->Fill();
  //  OutputFile->cd();
}

//--------------------------------------------------------------------------
void ROOTFileManager::Save()
{
  //  OutputFile->Write();
  //  OutputFile->Close();
}

//-------------------------------------------------------------------------- 
// Get event. The number of run is in case we have several events in the same file, so we can loop over them.
void ROOTFileManager::GetEvent(int nRun)
{
  input_tree->GetEntry(nRun);
}//Get_Event

//--------------------------------------------------------------------------
//Get event in case we dont have an argument. It is for information  on the intial and primary particles
void ROOTFileManager::GetEvent()
{
    if(TreeName_=="InitialParticle")
    {
      for(int i = 0; i < input_tree->GetEntries(); i++)
      {
        input_tree->GetEntry(i);
        InitialPDG.push_back(InitialPDG_);
        InitialEnergy->push_back(InitialEnergy_);
        InteractionTime->push_back(InteractionTime_);
        InteractionVertexX->push_back(InteractionVertexX_);
        InteractionVertexY->push_back(InteractionVertexY_);
        InteractionVertexZ->push_back(InteractionVertexZ_);
      }      
  
    }
    else if(TreeName_=="PrimaryParticle") 
    {
      for(int i = 0; i < input_tree->GetEntries(); i++)
      {
        input_tree->GetEntry(i);
        PrimaryPDG.push_back(PrimaryPDG_);
        PrimaryEnergy->push_back(PrimaryEnergy_);
        PrimaryPx->push_back(PrimaryPx_);
        PrimaryPy->push_back(PrimaryPy_);
        PrimaryPz->push_back(PrimaryPz_);
      }
    }
    else if(TreeName_=="Background") 
    {
      for(int i = 0; i < input_tree->GetEntries(); i++)
      {
        input_tree->GetEntry(i);
        BackgroundDecayTime->push_back(BackgroundDecayTime_);
        BackgroundAtomicNumber->push_back(BackgroundAtomicNumber_);
        BackgroundAtomicMass->push_back(BackgroundAtomicMass_);
        BackgroundVertexX->push_back(BackgroundVertexX_);
        BackgroundVertexY->push_back(BackgroundVertexY_);
        BackgroundVertexZ->push_back(BackgroundVertexZ_);
      }
    }
}//Get_Event

//--------------------------------------------------------------------------
// Reset the event
 void ROOTFileManager::EventReset()
 {
    int run=-1;
    hitX_start->clear(); 
    hitX_end->clear();
    hitY_start->clear();
    hitY_end->clear();
    hitZ_start->clear();
    hitZ_end->clear();
    time_start->clear();
    time_end->clear();
    edep->clear();
    length->clear();
    pdg->clear();

    InitialPDG.clear();
    InteractionTime->clear();
    InitialEnergy->clear();
    InteractionVertexX->clear();
    InteractionVertexY->clear();
    InteractionVertexZ->clear();

    PrimaryPDG.clear();
    PrimaryEnergy->clear();
    PrimaryPx->clear();
    PrimaryPy->clear();
    PrimaryPz->clear();

    BackgroundDecayTime->clear();
    BackgroundAtomicMass->clear();
    BackgroundAtomicNumber->clear();
    BackgroundVertexX->clear();
    BackgroundVertexY->clear();
    BackgroundVertexZ->clear();
 }
 //--------------------------------------------------------------------------
// Get number of entries
 int ROOTFileManager::NEntries()
 {
    return input_tree->GetEntries();
 }

//--------------------------------------------------------------------------
// Close input file
 void ROOTFileManager::CloseInput()
 {
    input_file->Close();
    delete input_file;
 }


