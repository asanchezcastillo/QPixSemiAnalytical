// -----------------------------------------------------------------------------
//  ROOTFileManager.cpp
//
//  Class definition of the ROOT file manager
//   * Author: Everybody is an author!
//   * Creation date: 31 August 2020
// -----------------------------------------------------------------------------

#include "ROOTFileManager.h"

// ROOT includes
#include "TObject.h"

// math includes
#include <math.h>

//--------------------------------------------------------------------------
ROOTFileManager::ROOTFileManager(char* inputFileName, char* outputFileName)
{
    Initialize(inputFileName, outputFileName);
}

//--------------------------------------------------------------------------
ROOTFileManager::~ROOTFileManager()
{}

//--------------------------------------------------------------------------

void ROOTFileManager::Initialize(char* inputFileName, char* outputFileName)
{
    //Input file:
    input_file = new TFile(inputFileName);
    input_tree = (TTree*)input_file->Get("event");
    //TFile *OutputFile = TFile::Open(outputFileName, "RECREATE");  /TO BE DONE: USE RFM TO FILL OUTPUT FILE
    set_branch_addresses(input_tree);
}

//--------------------------------------------------------------------------
void ROOTFileManager::set_branch_addresses(TTree *input_tree)
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
}

void ROOTFileManager::EventFill()
{
    OutputTree->Branch("eventID", &runID);
    OutputTree->Branch("SimPhotonsperOpChVUV",&SavePhotons);
    OutputTree->Branch("GeneratedPhotons",&generated_counter);
    OutputTree->Branch("DetectedPhotons",&nPhotons);
    OutputTree->Fill();
    OutputFile->cd();
}

//--------------------------------------------------------------------------
void ROOTFileManager::Save()
{
    OutputFile->Write();
    OutputFile->Close();
}

//--------------------------------------------------------------------------
// gets the event from the file and tunrs it into electrons
void ROOTFileManager::GetEvent(int nRun)
{
    input_tree->GetEntry(nRun);
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



