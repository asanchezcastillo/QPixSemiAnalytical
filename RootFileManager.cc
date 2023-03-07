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
ROOTFileManager::ROOTFileManager(std::string const& input_file, std::string const& output_file)
{
    Initialize(input_file, output_file);
}

//--------------------------------------------------------------------------
ROOTFileManager::~ROOTFileManager()
{}

//--------------------------------------------------------------------------


void ROOTFileManager::Initialize(std::string const& input_file, std::string const& output_file)
{
    
    //Input file:
    TFile input_file("banger_100MeVKE.root");
    TTree *input_tree = (TTree*)input_file.Get("event_tree;4");

    TFile *OutputFile = TFile::Open("Output.root", "RECREATE");

    set_branch_addresses();
    
}


//--------------------------------------------------------------------------
void ROOTFileManager::set_branch_addresses()
{
    run=-1;
    hitX_start = new std::vector<double>();
    hitX_end = new std::vector<double>();
    hitY_start = new std::vector<double>();
    hitY_end = new std::vector<double>();
    hitZ_start = new std::vector<double>();
    hitZ_end = new std::vector<double>();
    time_start = new std::vector<double>();
    time_end = new std::vector<double>();
    edep = new std::vector<double>();
    length = new std::vector<double>();

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

/*    
//--------------------------------------------------------------------------
unsigned int ROOTFileManager::NumberEntries()
{
    if (ttree_) return ttree_->GetEntries();
    return -1;
}

//--------------------------------------------------------------------------
*/
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
void ROOTFileManager::EventReset(int fNOpChannels)
{   
    runID=-1;
    generated_counter=0;
    nPhotons=0;
    SavePhotons = new std::vector<std::vector<double>>();
    SavePhotons.resize(fNOpChannels);
    TTree *OutputTree = new TTree("event", "event");

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
    input_tree->input_tree->GetEntry(nRun);;

}//Get_Event




