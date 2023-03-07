/* main.cc */

#include <iostream>

#include "SemiAnalyticalModel.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>

//
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "SimPhotons.h"
#include "PropagationTimeModel.h"
#include "EnergyDeposition.h"


#include <nlohmann/json.hpp>

using json = nlohmann::json;

using namespace std;

int main(int argc, char **argv)
{
  
 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  
  char* fileName=0;
  for (int i = 1; i < argc; i++)
  {  
   if (i + 1 != argc)
   {
     if (strcmp(argv[i], "-file") == 0) 
      {                 
       fileName = argv[i + 1];   
       std::cout << "fname " << fileName << std::endl; // Curioso que std::cout trata a char * fileName como un puntero a un string de C, 
                                                        // y hace el cout como si fuera string y no puntero.
        i++;   
      }
    }
  }
   
  if(fileName==0) 
  {
   std::cout << "No file name given, aborting execution" << std::endl; 
   //abort();
  }






  std::ifstream f("params.json");
  json OpParams = json::parse(f);


  //Input file:
  TFile input_file("banger_100MeVKE.root");
  TTree *input_tree = (TTree*)input_file.Get("event_tree;4");
  
  // Reading root file information. Some cleanup needed.
  int run;
  std::vector<double> *hitX_start = new std::vector<double>();
  std::vector<double> *hitX_end = new std::vector<double>();
  std::vector<double> *hitY_start = new std::vector<double>();
  std::vector<double> *hitY_end = new std::vector<double>();
  std::vector<double> *hitZ_start = new std::vector<double>();
  std::vector<double> *hitZ_end = new std::vector<double>();
  std::vector<double> *time_start = new std::vector<double>();
  std::vector<double> *time_end = new std::vector<double>();
  std::vector<double> *edep = new std::vector<double>();
  std::vector<double> *length = new std::vector<double>();

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

  //Output file:
  TFile *OutputFile = TFile::Open("Output_19.root", "RECREATE");

  std::unique_ptr<std::vector<SimPhotons>> photonCol{new std::vector<SimPhotons>{}};
  auto& photonHitCollection{*photonCol};
  unsigned int fNOpChannels = OpParams["nOpDet"];
  photonHitCollection.resize(fNOpChannels);
  for (size_t i = 0; i < fNOpChannels; ++i)
  {
    photonHitCollection[i].OpChannel = i;
  }
  std::vector<int> DetectedNum(fNOpChannels);
  std::vector<double> OpDetVisibilities;
  std::unique_ptr<SemiAnalyticalModel> semi;
  semi = std::make_unique<SemiAnalyticalModel>(OpParams); //Initialize SemiAnalyticalModel object
  std::unique_ptr<PropagationTimeModel> PropTime;
  PropTime = std::make_unique<PropagationTimeModel>(OpParams); //Initialize PropagationTimeModel object
  // Fill the OpDetVisibilities vector with the give scintillation point
  double cum_edep=0;
  int generated_counter=0;
  unsigned long runID;
  std::vector<std::vector<double>> SimPhotons;

  for (size_t nRun = 19; nRun < 20; nRun++)
  {    
    int time_call=0;
    SimPhotons.clear();
    SimPhotons.resize(fNOpChannels);
    TTree *OutputTree = new TTree("event", "event");
    std::cout << "Reading event number: " << nRun << std::endl;
    input_tree->GetEntry(nRun);
    unsigned int nPhotons=0;
    runID=nRun;
    for (size_t nHit = 0; nHit < hitX_start->size(); nHit++ )
     {
     // Initialize the energy deposition object with its StartPoint and the EndPoint:
     SemiAnalyticalModel::Point_t StartPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
     SemiAnalyticalModel::Point_t EndPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
     std::unique_ptr<EnergyDeposition> Edep;
     Edep = std::make_unique<EnergyDeposition>(OpParams, edep->at(nHit), StartPoint, EndPoint, time_start->at(nHit), time_end->at(nHit) ,length->at(nHit));
     SemiAnalyticalModel::Point_t ScintPoint{(Edep->MidPoint()).x, (Edep->MidPoint()).y, (Edep->MidPoint()).z};
     semi->detectedDirectVisibilities(OpDetVisibilities, ScintPoint);
     double nphot_fast=Edep->LArQL(); 

    // std::cout << "LightYield " <<  nphot_fast/Edep->Energy() << std::endl; 
     //double nphot_fast=Edep->Energy()*24000; 
     generated_counter = generated_counter + nphot_fast;
     // Fill the DetectedNum vector with the given OpDetVisibilities
     cum_edep=cum_edep+Edep->Energy();
    // std::cout << "Edep cum: " << cum_edep << std::endl; 
    // std::cout << "Generated total: " << generated_counter << std::endl; 
     semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot_fast);
     for (int channel = 0 ; channel< DetectedNum.size() ; channel++)
      {
       int n_detected = DetectedNum.at(channel);
       nPhotons=nPhotons+n_detected;
       std::vector<double> transport_time;
       transport_time.resize(n_detected);
       PropTime->propagationTime(transport_time, ScintPoint, channel);
       for (size_t i = 0; i < n_detected; ++i)
        {
         int time=0; 
         time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+PropTime->ScintTime() );
         ++photonHitCollection[channel].DetectedPhotons[time];
         time_call+=1;
        }
      }
    // std::cout << "Detected Photons "<<  nPhotons << std::endl; 
   }// end hits loop    
    // With the Photon
    for ( auto & fPhotons : (photonHitCollection) )
    {
      int opChannel = fPhotons.OpChannel;
      std::map<int, int> fPhotons_map = fPhotons.DetectedPhotons;
      for (auto fPhotons = fPhotons_map.begin(); fPhotons!= fPhotons_map.end(); fPhotons++)
      {       
         for(int i = 0; i < fPhotons->second ; i++)
         {
         SimPhotons.at(opChannel).push_back(fPhotons->first);
         }
      }
    }
    std::cout << "The number of time calls is : " << time_call << std::endl;
    OutputTree->Branch("eventID", &runID);
    OutputTree->Branch("SimPhotonsperOpChVUV",&SimPhotons);
    OutputTree->Branch("GeneratedPhotons",&generated_counter);
    OutputTree->Branch("DetectedPhotons",&nPhotons);
    OutputTree->Fill();
    OutputFile->cd();
  }// end event loop
    OutputFile->Write();
    OutputFile->Close();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000 << "[s]" << std::endl;
  return 0;
}


