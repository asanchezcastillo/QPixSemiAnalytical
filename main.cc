// SemiAnalyticalModel
//  - fast optical simulation of scintillation photons using semi-analytical model.

// March 2022 by A. Sánchez Castillo

#include <iostream>



#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <nlohmann/json.hpp>

//Root includes
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

//Class headers includes
#include "SemiAnalyticalModel.h"
#include "SimPhotons.h"
#include "PropagationTimeModel.h"
#include "EnergyDeposition.h"
#include"RootFileManager.h"
#include"RootFileManager.cc"



using json = nlohmann::json;

using namespace std;

int main(int argc, char **argv)
{
  
  char* inputFileName=0;
  char* outputFileName=0;
  for (int i = 1; i < argc; i++)
  {  
   if (i + 1 != argc)
   {
    if (strcmp(argv[i], "--file") == 0 || strcmp(argv[i], "-f") == 0 )  
    {                 
     inputFileName = argv[i + 1];   
     i++;   
    }
    if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0 )  
    {                 
     outputFileName = argv[i + 1];   
     i++;   
    }
   }
  }
   
  if(inputFileName==0) 
  {
   std::cout << "No input file name given, aborting execution" << std::endl; 
   abort();
  }
  if(outputFileName==0) 
  {
   std::cout << "No output file name given, aborting execution" << std::endl; 
   abort();
  }
  std::ifstream f("params.json");
  json OpParams = json::parse(f);

  std::unique_ptr<ROOTFileManager> rfm;
  rfm = std::make_unique<ROOTFileManager>(inputFileName, outputFileName); //Initialize rfm object
  int nEntries = rfm->NEntries();
   std::cout << "The number of events is: " << rfm->GetEntries() << std::endl;
  //Output file:
  TFile *OutputFile = TFile::Open(outputFileName, "RECREATE");
  // Initialize PhotonHitCollection to store simulated hits.
  std::unique_ptr<std::vector<SimPhotons>> photonCol{new std::vector<SimPhotons>{}};
  auto& photonHitCollection{*photonCol};
  unsigned int fNOpChannels = OpParams["nOpDet"];
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
  std::vector<std::vector<double>> SavePhotons;
  
  for (size_t nRun = 0; nRun < 20; nRun++)
  {    
    photonHitCollection.resize(fNOpChannels);
    for (size_t i = 0; i < fNOpChannels; ++i)
    {
      photonHitCollection[i].OpChannel = i;
    }
    SavePhotons.resize(fNOpChannels);
    rfm->GetEvent(nRun);
    // Reading root file information.
    std::vector<double> *hitX_start = rfm->GetXStart();
    std::vector<double> *hitX_end = rfm->GetXEnd();
    std::vector<double> *hitY_start = rfm->GetYStart();
    std::vector<double> *hitY_end = rfm->GetYEnd();
    std::vector<double> *hitZ_start = rfm->GetZStart();
    std::vector<double> *hitZ_end = rfm->GetZEnd();
    std::vector<double> *time_start = rfm->GetTimeStart();
    std::vector<double> *time_end = rfm->GetTimeEnd();
    std::vector<double> *edep = rfm->GetEdep();
    std::vector<double> *length = rfm->GetLength();

    //Start counting time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    TTree *OutputTree = new TTree("event", "event");
    std::cout << "Reading event number: " << nRun << std::endl;
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
     double nphot=Edep->LArQL(); 
     //double nphot=Edep->Energy()*24000; 
     generated_counter = generated_counter + nphot;
     // Fill the DetectedNum vector with the given OpDetVisibilities
     cum_edep=cum_edep+Edep->Energy(); //Cumulative energy deposition per step. Currently not stored.
    //  std::cout << "Cum edep: " << cum_edep << " nHit " << nHit <<std::endl;
    //  std::cout << "Number of photons: "<< nphot << std::endl;
     semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot);
     for (int channel = 0 ; channel< DetectedNum.size() ; channel++)
      {
       int n_detected = DetectedNum.at(channel);
       nPhotons+=n_detected;
       std::vector<double> transport_time;
       transport_time.resize(n_detected);
       PropTime->propagationTime(transport_time, ScintPoint, channel);
       for (size_t i = 0; i < n_detected; ++i)
        {
         int time=0; 
         time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+PropTime->ScintTime() );
         ++photonHitCollection[channel].DetectedPhotons[time];
        }
      }
   }// end hits loop 
    for ( auto & fPhotons : (photonHitCollection) )
    {
     int opChannel = fPhotons.OpChannel;
     std::map<int, int> fPhotons_map = fPhotons.DetectedPhotons;
     for (auto fPhotons = fPhotons_map.begin(); fPhotons!= fPhotons_map.end(); fPhotons++){       
       for(int i = 0; i < fPhotons->second ; i++){
       SavePhotons.at(opChannel).push_back(fPhotons->first);
       }
     }
    }
    OutputTree->Branch("eventID", &runID);
    OutputTree->Branch("SavedPhotons",&SavePhotons);
    OutputTree->Branch("GeneratedPhotons",&generated_counter);
    OutputTree->Branch("DetectedPhotons",&nPhotons);
    OutputTree->Fill();
    OutputFile->cd();
    rfm->EventReset();
    photonHitCollection.clear();
    SavePhotons.clear();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Event: " << nRun <<" Elapsed time: " <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000 << "[s]" << std::endl;
  }// end event loop
  OutputFile->Write();
  OutputFile->Close();
  return 0;
}


