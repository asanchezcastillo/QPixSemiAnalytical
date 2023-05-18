// SemiAnalyticalModel
//  - fast optical simulation of scintillation photons using semi-analytical model.
// March 2022 by A. Sánchez Castillo 
// asanchezcastillo@ugr.es

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
  int nMaxEvents=0;
  for (int i = 1; i < argc; i++)
  {  
   if (i + 1 != argc)
   {
    if (strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0 )  //Input file name
    {                 
     inputFileName = argv[i + 1];   
     i++;   
    }
    if (strcmp(argv[i], "-n") == 0 )  //Max number of events
    {                 
     nMaxEvents = atoi(argv[i + 1]);   
     i++;   
    }
   }
  }
  if(inputFileName==0) 
  {
   std::cout << "WARNING: No input file name given, aborting execution" << std::endl; 
   abort();
  }
  //Start file loop
  ifstream Traks_file2(inputFileName);
  if(!Traks_file2) cerr << "WARNING:  Failed to open file with Input file names"<< endl;
  Traks_file2.seekg(0);
  vector<string> names;
  string nombre;
  while(!(Traks_file2.eof())) {
    Traks_file2 >> nombre;
    names.push_back(nombre);
  }
  const int n_files = names.size() - 1;
  cout<< "Number of input files: "<<n_files<<endl;
  std::ifstream f("params.json");
  json OpParams = json::parse(f);

  std::unique_ptr<PropagationTimeModel> PropTime;
  PropTime = std::make_unique<PropagationTimeModel>(OpParams); //Initialize PropagationTimeModel object
  std::unique_ptr<SemiAnalyticalModel> semi;
  semi = std::make_unique<SemiAnalyticalModel>(OpParams); //Initialize SemiAnalyticalModel object

  for(int n=0; n<n_files; n++)
  {  
    const char* inputFile=0;
    inputFile = names.at(n).c_str();

    string outputname = "output_"+std::to_string(n)+".root";
    const char* outputFile = outputname.c_str();;

    std::unique_ptr<ROOTFileManager> rfm;
    std::cout << "Input file name: " << inputFile << " Output file name " << outputFile << std::endl;
    rfm = std::make_unique<ROOTFileManager>((char*)inputFile, (char*)outputFile); //Initialize rfm object

    if(nMaxEvents==0)
    {
      nMaxEvents = rfm->NEntries();
    }
    else
    {
      if(nMaxEvents>rfm->NEntries())
      {
        std::cout << "WARNING: The input file does not contain that many events! Aborting execution." << std::endl;
        abort();
      }
    }
    std::cout << "The number of events to be simulated is: " << nMaxEvents << std::endl;
    //Output file:

    TFile *OutputFile = TFile::Open((char*)outputFile, "RECREATE");
    // Initialize PhotonHitCollection to store simulated hits.
    std::unique_ptr<std::vector<SimPhotons>> photonCol{new std::vector<SimPhotons>{}};
    auto& photonHitCollection{*photonCol};
    unsigned int fNOpChannels = OpParams["nOpDet"];
    std::vector<int> DetectedNum(fNOpChannels);
    std::vector<double> OpDetVisibilities;
    unsigned long runID;
    std::vector<std::vector<double>> SavePhotons;
    
    for (size_t nRun = 0; nRun < nMaxEvents; nRun++)
    {    
      double cum_edep=0;
      unsigned int generated_counter=0;
      unsigned int nPhotons=0;
      photonHitCollection.resize(fNOpChannels);
      runID=nRun;
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

      //Compute event information (weighted drift distance)
      double event_x=0;
      double event_y=0;
      double event_z=0;

      double weighted_x=0;
      double weight_x=0;
      double weighted_y=0;
      double weight_y=0;
      double weighted_z=0;
      double weight_z=0;
      double sum = 0;
      double weight=0;
        for (int i=0; i<hitX_start->size(); i++)
        {
          if(edep->at(i)!= 0)
          {
            weight_x = weight_x+ 0.5*(hitX_start->at(i)+hitX_end->at(i)) * edep->at(i);
            weight_y = weight_y+ 0.5*(hitY_start->at(i)+hitY_end->at(i)) * edep->at(i);
            weight_z = weight_z+ 0.5*(hitZ_start->at(i)+hitZ_end->at(i)) * edep->at(i);
            sum = sum + edep->at(i);
          }
        }
      weighted_x = weight_x/sum;
      weighted_y = weight_y/sum;  
      weighted_z = weight_z/sum; 
      event_x=weighted_x;
      event_y=weighted_y;
      event_z=weighted_z; 
  
      //Start counting time
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      TTree *OutputTree = new TTree("event", "event");
      std::cout << "Reading event number: " << nRun << std::endl;
      std::vector<std::vector<double>> *distance = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *angle = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *photons_per_edep = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *channels = new std::vector<std::vector<double>>();

      distance->resize(fNOpChannels);
      angle->resize(fNOpChannels);
      photons_per_edep->resize(fNOpChannels);
      channels->resize(fNOpChannels);

      for (size_t nHit = 0; nHit < hitX_start->size(); nHit++ )
      {
      // Initialize the energy deposition object with its StartPoint and the EndPoint:
      SemiAnalyticalModel::Point_t StartPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      SemiAnalyticalModel::Point_t EndPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      std::unique_ptr<EnergyDeposition> Edep;
      Edep = std::make_unique<EnergyDeposition>(OpParams, edep->at(nHit), StartPoint, EndPoint, time_start->at(nHit), time_end->at(nHit) ,length->at(nHit));
      SemiAnalyticalModel::Point_t ScintPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      semi->detectedDirectVisibilities(OpDetVisibilities, ScintPoint);
      double nphot=Edep->Energy()*24000; // Number of photons computed from a constant LY. To be replaced with LArQL. 
      generated_counter = generated_counter + nphot;
      cum_edep+=Edep->Energy();
      double sum_of_elems=0;
      for(int i=0; i<OpDetVisibilities.size(); i++)
      sum_of_elems += OpDetVisibilities[i];
      semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot);
      
      for (int channel = 0 ; channel< DetectedNum.size() ; channel++)
        {
        int n_detected = DetectedNum.at(channel);
        if(n_detected>0)
        {
        distance->at(channel).push_back(semi->GetDistance(ScintPoint, channel));
        angle->at(channel).push_back(semi->GetAngle(ScintPoint, channel));
        photons_per_edep->at(channel).push_back(DetectedNum.at(channel));
        channels->at(channel).push_back(static_cast<double>(channel));
        }
        nPhotons+=n_detected;
        std::vector<double> transport_time;
        transport_time.resize(n_detected);
        PropTime->propagationTime(transport_time, ScintPoint, channel);
        for (size_t i = 0; i < n_detected; ++i)
          {
          int time=0; 
          time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+ PropTime->ScintTime() );
          ++photonHitCollection[channel].DetectedPhotons[time];
          }
        }// end channels loop
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
    //Loop to compute the weighted mean for each optical channel
    std::vector<double> *distance_average = new std::vector<double>();
    std::vector<double> *angle_average = new std::vector<double>();
    for (int channel=0; channel<fNOpChannels; channel++) 
    {
      double weighted_distance=0;
      double weight_distance=0;
      double weighted_angle=0;
      double weight_angle=0;
      double sum = 0;
      double weight=0;
        for (int j=0; j<distance->at(channel).size(); j++)
        {
            weight_distance = weight_distance + distance->at(channel).at(j) * photons_per_edep->at(channel).at(j);
            weight_angle = weight_angle + angle->at(channel).at(j) * photons_per_edep->at(channel).at(j);
            sum = sum + photons_per_edep->at(channel).at(j);
        }
      weighted_distance = weight_distance/sum;
      weighted_angle = weight_angle/sum;  
      distance_average->push_back(weighted_distance);
      angle_average->push_back(weighted_angle);
    }
    OutputTree->Branch("eventID", &runID);
    OutputTree->Branch("SavedPhotons",&SavePhotons);
    OutputTree->Branch("GeneratedPhotons",&generated_counter);
    OutputTree->Branch("DetectedPhotons",&nPhotons);
    OutputTree->Branch("TotalEdep", &cum_edep);
    OutputTree->Branch("Distance", &distance);
    OutputTree->Branch("Angle", &angle);
    OutputTree->Branch("Distance_average", &distance_average);
    OutputTree->Branch("Angle_average", &angle_average);
    OutputTree->Branch("Channels", &channels);
    OutputTree->Branch("photons_per_edep", &photons_per_edep);
    OutputTree->Branch("event_x", &event_x);
    OutputTree->Branch("event_y", &event_y);
    OutputTree->Branch("event_z", &event_z);
    OutputTree->Fill();
    OutputFile->cd();
    rfm->EventReset();
    photonHitCollection.clear();
    SavePhotons.clear();
    distance->clear();
    angle->clear();
    channels->clear();
    photons_per_edep->clear();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Event: " << nRun <<" Elapsed time: " <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000 << "[s]" << std::endl;
    }// end event loop
    OutputFile->Write();
    OutputFile->Close();
  }
  return 0;
}
