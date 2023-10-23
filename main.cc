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
    string nameroot = names.at(n);
    nameroot.resize(nameroot.size()-5);
    string outputname = nameroot+"_SemiAnalytical.root";
    const char* outputFile = outputname.c_str();
    std::unique_ptr<ROOTFileManager> rfm;
    std::cout << "Input file name: " << inputFile << " Output file name " << outputFile << std::endl;

    // Get information on the inital particle

    rfm = std::make_unique<ROOTFileManager>((char*)inputFile, (char*)outputFile, "InitialParticle"); //Initialize rfm object with the intial particle information
    rfm->GetEvent();
    rfm->CloseInput();
    std::vector<double> * InitialParticleEnergy = rfm->GetInitialEnergy(); 
    std::vector<double> * InteractionTime = rfm->GetInteractionTime(); 
    std::vector<int> InitialParticlePDG = rfm->GetInitialPDG();
    std::vector<double> * InteractionVertexX = rfm->GetInteractionVertexX();
    std::vector<double> * InteractionVertexY = rfm->GetInteractionVertexY();
    std::vector<double> * InteractionVertexZ = rfm->GetInteractionVertexZ();
    SemiAnalyticalModel::Point_t InteractionVertex;
    if(InteractionVertexX->size()>0)
    {
      InteractionVertex={InteractionVertexX->at(0), InteractionVertexY->at(0), InteractionVertexZ->at(0)};      
    }
    else
    {
      InteractionVertex={0.,0.,0.};
    }

    // Get information on the primary particles

    rfm = std::make_unique<ROOTFileManager>((char*)inputFile, (char*)outputFile, "PrimaryParticle"); //Initialize rfm object with the primary particle information
    rfm->GetEvent();
    rfm->CloseInput();
    std::vector<double> * PrimaryParticleEnergy = rfm->GetPrimaryEnergy(); 
    std::vector<int> PrimaryParticlePDG = rfm->GetPrimaryPDG();
    std::vector<double> * PrimaryParticlePx = rfm->GetPrimaryPx();
    std::vector<double> * PrimaryParticlePy = rfm->GetPrimaryPy();
    std::vector<double> * PrimaryParticlePz = rfm->GetPrimaryPz();

    // Get information on the background
    rfm = std::make_unique<ROOTFileManager>((char*)inputFile, (char*)outputFile, "Background"); //Initialize rfm object with the background information
    rfm->GetEvent();
    rfm->CloseInput();
    std::vector<double> * BackgroundDecayTime = rfm->GetBackgroundDecayTime(); 
    std::vector<int> * BackgroundAtomicNumber = rfm->GetBackgroundAtomicNumber();
    std::vector<int> * BackgroundAtomicMass = rfm->GetBackgroundAtomicMass();
    std::vector<double> * BackgroundVertexX = rfm->GetBackgroundVertexX(); 
    std::vector<double> * BackgroundVertexY = rfm->GetBackgroundVertexY(); 
    std::vector<double> * BackgroundVertexZ = rfm->GetBackgroundVertexZ(); 

    // Get information on the input to the event

    rfm = std::make_unique<ROOTFileManager>((char*)inputFile, (char*)outputFile, "Hits"); //Initialize rfm object with hits information

    if(nMaxEvents==0)
    {
      nMaxEvents = rfm->NEntries(); //Max number of events is the number of entries in our input tree. Each entry contains a vector with all the information on the energy depositions. 
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
      std::vector<double> * LightYield = new std::vector<double>();
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
      std::vector<int> *pdg = rfm->GetPDG();
      rfm->CloseInput();
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
      
      std::vector<int> * num_photons = new std::vector<int>();
      std::vector<int> * num_electrons = new std::vector<int>();
      std::vector<double> * energy_deposited = new std::vector<double>();

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
      TTree *PhotonsTree = new TTree("Photons", "Photons"); // Tree containing photons information
      TTree *EventTree = new TTree("Event", "Event"); // Tree containing information on the event 
      TTree *GeometryTree = new TTree("Geometry", "Geometry"); // Tree containing information on the geoemtry of the event
      TTree *BackgroundTree = new TTree("Brackground", "Brackground"); // Tree containing information on the geoemtry of the event

      std::cout << "Reading event number: " << nRun << std::endl;
      std::vector<std::vector<double>> *distance = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *angle = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *photons_per_edep = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *channels = new std::vector<std::vector<double>>();
      std::vector<std::vector<double>> *visibility_vector = new std::vector<std::vector<double>>();
      std::vector<double> visibility_at_vertex;

      std::vector<std::vector<double>> *BackgroundVisibility = new std::vector<std::vector<double>>();
      BackgroundVisibility->resize(BackgroundVertexX->size());

      distance->resize(fNOpChannels);
      angle->resize(fNOpChannels);
      photons_per_edep->resize(fNOpChannels);
      channels->resize(fNOpChannels);
      visibility_vector->resize(hitX_start->size());
      semi->detectedDirectVisibilities(visibility_at_vertex, InteractionVertex);
      for(size_t i = 0; i < BackgroundVertexX->size(); i++)
      {
        SemiAnalyticalModel::Point_t BackgroundVertex={BackgroundVertexX->at(i), BackgroundVertexY->at(i), BackgroundVertexZ->at(i)};
        semi->detectedDirectVisibilities(BackgroundVisibility->at(i), BackgroundVertex);
      }
      for (size_t nHit = 0; nHit < hitX_start->size(); nHit++)
      {
      // Initialize the energy deposition object with its StartPoint and the EndPoint:
      SemiAnalyticalModel::Point_t StartPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      SemiAnalyticalModel::Point_t EndPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      std::unique_ptr<EnergyDeposition> Edep;
      Edep = std::make_unique<EnergyDeposition>(OpParams, edep->at(nHit), StartPoint, EndPoint, time_start->at(nHit), time_end->at(nHit) ,length->at(nHit), pdg->at(nHit));
      SemiAnalyticalModel::Point_t ScintPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
      semi->detectedDirectVisibilities(OpDetVisibilities, ScintPoint);
      //const int nphot= round(Edep->Energy()*24000); // Number of photons computed from a constant LY. 
      Edep->LArQL();
      const int nphot = round(Edep->GetNumberPhotons());
      num_photons->push_back(nphot);
      num_electrons->push_back(round(Edep->GetNumberElectrons()));      
      energy_deposited->push_back(Edep->Energy());
      LightYield->push_back(Edep->GetLightYield());
      generated_counter = generated_counter+nphot;
      cum_edep+=Edep->Energy();
      semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot);
      for(int channel =0; channel< fNOpChannels; channel++)
      {
        visibility_vector->at(nHit).push_back(OpDetVisibilities.at(channel));
      }
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
          double ScintTime = PropTime->ScintTime(pdg->at(nHit));
          time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+ ScintTime );
          ++photonHitCollection[channel].DetectedPhotons[time];
          }
        }// end channels loop
      }// end hits loop 
      for ( auto & fPhotons : (photonHitCollection) )
      {
        int opChannel = fPhotons.OpChannel;
        std::map<int, int> fPhotons_map = fPhotons.DetectedPhotons;
        for (auto fPhotons = fPhotons_map.begin(); fPhotons!= fPhotons_map.end(); fPhotons++){       
          for(int i = 0; i < fPhotons->second ; i++)
          {
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
      std::vector<double> *visibility_average = new std::vector<double>();
      for (int channel=0; channel<fNOpChannels; channel++)
      {
        double weighted_visibility=0;
        double weight_visibility=0;
        double sum = 0;
        for (int j = 0; j<edep->size(); j++ )
        {
          weight_visibility = weight_visibility + visibility_vector->at(j).at(channel) * edep->at(j);
          sum = sum + edep->at(j);
        }
        weighted_visibility=weight_visibility/sum;
        visibility_average->push_back(weighted_visibility);
      }

      //Visibility from vertex. 

      std::cout << "The sum of all the visbilities is: " << accumulate(visibility_average->begin(),visibility_average->end(),0.000);
      PhotonsTree->Branch("eventID", &runID);
      PhotonsTree->Branch("SavedPhotons",&SavePhotons);
      PhotonsTree->Branch("GeneratedPhotons",&generated_counter);
      PhotonsTree->Branch("DetectedPhotons",&nPhotons);
      PhotonsTree->Branch("PhotonsPerEdep", &photons_per_edep);
      GeometryTree->Branch("Distance", &distance);
      GeometryTree->Branch("Angle", &angle);
      GeometryTree->Branch("DistanceAverage", &distance_average);
      GeometryTree->Branch("AngleAverage", &angle_average);
      GeometryTree->Branch("Channels", &channels);
      GeometryTree->Branch("VisibilityVector", &visibility_average);
      GeometryTree->Branch("VisibilityAtVertex", &visibility_at_vertex);
      EventTree->Branch("EventMeanX", &event_x);
      EventTree->Branch("EventMeanY", &event_y);
      EventTree->Branch("EventMeanZ", &event_z);
      EventTree->Branch("InitialParticleEnergy", &InitialParticleEnergy);
      EventTree->Branch("TotalEdep", &cum_edep);
      EventTree->Branch("InitialParticlePDG", &InitialParticlePDG);
      EventTree->Branch("InteractionTime", &InteractionTime);
      EventTree->Branch("PrimaryParticleEnergy", &PrimaryParticleEnergy);
      EventTree->Branch("PrimaryParticlePDG", &PrimaryParticlePDG);
      EventTree->Branch("PrimaryParticlePx", &PrimaryParticlePx);
      EventTree->Branch("PrimaryParticlePy", &PrimaryParticlePy);
      EventTree->Branch("PrimaryParticlePz", &PrimaryParticlePz);
      EventTree->Branch("InteractionVertexX", &InteractionVertexX);
      EventTree->Branch("InteractionVertexY", &InteractionVertexY);
      EventTree->Branch("InteractionVertexZ", &InteractionVertexZ);
      EventTree->Branch("LightYield", &LightYield);
      EventTree->Branch("NumberOfPhotons", &num_photons);
      EventTree->Branch("NumberOfElectrons", &num_electrons);
      EventTree->Branch("DepositedEnergy", &energy_deposited);
      BackgroundTree->Branch("BackgroundDecayTime", &BackgroundDecayTime);
      BackgroundTree->Branch("BackgroundAtomicNumber", &BackgroundAtomicNumber);
      BackgroundTree->Branch("BackgroundAtomicMass", &BackgroundAtomicMass);
      BackgroundTree->Branch("BackgroundVertexX", &BackgroundVertexX);
      BackgroundTree->Branch("BackgroundVertexY", &BackgroundVertexY);
      BackgroundTree->Branch("BackgroundVertexZ", &BackgroundVertexZ);
      BackgroundTree->Branch("BackgroundVisibility", &BackgroundVisibility);
      PhotonsTree->Fill();
      GeometryTree->Fill();
      EventTree->Fill();
      BackgroundTree->Fill();
      OutputFile->cd();
      rfm->EventReset();
      photonHitCollection.clear();
      SavePhotons.clear();
      distance->clear();
      angle->clear();
      channels->clear();
      photons_per_edep->clear();
      num_photons->clear();
      num_electrons->clear();
      energy_deposited->clear();
      InteractionVertexX->clear();
      InteractionVertexY->clear();
      InteractionVertexZ->clear();
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      std::cout << "Event: " << nRun <<" Elapsed time: " <<  std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000 << "[s]" << std::endl;
    }// end event loop
    OutputFile->Write();
    delete OutputFile;
  }
  return 0;
}
