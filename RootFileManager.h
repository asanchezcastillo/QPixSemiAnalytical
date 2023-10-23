// -----------------------------------------------------------------------------
//  ROOTFileManager.h
//
//  Class definition of the ROOT file manager
//   * Author: Everybody is an author!
//   * Creation date: 31 August 2020
// -----------------------------------------------------------------------------

#ifndef ROOTFileManager_h
#define ROOTFileManager_h 

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

// Include for Point_t
#include "SemiAnalyticalModel.h"

// C++ includes
#include <map>
#include <set>

class ROOTFileManager {

    public:
        ROOTFileManager(char* inputFileName, char* outputFileName, std::string TreeName);
        ~ROOTFileManager();
        
        unsigned int NumberEntries();
        void EventFill();
        void Save();
        void GetEvent(int);
        void GetEvent();
        void EventReset();
        int NEntries();
        void CloseInput();
        int GetRun(){return run;}
        std::vector<double> * GetInteractionVertexX(){return InteractionVertexX;}
        std::vector<double> * GetInteractionVertexY(){return InteractionVertexY;}
        std::vector<double> * GetInteractionVertexZ(){return InteractionVertexZ;}
        std::vector<double> * GetXStart(){return hitX_start;}
        std::vector<double> * GetXEnd(){return hitX_end;}
        std::vector<double> * GetYStart(){return hitY_start;}
        std::vector<double> * GetYEnd(){return hitY_end;}
        std::vector<double> * GetZStart(){return hitZ_start;}
        std::vector<double> * GetZEnd(){return hitZ_end;}
        std::vector<double> * GetTimeStart(){return time_start;}
        std::vector<double> * GetTimeEnd(){return time_end;}
        std::vector<double> * GetEdep(){return edep;}
        std::vector<double> * GetLength(){return length;}
        std::vector<int> * GetPDG(){return pdg;}
        std::vector<int> GetInitialPDG(){return InitialPDG;} // Vector containing pdgs of all intial particles
        std::vector<int> GetPrimaryPDG(){return PrimaryPDG;} // Vector containing pdgs of all primaries particles
        std::vector<double> * GetInteractionTime(){return InteractionTime;} // Vector containing interaction times of all intial particles
        std::vector<double> * GetInitialEnergy(){return InitialEnergy;} // Vector containing energies of all intial particles
        std::vector<double> * GetPrimaryEnergy(){return PrimaryEnergy;} // Vector containing energies of all primary particles
        std::vector<double> * GetPrimaryPx(){return PrimaryPx;} // Vector containing px of all primary particles
        std::vector<double> * GetPrimaryPy(){return PrimaryPy;} // Vector containing py of all primary particles
        std::vector<double> * GetPrimaryPz(){return PrimaryPz;} // Vector containing pz of all primary particles
        std::vector<double> * GetBackgroundDecayTime(){return BackgroundDecayTime;} // Vector containing decay time of all background particles
        std::vector<int> * GetBackgroundAtomicNumber(){return BackgroundAtomicNumber;} // Vector containing atomic number of all background particles
        std::vector<int> * GetBackgroundAtomicMass(){return BackgroundAtomicMass;} // Vector containing atomic mass of all background particles
        std::vector<double> * GetBackgroundVertexX(){return BackgroundVertexX;} // Vector containing interaction X position of all background particles
        std::vector<double> * GetBackgroundVertexY(){return BackgroundVertexY;} // Vector containing interaction Y position of all background particles
        std::vector<double> * GetBackgroundVertexZ(){return BackgroundVertexZ;} // Vector containing interaction Z position of all background particles



        TTree * GetInputTree(){return input_tree;}

    private:

        TFile * input_file;
        TFile * OutputFile;
        TTree * input_tree;
        TTree * OutputTree;

        //--------------------------------------------------
        // new branch variables
        //--------------------------------------------------
        int runID;
        std::vector<std::vector<double>> SavePhotons;
        int generated_counter;
        int nPhotons;

        //--------------------------------------------------
        // existing branch variables
        //--------------------------------------------------

        std::string TreeName_;

        // Vectors containing information of the hits in the event 
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
        std::vector<int> *pdg = new std::vector<int>();


        // Vectors containing information of the initial/primary particles in the event 
        std::vector<int> InitialPDG; 
        std::vector<int> PrimaryPDG;
        std::vector<double> * InteractionTime = new std::vector<double>(); 
        std::vector<double> * InitialEnergy = new std::vector<double>();
        std::vector<double> * PrimaryEnergy = new std::vector<double>();
        std::vector<double> * PrimaryPx = new std::vector<double>();
        std::vector<double> * PrimaryPy = new std::vector<double>();
        std::vector<double> * PrimaryPz = new std::vector<double>();
        std::vector<double> * InteractionVertexX = new std::vector<double>();
        std::vector<double> * InteractionVertexY = new std::vector<double>();
        std::vector<double> * InteractionVertexZ = new std::vector<double>();

        

        std::vector<double> * BackgroundDecayTime = new std::vector<double>(); 
        std::vector<int> * BackgroundAtomicNumber = new std::vector<int>(); 
        std::vector<int> * BackgroundAtomicMass = new std::vector<int>(); 
        std::vector<double> * BackgroundVertexX = new std::vector<double>(); 
        std::vector<double> * BackgroundVertexY = new std::vector<double>(); 
        std::vector<double> * BackgroundVertexZ = new std::vector<double>(); 

        int InitialPDG_;
        int PrimaryPDG_;
        double InteractionTime_;
        double InitialEnergy_;
        double PrimaryEnergy_;
        double PrimaryPx_;
        double PrimaryPy_;
        double PrimaryPz_;
        double BackgroundDecayTime_;
        int BackgroundAtomicNumber_;
        int BackgroundAtomicMass_;
        double InteractionVertexX_;
        double InteractionVertexY_;
        double InteractionVertexZ_;
        double BackgroundVertexX_;
        double BackgroundVertexY_;
        double BackgroundVertexZ_;

        //--------------------------------------------------
        // set branch addresses
        //--------------------------------------------------
        void set_branch_addresses(TTree *);
        void Initialize(char*, char*, std::string);

};

#endif
