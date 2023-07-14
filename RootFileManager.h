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
        std::vector<int> GetInitialPDG(){return InitialPDG;} // Vector containing pdgs of all intial particles
        std::vector<int> GetPrimaryPDG(){return PrimaryPDG;} // Vector containing pdgs of all primaries particles
        std::vector<double> * GetInteractionTime(){return InteractionTime;} // Vector containing interaction times of all intial particles
        std::vector<double> * GetInitialEnergy(){return InitialEnergy;} // Vector containing energies of all intial particles
        std::vector<double> * GetPrimaryEnergy(){return PrimaryEnergy;} // Vector containing energies of all primary particles
        std::vector<double> * GetPrimaryPx(){return PrimaryPx;} // Vector containing px of all primary particles
        std::vector<double> * GetPrimaryPy(){return PrimaryPy;} // Vector containing py of all primary particles
        std::vector<double> * GetPrimaryPz(){return PrimaryPz;} // Vector containing pz of all primary particles
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


        // Vectors containing information of the initial/primary particles in the event 
        std::vector<int> InitialPDG; 
        std::vector<int> PrimaryPDG;
        std::vector<double> * InteractionTime = new std::vector<double>(); 
        std::vector<double> * InitialEnergy = new std::vector<double>();
        std::vector<double> * PrimaryEnergy = new std::vector<double>();
        std::vector<double> * PrimaryPx = new std::vector<double>();
        std::vector<double> * PrimaryPy = new std::vector<double>();
        std::vector<double> * PrimaryPz = new std::vector<double>();

        int InitialPDG_;
        int PrimaryPDG_;
        double InteractionTime_;
        double InitialEnergy_;
        double PrimaryEnergy_;
        double PrimaryPx_;
        double PrimaryPy_;
        double PrimaryPz_;


        //--------------------------------------------------
        // set branch addresses
        //--------------------------------------------------
        void set_branch_addresses(TTree *);
        void Initialize(char*, char*, std::string);

};

#endif
