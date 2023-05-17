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
        ROOTFileManager(char* inputFileName, char* outputFileName);
        ~ROOTFileManager();
        
        unsigned int NumberEntries();
        void EventFill();
        void Save();
        void GetEvent(int);
        void EventReset();
        int NEntries();
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
        TTree * GetTree(){return input_tree;}

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
        
        //--------------------------------------------------
        // set branch addresses
        //--------------------------------------------------
        void set_branch_addresses(TTree *);
        void Initialize(char*, char*);

};

#endif
