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

        ROOTFileManager(std::string const& input_file, std::string const& output_file);
        ~ROOTFileManager();
        
        unsigned int NumberEntries();
        void EventFill();
        void EventReset(int fNOpChannels);
        void Save();
        void GetEvent(int);

    
    private:

        
        TFile * tfile_;
        TTree * ttree_;

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
        // metadata
        //--------------------------------------------------



        //--------------------------------------------------
        // initialize
        //--------------------------------------------------
        void Initialize(std::string const&, std::string const&);

        //--------------------------------------------------
        // set branch addresses
        //--------------------------------------------------
        void set_branch_addresses(TTree *);

};

#endif
