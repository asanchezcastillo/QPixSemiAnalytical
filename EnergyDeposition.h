
    #ifndef SIMENERGYDEPOSIT_H
    #define SIMENERGYDEPOSIT_H

    //  EnergyDeposition
    //  - Class used to store energy depositions. Has functions to return the number of scintillation photons produced and some properties of the energy deposition.

    // March 2023 by A. Sánchez Castillo

    #include "SemiAnalyticalModel.h"

    #include <vector>

    #include <nlohmann/json.hpp>

    using json = nlohmann::json;

    class EnergyDeposition {

    public:

    EnergyDeposition(json OpParams,
                        double e,
                        SemiAnalyticalModel::Point_t position_start,
                        SemiAnalyticalModel::Point_t position_end,
                        double time_start,
                        double time_end,
                        double edepdx)
        : edep(e)
        , start_position{position_start}
        , end_position{position_end}
        , start_time{time_start}
        , end_time{time_end}
        , dx{edepdx}
        , fWion{OpParams["fWion"]}
        , fWph{OpParams["fWph"]}
        , fRecombA{OpParams["fRecombA"]}
        , fRecombk{OpParams["fRecombk"]}
        , fLarqlChi0A{OpParams["fLarqlChi0A"]}
        , fLarqlChi0B{OpParams["fLarqlChi0B"]}
        , fLarqlChi0C{OpParams["fLarqlChi0C"]}
        , fLarqlChi0D{OpParams["fLarqlChi0D"]}
        , fLarqlAlpha{OpParams["fLarqlAlpha"]}
        , fLarqlBeta{OpParams["fLarqlBeta"]}
        , EF{OpParams["EF"]}
    {}


    double Energy() const { return edep; }

    double LArQL() const
    {

        double num_quanta;
        double num_ions;
        double num_electrons;
        double recomb;
        double num_photons;
        double EscapingFraction;
        double FieldCorrection;
        
        num_quanta = edep / fWph; 
        num_ions = edep / fWion;
        recomb = fRecombA / (1. + (edep/dx) * fRecombk / EF);
        EscapingFraction = fLarqlChi0A / (fLarqlChi0B + std::exp(fLarqlChi0C + fLarqlChi0D * (edep/dx)) );
        FieldCorrection = std::exp(-EF / (fLarqlAlpha * std::log(edep/dx) + fLarqlBeta));
        recomb = recomb + EscapingFraction*FieldCorrection;
        num_electrons = num_ions*recomb;
        num_photons = num_quanta - num_electrons;



        if(num_photons<0){   
         std::cout << "WARNING: the number of photons is < 0, ignoring this hit" << std::endl;;
         return 0;
         }
        return num_photons ;
    }

    SemiAnalyticalModel::Point_t PositionStart() const { return {start_position.x, start_position.y, start_position.z}; }
    SemiAnalyticalModel::Point_t PositionEnd() const { return {end_position.x, end_position.y, end_position.z}; }
    double TimeEnd() const { return start_time ; }
    double TimeStart() const { return end_time ; }
    double Dx() const { return dx ; }
    SemiAnalyticalModel::Point_t MidPoint() const {return { (end_position.x+start_position.x)/2, (end_position.y+start_position.y)/2, (end_position.z + start_position.z)/2} ;}

    private:

    float edep;            ///< energy deposition (MeV)
    SemiAnalyticalModel::Point_t start_position; ///< positions in (cm)
    SemiAnalyticalModel::Point_t end_position; ///< positions in (cm)
    double start_time; ///< (ns)
    double end_time; ///< (ns)
    double dx;

 
    // LArQL Parameters.
    double fWion;
    double fWph;
    double fRecombk;
    double fRecombA;
    double fLarqlChi0A;
    double fLarqlChi0B;
    double fLarqlChi0C;
    double fLarqlChi0D;
    double fLarqlAlpha;
    double fLarqlBeta;
    double EF;

    };

    #endif 