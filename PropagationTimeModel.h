#ifndef PROPAGATIONTIMEMODEL_H
#define PROPAGATIONTIMEMODEL_H

// PropagationTimeModel
//  - Contains functions to calculate the propagation time of scintillation photons

// March 2023 by A. Sánchez Castillo

#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <nlohmann/json.hpp>

#include "SemiAnalyticalModel.h"

#include "TVector3.h"
#include "TF1.h"
#include "TRandom.h"

#include <array>
#include <vector>

using json = nlohmann::json;


class PropagationTimeModel {

public:
  // constructor
  PropagationTimeModel(json OpParams);

  // propagation time
 void propagationTime(std::vector<double>& arrival_time_dist,
                       SemiAnalyticalModel::Point_t const& x0,
                       const size_t OpChannel);
  // Scintillation time
  double ScintTime();


private:
  // parameter and geometry initialization
  void Initialization(json OpParams);

  // direct / VUV light
  void getVUVTimes(std::vector<double>& arrivalTimes,
                   const double distance_in_cm,
                   const size_t angle_bin);

  //Set scintillation properties
  void SetScintillation();

  void generateParam(const size_t index, const size_t angle_bin);

  // reflected / visible light
  void getVISTimes(std::vector<double>& arrivalTimes,
                   const TVector3& ScintPoint,
                   const TVector3& OpDetPoint);

  // utility functions
  double fast_acos(double x) const;

  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     double x,
                     bool extrapolate,
                     size_t i = 0) const;

  void interpolate3(std::array<double, 3>& inter,
                    const std::vector<double>& xData,
                    const std::vector<double>& yData1,
                    const std::vector<double>& yData2,
                    const std::vector<double>& yData3,
                    double x,
                    bool extrapolate);
      
  static double finter_d(const double* x, const double* par);

  static double timing_model(const double* x, const double* par);

  static double ScintillationFunction(const double* x, const double* par);

  // photodetector geometry properties
  const std::string geometryfile;
  size_t nOpDets;
  std::vector<SemiAnalyticalModel::Point_t> fOpDetCenter;
  std::vector<double> fOpDetLength;
  std::vector<double> fOpDetHeight;

  // For VUV propagation time parametrization
  double fstep_size, fmax_d, fmin_d, fvuv_vgroup_mean, fvuv_vgroup_max, fangle_bin_timing_vuv;
  std::vector<std::vector<double>> fparameters[7];
  // vector containing generated VUV timing parameterisations
  std::vector<std::vector<TF1>> fVUV_timing;
  // vector containing min and max range VUV timing parameterisations are sampled to
  std::vector<std::vector<double>> fVUV_max;
  std::vector<std::vector<double>> fVUV_min;

  //Scintillation properties
  double tau_fast;
  double tau_slow;
  double scint_ratio;
  TF1 *ScintFunct;
  };

#endif
