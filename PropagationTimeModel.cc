#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <nlohmann/json.hpp>

#include "PropagationTimeModel.h"
#include "SemiAnalyticalModel.h"

#include "TMath.h"

#include <iostream>
#include <chrono>
#include <fstream>

#include "TF1.h"

// constructor
PropagationTimeModel::PropagationTimeModel(json OpParams)
  :   nOpDets{OpParams["nOpDet"]},
  fstep_size{OpParams["StepSize"]},
  fmax_d{OpParams["max_d"]},
  fmin_d{OpParams["min_d"]},
  tau_fast{OpParams["tau_fast"]},
  tau_slow{OpParams["tau_slow"]},
  scint_ratio{OpParams["scint_ratio"]},
  fvuv_vgroup_mean{OpParams["VUVGroupMean"]},
  fvuv_vgroup_max{OpParams["VUVGroupMax"]},
  finflexion_point_distance{OpParams["inflexion_point_distance"]},
  fangle_bin_timing_vuv{OpParams["angle_bin_timing"]},
  geometryfile{OpParams["geometry_file"]}
{
  // initialise parameters and geometry
  Initialization(OpParams);
}

// initialization
void PropagationTimeModel::Initialization(json OpParams)
{

  // Initialize the parameters for the Propagation time. I'm sure there's a way more efficient way of doing it.
   for (size_t i = 0; i < std::size(OpParams["Distances_landau"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Distances_landau"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Distances_landau"][i][j]));
   }
   fparameters[0].push_back(row_params);
  } //end par0

    for (size_t i = 0; i < std::size(OpParams["Norm_over_entries"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Norm_over_entries"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Norm_over_entries"][i][j]));
   }
   fparameters[1].push_back(row_params);
  } //end par1

    for (size_t i = 0; i < std::size(OpParams["Mpv"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Mpv"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Mpv"][i][j]));
   }
   fparameters[2].push_back(row_params);
  } //end par2

    for (size_t i = 0; i < std::size(OpParams["Width"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Width"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Width"][i][j]));
   }
   fparameters[3].push_back(row_params);
  } //end par3

    for (size_t i = 0; i < std::size(OpParams["Distances_exp"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Distances_exp"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Distances_exp"][i][j]));
   }
   fparameters[4].push_back(row_params);
  } //end par2


    for (size_t i = 0; i < std::size(OpParams["Slope"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Slope"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Slope"][i][j]));
   }
   fparameters[5].push_back(row_params);
  } //end par5

    for (size_t i = 0; i < std::size(OpParams["Expo_over_Landau_norm"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Expo_over_Landau_norm"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Expo_over_Landau_norm"][i][j]));
   }
   fparameters[6].push_back(row_params);
  } //end par6
    
    // create vector of empty TF1s that will be replaced with the sampled
    // parameterisations that are generated as they are required
    const size_t num_params =
      (fmax_d - fmin_d) /
      fstep_size; // for d < fmin_d, no parameterisaton, a delta function is used instead
    const size_t num_angles = std::round(90 / fangle_bin_timing_vuv);
    fVUV_timing = std::vector(num_angles, std::vector(num_params, TF1()));

    // initialise vectors to contain range parameterisations sampled to in each case
    // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling
    // is regenerated, this is the slow part!
    fVUV_max = std::vector(num_angles, std::vector(num_params, 0.0));
    fVUV_min = std::vector(num_angles, std::vector(num_params, 0.0));

    std::cout << "Generating timing parameters: " << std::endl;

    // generate VUV parameters
    
    for (size_t angle_bin = 0; angle_bin < num_angles; ++angle_bin) {
      for (size_t index = 0; index < num_params; ++index) {
        std::cout << "index " << index << "angle bin " << angle_bin<<std::endl;
        generateParam(index, angle_bin);

      }
    }

 // get PDS information
 fOpDetCenter.reserve(nOpDets);
 fOpDetLength.reserve(nOpDets);
 fOpDetHeight.reserve(nOpDets);
 std::ifstream infile;
 std::vector<int> OpChannel;
 std::vector<double> pos_x;
 std::vector<double> pos_y;
 std::vector<double> pos_z;
 std::vector<double> detector_w;
 std::vector<double> detector_h;
 int ID;
 double posx;
 double posy;
 double posz;
 double detectorw;
 double detectorh;
 infile.open(geometryfile);
 if(infile.fail()) // checks to see if file opended 
  {
    std::cout << "Error: the file " << geometryfile << " could not be opened" << std::endl;
  }
    while(infile>>ID>>posx>>posy>>posz>>detectorw>>detectorh) 
  {
    OpChannel.push_back(ID);
    pos_x.push_back(posx);
    pos_y.push_back(posy);
    pos_z.push_back(posz);
    detector_w.push_back(detectorw);
    detector_h.push_back(detectorh);
  }
    infile.close();

  for (size_t i=0 ; i < nOpDets; i++) 
  {
    // Get detector information of all detectors (position and sizes)
    SemiAnalyticalModel::Point_t center{pos_x.at(i),pos_y.at(i),pos_z.at(i)};
    fOpDetCenter.push_back(center);
    fOpDetLength.push_back(detector_w.at(i));
    fOpDetHeight.push_back(detector_h.at(i));
  }

  SetScintillation();
 
}


//......................................................................
// Propagation time calculation function
void PropagationTimeModel::propagationTime(std::vector<double>& arrival_time_dist,
                                           SemiAnalyticalModel::Point_t const& x0,
                                           const size_t OpChannel)
{
  // Get VUV photons transport time distribution from the parametrization
  SemiAnalyticalModel::Point_t const& opDetCenter = fOpDetCenter[OpChannel];
  double distance =
    std::hypot(x0.x - opDetCenter.x, x0.y - opDetCenter.y, x0.z - opDetCenter.z);
  double cosine;
  cosine = std::abs(x0.y - opDetCenter.y) / distance;
  double theta = fast_acos(cosine) * 180. / 3.1415;
  int angle_bin = theta / fangle_bin_timing_vuv;
  getVUVTimes(arrival_time_dist, distance, angle_bin); // in ns
}

//......................................................................
// VUV propagation times calculation function
void PropagationTimeModel::getVUVTimes(std::vector<double>& arrivalTimes,
                                       const double distance,
                                       const size_t angle_bin)
{
  if (distance < fmin_d) {
    // times are fixed shift i.e. direct path only
    double t_prop_correction = distance / fvuv_vgroup_mean;
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] = t_prop_correction;
    }
  }
  else {
    // determine nearest parameterisation in discretisation
    int index = std::round((distance - fmin_d) / fstep_size);
    // randomly sample parameterisation for each
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
      arrivalTimes[i] = fVUV_timing[angle_bin][index].GetRandom(fVUV_min[angle_bin][index],
                                                                fVUV_max[angle_bin][index]);
    }
  }
}

//......................................................................
// VUV propagation times parameterization generation function
void PropagationTimeModel::generateParam(const size_t index, const size_t angle_bin)
{

  // get distance
  double distance_in_cm = (index * fstep_size) + fmin_d;

  // time range
  const double signal_t_range = 5000.;

  // parameterisation TF1
  TF1 VUVTiming;

  // direct path transport time
  double t_direct_mean = distance_in_cm / fvuv_vgroup_mean;
  double t_direct_min = distance_in_cm / fvuv_vgroup_max;

  // Defining the model function(s) describing the photon transportation timing vs distance
  // Getting the landau parameters from the time parametrization

  std::array<double, 3> pars_landau;
  interpolate3(pars_landau,
               fparameters[0][0],
               fparameters[2][angle_bin],
               fparameters[3][angle_bin],
               fparameters[1][angle_bin],
               distance_in_cm,
               true);
  // Deciding which time model to use (depends on the distance)
  // defining useful times for the VUV arrival time shapes
  if (distance_in_cm >= finflexion_point_distance) {
    double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
    // Set model: Landau
    VUVTiming = TF1("VUVTiming", model_far, 0, signal_t_range, 4);
    VUVTiming.SetParameters(pars_far);
  }
  else {
    // Set model: Landau + Exponential
    VUVTiming = TF1("VUVTiming", model_close, 0, signal_t_range, 7);
    // Exponential parameters
    double pars_expo[2];
    // Getting the exponential parameters from the time parametrization
    pars_expo[1] = interpolate(fparameters[4][0], fparameters[5][angle_bin], distance_in_cm, true);
    pars_expo[0] = interpolate(fparameters[4][0], fparameters[6][angle_bin], distance_in_cm, true);
    pars_expo[0] *= pars_landau[2];
    pars_expo[0] = std::log(pars_expo[0]);
    // this is to find the intersection point between the two functions:
    TF1 fint = TF1("fint", finter_d, pars_landau[0], 4 * t_direct_mean, 5);
    double parsInt[5] = {
      pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
    fint.SetParameters(parsInt);
    double t_int = fint.GetMinimumX();
    double minVal = fint.Eval(t_int);
    // the functions must intersect - output warning if they don't
    if (minVal > 0.015) {
    
      std::cout << "WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
    }

 
    double parsfinal[7] = {t_int,
                           pars_landau[0],
                           pars_landau[1],
                           pars_landau[2],
                           pars_expo[0],
                           pars_expo[1],
                           t_direct_min};
    VUVTiming.SetParameters(parsfinal);
  }

  // set the number of points used to sample parameterisation
  // for shorter distances, peak is sharper so more sensitive sampling required
  int fsampling;
  if (distance_in_cm < 50)
    fsampling = 10000;
  else if (distance_in_cm < 100)
    fsampling = 5000;
  else
    fsampling = 1000;
  VUVTiming.SetNpx(fsampling);

  // calculate max and min distance relevant to sample parameterisation
  // max
  const size_t nq_max = 1;
  double xq_max[nq_max];
  double yq_max[nq_max];
  xq_max[0] = 0.975; // include 97.5% of tail
  VUVTiming.GetQuantiles(nq_max, yq_max, xq_max);
  double max = yq_max[0];
  // min
  double min = t_direct_min;

  // store TF1 and min/max, this allows identical TF1 to be used every time sampling
  // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
  // all subsequent calls check if it has been generated previously and are ~100+ times quicker
  fVUV_timing[angle_bin][index] = VUVTiming;
  fVUV_max[angle_bin][index] = max;
  fVUV_min[angle_bin][index] = min;
}

//......................................................................

double PropagationTimeModel::fast_acos(double x) const
{
  double negate = double(x < 0);
  x = std::abs(x);
  x -= double(x > 1.0) * (x - 1.0); // <- equivalent to min(1.0,x), but faster
  double ret = -0.0187293;
  ret = ret * x;
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * std::sqrt(1.0 - x);
  ret = ret - 2. * negate * ret;
  return negate * 3.14159265358979 + ret;
}

//......................................................................
// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly
// monotonic increasing boolean argument extrapolate determines behaviour
// beyond ends of array (if needed)
double PropagationTimeModel::interpolate(const std::vector<double>& xData,
                                         const std::vector<double>& yData,
                                         double x,
                                         bool extrapolate,
                                         size_t i) const
{
  if (i == 0) {
    size_t size = xData.size();
    if (x >= xData[size - 2]) { // special case: beyond right end
      i = size - 2;
    }
    else {
      while (x > xData[i + 1])
        i++;
    }
  }
  double xL = xData[i];
  double xR = xData[i + 1];
  double yL = yData[i];
  double yR = yData[i + 1]; // points on either side (unless beyond ends)
  if (!extrapolate) {       // if beyond ends of array and not extrapolating
    if (x < xL) return yL;
    if (x > xR) return yL;
  }
  const double dydx = (yR - yL) / (xR - xL); // gradient
  return yL + dydx * (x - xL);               // linear interpolation
}

//......................................................................
void PropagationTimeModel::interpolate3(std::array<double, 3>& inter,
                                        const std::vector<double>& xData,
                                        const std::vector<double>& yData1,
                                        const std::vector<double>& yData2,
                                        const std::vector<double>& yData3,
                                        double x,
                                        bool extrapolate)
{
  size_t size = xData.size();
  size_t i = 0;               // find left end of interval for interpolation
  if (x >= xData[size - 2]) { // special case: beyond right end
    i = size - 2;
  }
  else {
    while (x > xData[i + 1])
      i++;
  }
  double xL = xData[i];
  double xR = xData[i + 1]; // points on either side (unless beyond ends)
  double yL1 = yData1[i];
  double yR1 = yData1[i + 1];
  double yL2 = yData2[i];
  double yR2 = yData2[i + 1];
  double yL3 = yData3[i];
  double yR3 = yData3[i + 1];

  if (!extrapolate) { // if beyond ends of array and not extrapolating
    if (x < xL) {
      inter[0] = yL1;
      inter[1] = yL2;
      inter[2] = yL3;
      return;
    }
    if (x > xR) {
      inter[0] = yL1;
      inter[1] = yL2;
      inter[2] = yL3;
      return;
    }
  }
  const double m = (x - xL) / (xR - xL);
  inter[0] = m * (yR1 - yL1) + yL1;
  inter[1] = m * (yR2 - yL2) + yL2;
  inter[2] = m * (yR3 - yL3) + yL3;
}

//......................................................................
double PropagationTimeModel::finter_d(const double* x, const double* par)
{
  double y1 = par[2] * TMath::Landau(x[0], par[0], par[1]);
  double y2 = TMath::Exp(par[3] + x[0] * par[4]);

  return TMath::Abs(y1 - y2);
}

//......................................................................
double PropagationTimeModel::model_close(const double* x, const double* par)
{
  // par0 = joining point
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par4 = Expo cte
  // par5 = Expo tau
  // par6 = t_min

  double y1 = par[3] * TMath::Landau(x[0], par[1], par[2]);
  double y2 = TMath::Exp(par[4] + x[0] * par[5]);
  if (x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
  if (x[0] < par[0]) y2 = 0.;

  return (y1 + y2);
}

//......................................................................
double PropagationTimeModel::model_far(const double* x, const double* par)
{
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par0 = t_min

  double y = par[3] * TMath::Landau(x[0], par[1], par[2]);
  if (x[0] <= par[0]) y = 0.;

  return y;
}

//......................................................................
//Sets parameter for the ScintillationTime function
void PropagationTimeModel::SetScintillation()
{
ScintFunct=new TF1("ScintFunc", "([0]/[1])*exp(-(x/[1]))+((1-[0])/[2])*exp(-(x/[2]))" , 0.001, 5000.);
double SlowComp = scint_ratio/(1+scint_ratio);
ScintFunct->SetParameter(0, SlowComp);
ScintFunct->SetParameter(1, tau_slow);
ScintFunct->SetParameter(2, tau_fast);
return;
}

//......................................................................
// Returns random number from ScintillationTime distribution
double PropagationTimeModel::ScintTime()
{
gRandom->SetSeed(0);
return ScintFunct->GetRandom ();
}



