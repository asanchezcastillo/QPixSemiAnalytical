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

#include <TApplication.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TFormula.h"
#include "TLegend.h"
#include "TRootCanvas.h"

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

    for (size_t i = 0; i < std::size(OpParams["Norm_over_entries1"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Norm_over_entries1"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Norm_over_entries1"][i][j]));
   }
   fparameters[1].push_back(row_params);
  } //end par1

    for (size_t i = 0; i < std::size(OpParams["Mpv1"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Mpv1"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Mpv1"][i][j]));
   }
   fparameters[2].push_back(row_params);
  } //end par2

    for (size_t i = 0; i < std::size(OpParams["Width1"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Width1"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Width1"][i][j]));
   }
   fparameters[3].push_back(row_params);
  } //end par3

    for (size_t i = 0; i < std::size(OpParams["Norm_over_entries2"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Norm_over_entries2"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Norm_over_entries2"][i][j]));
   }
   fparameters[4].push_back(row_params);
  } //end par2


    for (size_t i = 0; i < std::size(OpParams["Mpv2"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Mpv2"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Mpv2"][i][j]));
   }
   fparameters[5].push_back(row_params);
  } //end par5

    for (size_t i = 0; i < std::size(OpParams["Width2"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["Width2"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["Width2"][i][j]));
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
  else if (distance > (fmax_d-fmin_d))
  {
        // times are fixed shift i.e. direct path only
    double t_prop_correction = (distance-fmax_d) / fvuv_vgroup_mean;
    int max_index = std::round((fmax_d - fmin_d) / fstep_size)-1;
    for (size_t i = 0; i < arrivalTimes.size(); ++i) {
    arrivalTimes[i] = fVUV_timing[angle_bin][max_index].GetRandom(fVUV_min[angle_bin][max_index],
                                                                fVUV_max[angle_bin][max_index]) + t_prop_correction;
    }
  }
  else{
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

  //distance 


//distance 
fparameters[0]= {{20,  40,  60,  80, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780}},
//Norm1
fparameters[1]= {{2.30207, 1.43599, 1.06465 ,0.800104 , 0.60344741, 0.4524691 , 0.35670544, 0.287219, 0.2420945 , 0.20189513, 0.17357464, 0.15542911, 0.1417068 , 0.12901386, 0.11973669, 0.10995821, 0.10727, 0.0998892, 0.0987589, 0.08801661, 0.08497633, 0.0844696, 0.07927445, 0.0769224 , 0.07389544, 0.07135388, 0.06832368, 0.0676256 , 0.06552702, 0.06465223, 0.06228655, 0.06153403, 0.06135376, 0.06598879, 0.05748598, 0.0546869 , 0.05600351, 0.05299382, 0.05503482},
				{1.95656569, 1.51750074, 1.14246486, 0.82002433, 0.60344741, 0.4524691 , 0.35670544, 0.28529065, 0.2420945 , 0.20189513, 0.17770681, 0.15840516, 0.1417068 , 0.12901386, 0.11973669, 0.11203132, 0.11203132, 0.09965533, 0.09407799, 0.09407799, 0.08636203, 0.08292145, 0.07995421, 0.07732617, 0.07389544, 0.07135388, 0.06895246, 0.0676256 , 0.06552702, 0.06484652, 0.06228655, 0.06153403, 0.06135376, 0.06598879, 0.05748598, 0.0546869 , 0.05600351, 0.05337308, 0.05515088}} ;       

//MPV1
fparameters[2]= {{2.42653,  3.47639 ,   4.56213,  5.6563,  6.88886888, 8.17701693,  9.40446188, 10.5 , 12.00399106, 13.26065051, 14.43164457, 16.14445637, 18.00247545, 19.8008684 , 21.59126955, 22.70877703, 24.8765 , 26.8713, 28.4459 , 29.22631542, 31.29810915, 33.68749749, 34.49888527, 36.87638631, 38.83061361, 40.69070441, 40.59878376, 43.75153127, 43.73158678, 47.01081755, 48.01008117, 49.8888536 , 57.89347499, 67.2153492 , 54.49385121, 62.54398596, 59.49488718, 56.93540153, 62.0889571},	
     {2.57530694,  3.65498052,  4.76693438,  5.62376689,  6.88886888, 8.17701693,  9.40446188, 10.55087528, 12.00399106, 13.26065051, 14.96571284, 16.40432517, 18.00247545, 19.8008684 , 21.59126955, 23.69257219, 23.69257219, 27.39929068, 29.19801432, 29.19801432, 32.52865135, 33.68749749, 34.96707396, 37.29929936, 38.83061361, 40.69070441, 40.96970689, 43.75153127, 43.73158678, 47.30142044, 48.01008117, 49.8888536 , 57.89347499, 67.2153492 , 54.49385121, 62.54398596, 59.49488718, 57.44212267, 63.55451455}};
//Width1
fparameters[3]= {{ 0.405132, 0.715057  , 0.794829,  0.801921,  0.91809635, 1.08362321,  1.18470236,  1.24019,  1.49942309,  1.64501386, 1.73148413,  2.11828066,  2.60549293,  3.03728779,  3.44599235, 3.48106687,  4.04596,  4.62975 ,  4.67271,  4.68289077, 5.24792862,  5.64741,  5.81401534,  6.56762674,  7.15636693, 7.67559989,  6.89769476,  8.07012022,  7.44043407,  8.75380652, 8.7608927 ,  9.08966054, 13.52909985, 17.75545632, 10.03269989, 14.74869569, 11.80014622,  9.62879362, 11.8624133 }, 
      {0.405132,   0.715057  ,   0.794829,  0.79364367,  0.91809635, 1.08362321,  1.18470236,  1.2584553 ,  1.49942309,  1.64501386, 2.0376486 ,  2.24672272,  2.60549293,  3.03728779,  3.44599235, 4.07610696,  4.07610696,  5.0151109 ,  5.47651764,  5.47651764, 5.9421471 ,  6.00607494,  6.07154964,  6.81048944,  7.15636693, 7.67559989,  7.08439575,  8.07012022,  7.44043407,  8.91719706, 8.7608927 ,  9.08966054, 13.52909985, 17.75545632, 10.03269989, 14.74869569, 11.80014622,  9.9047782 , 12.74704914}};
//Norm2
fparameters[4]= {{ 4.33424, 2.54351, 2.2766, 1.58996, 0.95430689, 0.57181117, 0.46803775, 0.393179, 0.29903106, 0.24370628, 0.20063657, 0.16672473, 0.14433451, 0.13006406, 0.1205768 , 0.11263477, 0.10367, 0.100851, 0.0951923, 0.08834885, 0.08269245, 0.0807862, 0.07580457, 0.0711996 , 0.07340463, 0.07070518, 0.06563322, 0.06600685, 0.06200836, 0.05804076, 0.05812742, 0.05565013, 0.05891011, 0.05641617, 0.05265762, 0.05109714, 0.05400902, 0.05251859, 0.05042397}, 
			{4.33424, 2.49891227, 2.27300607, 1.30869072, 0.95430689, 0.57181117, 0.46803775, 0.39121557, 0.29903106, 0.24370628, 0.18929091, 0.16390257, 0.14433451, 0.13006406, 0.1205768 , 0.11313906, 0.11313906, 0.10070831, 0.09453016, 0.09453016, 0.08420473, 0.08174702, 0.07516402, 0.07005652, 0.07340463, 0.07070518, 0.06499705, 0.06600685, 0.06200836, 0.05741199, 0.05812742, 0.05565013, 0.05891011, 0.05641617, 0.05265762, 0.05109714, 0.05400902, 0.05225252, 0.05151845}};      
//MPV2
fparameters[5]= {{1.03714,  2.00353,  2.06978 ,  2.38434,  2.99998031, 2.98329836,  2.98016387,  3.20985,  4.19367485,  6.18020639, 8.81552314, 12.14976923, 15.58995328, 18.31715997, 20.14110958, 22.20097885, 24.6279, 25.9483, 28.3386, 31.13365044, 34.32510721, 36.0084, 39.07766887, 42.28579215, 40.35914617, 42.43898711, 47.61783458, 47.52557665, 52.25768809, 57.17322067, 58.15873788, 63.24926596, 57.96764431, 60.77305596, 70.26579503, 73.83680514, 66.97572519, 70.66642167, 75.20285239}, 	
        {1.03714,   2.00353,  2.35654548,  2.99841883,  2.99998031, 2.98329836,  2.98016387,  3.00355199,  4.19367485,  6.18020639, 10.15218086, 13.11891358, 15.58995328, 18.31715997, 20.14110958, 22.29081567, 22.29081567, 25.92787348, 28.50242642, 28.50242642, 32.47647839, 34.06688792, 39.55717389, 43.03439706, 40.35914617, 42.43898711, 48.35428839, 47.52557665, 52.25768809, 58.12383936, 58.15873788, 63.24926596, 57.96764431, 60.77305596, 70.26579503, 73.83680514, 66.97572519, 71.21265562, 73.03332708}}  ;     
//Width2
fparameters[6]= {{0.437868,  0.618856 ,  0.865159,  1.26165,  1.92012288, 3.15970414,  4.02486644,  4.87484,  6.17726712,  7.17385398, 8.16326006,  9.1261957 ,  9.88137648, 10.62151363, 11.41571545, 12.06539693, 12.9497, 13.368, 14.0301, 14.78020489, 15.35464201, 15.6605, 16.37449006, 17.08670431, 17.65828025, 18.24225165, 18.61054478, 19.13594809, 19.34770053, 19.76814248, 19.94755511, 19.85776546, 20.84869415, 21.66065011, 20.539715  , 20.34528892, 22.00190503, 21.97688903, 21.94240403},
 				{0.26490919,  0.55768921,  0.80338542,  1.30258629,  1.92012288, 3.15970414,  4.02486644,  4.95934088,  6.17726712,  7.17385398, 8.20359991,  8.94174676,  9.88137648, 10.62151363, 11.41571545, 11.96574798, 11.96574798, 13.37161568, 13.99810043, 13.99810043, 15.65199129, 16.16268843, 16.36651 , 17.12778883, 17.65828025, 18.24225165, 18.56775829, 19.13594809, 19.34770053, 19.67772715, 19.94755511, 19.85776546, 20.84869415, 21.66065011, 20.539715  , 20.34528892, 22.00190503, 21.90825038, 21.96136561}};


  std::array<double, 3> parsLandau1;
  std::array<double, 3> parsLandau2;
  interpolate3(parsLandau1,
               fparameters[0][0],
               fparameters[2][angle_bin],
               fparameters[3][angle_bin],
               fparameters[1][angle_bin],
               distance_in_cm,
               true);

  interpolate3(parsLandau2,
               fparameters[0][0],
               fparameters[5][angle_bin],
               fparameters[6][angle_bin],
               fparameters[4][angle_bin],
               distance_in_cm,
               true);

  // Deciding which time model to use (depends on the distance)
  // defining useful times for the VUV arrival time shapes

    // Set model: Landau + Landau
    VUVTiming = TF1("VUVTiming", timing_model, 0, signal_t_range, 8);

    // this is to find the intersection point between the two functions:
    TF1 fint = TF1("fint", finter_d, parsLandau1[0], 4 * t_direct_mean, 6);
    double parsInt[6] = {
      parsLandau1[0], parsLandau1[1], parsLandau1[2], parsLandau2[0], parsLandau2[1], parsLandau2[2]};
    fint.SetParameters(parsInt);

    // Intersection point for both landaus
    double t_int=fint.GetMinimumX(parsLandau1[0],parsLandau1[0]+18);
    double minVal = fint.Eval(t_int);

    // the functions must intersect - output warning if they don't
    if (minVal > 0.015) {
      std::cout << "WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
    }

    double parsfinal[8] = {t_int,
                           parsLandau1[0],
                           parsLandau1[1],
                           parsLandau1[2],
                           parsLandau2[0],
                           parsLandau2[1],
                           parsLandau2[2],
                           t_direct_min};
    VUVTiming.SetParameters(parsfinal);
  
  // set the number of points used to sample parameterisation
  // for shorter distances, peak is sharper so more sensitive sampling required
  int fsampling;
  if (distance_in_cm < 50)
    fsampling = 5000;
  else if (distance_in_cm < 100)
    fsampling = 5000;
  else
    fsampling = 5000;
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
  Double_t y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
  Double_t y2 = par[5]*TMath::Landau(x[0],par[3],par[4]);
  
  return TMath::Abs(y1 - y2);
}

double PropagationTimeModel::timing_model(const double* x, const double* par)
{
  // par0 = joining point
  // par1 = Landau1 MPV
  // par2 = Landau1 width
  // par3 = normalization 1
  // par4 = Landau2 MPV
  // par5 = Landau2 width
  // par6 = normalization 2
  // par7 = t_direct_min
 
  Double_t y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
  Double_t y2 = par[6]*TMath::Landau(x[0],par[4],par[5]);
 
  if (x[0] <= par[7] || x[0] > par[0]) y1 = 0.;
  if (x[0] < par[7] || x[0]< par[0]) y2 = 0.;

  return (y1 + y2);
}

//......................................................................
//Sets parameter for the ScintillationTime function
void PropagationTimeModel::SetScintillation()
{
ScintFunct=new TF1("ScintFunc", "([0]/[1])*exp(-(x/[1]))+((1-[0])/[2])*exp(-(x/[2]))" , 0.001, 5000.);
double SlowComp = 1/(1+scint_ratio);
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
