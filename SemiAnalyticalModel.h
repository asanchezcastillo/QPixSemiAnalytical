#ifndef SEMIANALYTICALMODEL_H
#define SEMIANALYTICALMODEL_H

// SemiAnalyticalModel
//  - Contains functions to calculate the number of photons reaching each photon detector, along with the required utilities.

// March 2023 by A. Sánchez Castillo

#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <nlohmann/json.hpp>
#include "TRandom3.h"

using json = nlohmann::json;

class SemiAnalyticalModel {

public:
  // constructor
  SemiAnalyticalModel(json params);

  struct Point_t {
    double x; // x coord
    double y; // y coord
    double z; // z coord
  };

  void detectedDirectVisibilities(std::vector<double>& DetectedVisibilities, Point_t const& ScintPoint) const;
  void detectedNumPhotons(std::vector<int>& DetectedNumPhotons, const std::vector<double>& OpDetVisibilities, const int NumPhotons) const;

    double GetDistance(Point_t const& ScintPoint, int opDet) const;

  double GetAngle(Point_t const& ScintPoint, int opDet) const;

private:

  // parameter and geometry initialization
  void Initialization(json OpParams);

  // structure for rectangular solid angle calculation
  struct Dims {
    double h, w; // height, width
  };

  // Vector structure
  struct Vector_t {
    double x; // x coord
    double y; // y coord
    double z; // z coord
  };

  // structure for optical detector information
  struct OpticalDetector {
  double h; // height
  double w; // width
  Point_t OpDetPoint;
  };

  //Function to check overlaps in the geometry
  void CheckOverlaps() const;

  // direct light photo-detector visibility calculation
  double VUVVisibility(Point_t const& ScintPoint, OpticalDetector const& opDet) const;

  // Gaisser-Hillas
  double Gaisser_Hillas(const double x, const double* par) const;



  // solid angle calculations
  // rectangular aperture
  double Rectangle_SolidAngle(const double a, const double b, const double d) const;
  double Rectangle_SolidAngleXZ(Dims const& o,Vector_t const& v) const;
  double Rectangle_SolidAngleYZ(Dims const& o,Vector_t const& v) const;

  double fast_acos(double x) const;

  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     double x,
                     bool extrapolate,
                     size_t i = 0) const;

  double interpolate2(const std::vector<double>& xDistances,
                      const std::vector<double>& rDistances,
                      const std::vector<std::vector<std::vector<double> > >& parameters,
                      const double x,
                      const double r,
                      const size_t k) const;

  const Point_t fanode_centre;
  const std::string geometryfile;

  // photodetector geometry properties
  const bool fDetectorPlaneXZ;
  const size_t nOpDets;
  Dims fcathode_plane;
  Dims fanode_plane;
  std::vector<Point_t> fOpDetCenter;
  std::vector<double> fOpDetLength;
  std::vector<double> fOpDetHeight;

  const int fvuv_absorption_length;

  // For VUV semi-analytic hits
  double fdelta_angulo_vuv;

  std::vector<std::vector<double > > fGHparams;
  std::vector<std::vector<double>> fGHparams_border_angle;
  std::vector<std::vector<double> > fGHparams_border;

  // correction parameters for VIS Nhits estimation
  double fdelta_angulo_vis;
  double fAnodeReflectivity;
  // flat PDs
  std::vector<double> fvis_distances_x_flat;
  std::vector<double> fvis_distances_r_flat;
  std::vector<std::vector<std::vector<double> > > fvispars_flat;
};

#endif

