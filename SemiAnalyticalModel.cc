#include "SemiAnalyticalModel.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>

// constructor
SemiAnalyticalModel::SemiAnalyticalModel(json OpParams)
    : fanode_centre{OpParams["anode_centre"][0],OpParams["anode_centre"][1],OpParams["anode_centre"][2]},
      fvuv_absorption_length{OpParams["vuv_absorption_length"]},
      nOpDets{OpParams["nOpDet"]},
      geometryfile{OpParams["geometry_file"]},
      fDetectorPlaneXZ{OpParams["DetectorPlaneXZ"]},
      fdelta_angulo_vuv{OpParams["DeltaAngle"]}

{
  Initialization(OpParams);
}

// initialization
void SemiAnalyticalModel::Initialization(json OpParams)
{
  //Read GH parameters from json file. (I suspect there might be a more optimal way of doing this)
  for (size_t i = 0; i < std::size(OpParams["GH_params"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["GH_params"][0]); j++ )
   {
     row_params.push_back( static_cast< double >(OpParams["GH_params"][i][j]));
   }
   fGHparams.push_back(row_params);
  }
  //Read border effect GH params. 
  for (size_t i = 0; i < std::size(OpParams["GH_border"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["GH_border"][0]); j++ )
   {
    row_params.push_back( static_cast< double >(OpParams["GH_border"][i][j]));
   }
   fGHparams_border.push_back(row_params);
  }
  //Read angle bin for border effects. 
  for (size_t i = 0; i < std::size(OpParams["GH_border_angle"]); i++) 
  {
   std::vector<double> row_params;
   for (size_t j=0; j<std::size(OpParams["GH_border_angle"][0]); j++ )
   {
    row_params.push_back( static_cast< double >(OpParams["GH_border_angle"][i][j]));
   }
   fGHparams_border_angle.push_back(row_params);
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
    while(infile>>ID>>posx>>posy>>posz>>detectorh>>detectorw) 
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
    Point_t center{pos_x.at(i),pos_y.at(i),pos_z.at(i)};
    fOpDetCenter.push_back(center);
    fOpDetLength.push_back(detector_w.at(i));
    fOpDetHeight.push_back(detector_h.at(i));
  }

  CheckOverlaps();
}
//......................................................................
// VUV semi-analytical model visibility calculation

void SemiAnalyticalModel::detectedDirectVisibilities(std::vector<double>& DetectedVisibilities, Point_t const& ScintPoint) const
{
 // double visibility_sum=0; 
  DetectedVisibilities.resize(nOpDets);
  for (size_t OpDet=0; OpDet<nOpDets; OpDet++)
  {
   const SemiAnalyticalModel::OpticalDetector op{fOpDetHeight[OpDet],
                                                fOpDetLength[OpDet],
                                                fOpDetCenter[OpDet]};

    DetectedVisibilities[OpDet] = VUVVisibility(ScintPoint, op);
  }

}

double SemiAnalyticalModel::VUVVisibility( Point_t const& ScintPoint, OpticalDetector const& opDet) const
{ 
  // distance and angle between ScintPoint and OpDetPoint
  Vector_t const relative = {ScintPoint.x - opDet.OpDetPoint.x, ScintPoint.y- opDet.OpDetPoint.y, ScintPoint.z - opDet.OpDetPoint.z };
  
  double distance = sqrt(pow(relative.x,2)+ pow(relative.y,2)+ pow(relative.z,2)); 
  double cosine;
  cosine = std::abs(relative.y) / distance;
  const double theta = fast_acos(cosine) * 180. / 3.1416;
  double solid_angle = 0.;

 // get scintillation point coordinates relative to detector window centre
  Vector_t const abs_relative{
  std::abs(relative.x), std::abs(relative.y), std::abs(relative.z)};

 if(fDetectorPlaneXZ) 
 {
  solid_angle = Rectangle_SolidAngleXZ(Dims{opDet.h, opDet.w}, abs_relative); 
 }

 else
 {
  solid_angle = Rectangle_SolidAngleYZ(Dims{opDet.h, opDet.w}, abs_relative);
 }
 // calculate visibility by geometric acceptance
 // accounting for solid angle and LAr absorbtion length

  double visibility_geo =
    std::exp(-1. * distance / fvuv_absorption_length) * (solid_angle / (4 * 3.1416));
  // apply Gaisser-Hillas correction for Rayleigh scattering distance
  // and angular dependence offset angle bin

  const size_t j = (theta / fdelta_angulo_vuv);

  // determine GH parameters, accounting for border effects
  
  double r = std::hypot(ScintPoint.x - fanode_centre.x, ScintPoint.z - fanode_centre.z);  // ONLY FOR DETECTOR XZ PLANE, HAS TO BE CHANGED
  
  double pars_ini[4] = {0, 0, 0, 0};

  double s1 = 0;
  double s2 = 0;
  double s3 = 0;
  double s4 = 0;

  pars_ini[0] = fGHparams[0][j];
  pars_ini[1] = fGHparams[1][j];
  pars_ini[2] = fGHparams[2][j];
  pars_ini[3] = fGHparams[3][j];
 
  s1 = interpolate(fGHparams_border_angle[0], fGHparams_border[0], theta, true);
  s2 = interpolate(fGHparams_border_angle[0], fGHparams_border[1], theta, true);
  s3 = interpolate(fGHparams_border_angle[0], fGHparams_border[2], theta, true);
  s4 = interpolate(fGHparams_border_angle[0], fGHparams_border[3], theta, true);

  // add border correction to parameters
  pars_ini[0] = pars_ini[0] + s1 * r;
  pars_ini[1] = pars_ini[1] + s2 * r;
  pars_ini[2] = pars_ini[2] + s3 * r;
  pars_ini[3] = pars_ini[3] + s4 * r;
  
  
  // calculate correction
  double GH_correction = Gaisser_Hillas(distance, pars_ini);

  return GH_correction * visibility_geo / cosine;
}
//......................................................................
// Gaisser-Hillas function definition
double SemiAnalyticalModel::Gaisser_Hillas(const double x, const double* par) const
{
  
  double X_mu_0 = par[3];
  double Normalization = par[0];
  double Diff = par[1] - X_mu_0;
  double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
  double Exponential = std::exp((par[1] - x) / par[2]);

  return (Normalization * Term * Exponential);
  
}

//......................................................................
// solid angle of rectangular aperture
double SemiAnalyticalModel::Rectangle_SolidAngle(const double a,
                                                 const double b,
                                                 const double d) const
{
  double aa = a / (2. * d);
  double bb = b / (2. * d);
  double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
  return 4. * fast_acos(std::sqrt(aux));
}
//................... 
//Solid angle for detector in the plane YZ (VD geometry)

double SemiAnalyticalModel::Rectangle_SolidAngleXZ(SemiAnalyticalModel::Dims const& out, SemiAnalyticalModel::Vector_t const& v) const {


  if( v.x==0.0 && v.z==0.0){
    return Rectangle_SolidAngle(out.w,out.h,v.y);
  }

  if( (std::abs(v.x) > out.w/2.0) && (std::abs(v.z) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.x)-out.w/2.0;
    B = std::abs(v.z)-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.x) <= out.w/2.0) && (std::abs(v.z) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.x)+out.w/2.0;
    B = -std::abs(v.z)+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.x) > out.w/2.0) && (std::abs(v.z) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.x)-out.w/2.0;
    B = -std::abs(v.z)+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.x) <= out.w/2.0) && (std::abs(v.z) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.x)+out.w/2.0;
    B = std::abs(v.z)-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

//......................................................................
//Solid angle for a detector in the XZ plane (HD geometry)

double SemiAnalyticalModel::Rectangle_SolidAngleYZ(SemiAnalyticalModel::Dims const& out, SemiAnalyticalModel::Vector_t const& v) const {


  if( v.y==0.0 && v.z==0.0){
    return Rectangle_SolidAngle(out.w,out.h,v.x);
  }

  if( (std::abs(v.y) > out.w/2.0) && (std::abs(v.z) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.y)-out.w/2.0;
    B = std::abs(v.z)-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.z);
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.y) <= out.w/2.0) && (std::abs(v.z) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.y)+out.w/2.0;
    B = -std::abs(v.z)+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.y) > out.w/2.0) && (std::abs(v.z) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.y)-out.w/2.0;
    B = -std::abs(v.z)+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.y) <= out.w/2.0) && (std::abs(v.z) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.y)+out.w/2.0;
    B = std::abs(v.z)-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.y);
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}


//......................................................................

double SemiAnalyticalModel::fast_acos(double x) const
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
// calculates number of photons detected given visibility and emitted number of photons
void SemiAnalyticalModel::detectedNumPhotons(std::vector<int>& DetectedNumPhotons,
                                      const std::vector<double>& OpDetVisibilities,
                                      const int NumPhotons) const
{

  TRandom *random = new TRandom3(0);

  for (size_t i = 0; i < OpDetVisibilities.size(); ++i) {
   Double_t seed = OpDetVisibilities[i]*(double)NumPhotons;
   DetectedNumPhotons[i] = random->Poisson(seed);
  }
  delete random;
}


//......................................................................
// Checks for overlaps in the geometry. (Length = beam direction)
void SemiAnalyticalModel::CheckOverlaps() const
{
  std::cout << "Checking overlaps in the geometry: " << std::endl;

  if(!fDetectorPlaneXZ)
  {
   for(size_t i = 0; i<fOpDetCenter.size(); i++)
   {
    for(size_t j=i+1; j<fOpDetCenter.size(); j++)
    {
      if( ( (fOpDetCenter.at(i).z-(fOpDetHeight.at(i)/2))<(fOpDetCenter.at(j).z+(fOpDetHeight.at(j)/2)) 
      && (fOpDetCenter.at(j).z-(fOpDetHeight.at(j)/2))<(fOpDetCenter.at(i).z+(fOpDetHeight.at(i)/2))) &&
      ( (fOpDetCenter.at(i).y-(fOpDetLength.at(i)/2))<(fOpDetCenter.at(j).y+(fOpDetLength.at(j)/2)) 
      && (fOpDetCenter.at(j).y-(fOpDetLength.at(j)/2))<(fOpDetCenter.at(i).y+(fOpDetLength.at(i)/2)) )
      )
      {
      std::cout << "There is overlapping between pixels " << i << " and " << j <<  ". Aborting the job... " << std::endl;
      abort();
      }
    }
   }
  }

  else
  {
   for(size_t i = 0; i<fOpDetCenter.size(); i++)
   {
    for(size_t j=i+1; j<fOpDetCenter.size(); j++)
    {
      if( ( (fOpDetCenter.at(i).z-(fOpDetHeight.at(i)/2))<(fOpDetCenter.at(j).z+(fOpDetHeight.at(j)/2)) 
      && (fOpDetCenter.at(j).z-(fOpDetHeight.at(j)/2))<(fOpDetCenter.at(i).z+(fOpDetHeight.at(i)/2))) &&
      ( (fOpDetCenter.at(i).x-(fOpDetLength.at(i)/2))<(fOpDetCenter.at(j).x+(fOpDetLength.at(j)/2)) 
      && (fOpDetCenter.at(j).x-(fOpDetLength.at(j)/2))<(fOpDetCenter.at(i).x+(fOpDetLength.at(i)/2)) )
      )
      {
      std::cout << "There is overlapping between pixels " << i << " and " << j <<  ". Aborting the job... " << std::endl;
      abort();
      }
    }
   } 
  }
  

}

//......................................................................
// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly
// monotonic increasing boolean argument extrapolate determines behaviour
// beyond ends of array (if needed)

double SemiAnalyticalModel::interpolate(const std::vector<double>& xData,
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



double SemiAnalyticalModel::interpolate2(
  const std::vector<double>& xDistances,
  const std::vector<double>& rDistances,
  const std::vector<std::vector<std::vector<double>>>& parameters,
  const double x,
  const double r,
  const size_t k) const
{
  // interpolate in x for each r bin, for angle bin k
  const size_t nbins_r = parameters[k].size();
  std::vector<double> interp_vals(nbins_r, 0.0);
  {
    size_t idx = 0;
    size_t size = xDistances.size();
    if (x >= xDistances[size - 2])
      idx = size - 2;
    else {
      while (x > xDistances[idx + 1])
        idx++;
    }
    for (size_t i = 0; i < nbins_r; ++i) {
      interp_vals[i] = interpolate(xDistances, parameters[k][i], x, false, idx);
    }
  }
  // interpolate in r
  double border_correction = interpolate(rDistances, interp_vals, r, false);
  return border_correction;
}



double SemiAnalyticalModel::GetDistance(Point_t const& ScintPoint,  int opDet) const
{


   const SemiAnalyticalModel::OpticalDetector op{fOpDetHeight[opDet],
                                                fOpDetLength[opDet],
                                                fOpDetCenter[opDet]};

  Vector_t const relative = {ScintPoint.x - op.OpDetPoint.x, ScintPoint.y- op.OpDetPoint.y, ScintPoint.z - op.OpDetPoint.z };
  double distance = sqrt(pow(relative.x,2)+ pow(relative.y,2)+ pow(relative.z,2)); 
  return distance;

}


double SemiAnalyticalModel::GetAngle(Point_t const& ScintPoint, int opDet) const
{

   const SemiAnalyticalModel::OpticalDetector op{fOpDetHeight[opDet],
                                                fOpDetLength[opDet],
                                                fOpDetCenter[opDet]};


  Vector_t const relative = {ScintPoint.x - op.OpDetPoint.x, ScintPoint.y- op.OpDetPoint.y, ScintPoint.z - op.OpDetPoint.z };
  double distance = sqrt(pow(relative.x,2)+ pow(relative.y,2)+ pow(relative.z,2)); 

  double cosine;
  cosine = std::abs(relative.y) / distance;
  double theta = fast_acos(cosine) * 180. / 3.1416;
  return theta;

}