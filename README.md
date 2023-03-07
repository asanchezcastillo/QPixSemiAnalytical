# QPixSemiAnalytical

C++ implementation of SemiAnalytical model for QPix-Light Simulation.

An introduction
---------------

This code is a C++ implementation of the [semi-analytical model](https://arxiv.org/abs/2010.00324) for scintillation light simulation in a LArTPC. 
The code takes as an input a root file produced by qpixg4 code containing the information related to the energy depositions for an event and simulates the propagation and detection of the scintillation photons.

Requirements
---------------
Before building the project you will need to meet the following requirements:
```
cmake version 3.23.2
GNU Make 3.81
ROOT Version: 6.26/06
nlohmann-json 3.11.2 Library to manage json files. Can be downloded from (https://github.com/nlohmann/json)
gcc-12
g++-12
libomp 15.0.7
```
WARNING: Root needs to be configured and built with the same C++ standard as the program that will make use of it. This code requires std-c++17 so you will need to compile root with this version.

Building
---------------
First we need to clone the repository. Open a new terminal, go to the directory you want to have the code in and run:

```
git clone https://github.com/asanchezcastillo/QPixSemiAnalytical.git
```

Now you have to export the path to you gcc and g++ compilers by sourcing setup.sh. Note that you might need to change the paths inside the file: 

```
source setup.sh
```

Now you can make the project and create the executable file. 

```
mkdir build
cd build
cmake ..
make
```

Parameters file
---------------
Is is required to provide a "params.json" file containing all the relevant input parameters for the simulation. The required parameters are:
```
"DetectorPlaneXZ": whether the detector plane is the XZ plane. 
geometry_file: name of the geometry file. 
nOpDet: number of optical detectors.
vuv_absorption_length: LAr VUV absorption length in cm.
GH_params: Array of arrays containing GH parameters for each angle bin.
StepSize: Step size in cm to generate the timing parameterization functions.
max_d:
min_d:
inflexion_point_distance: Distance at which the timing sampling switches from a landau + exponential distribution to only an exponential.
VUVGroupMean: VUV photons mean group velocity in cm/ns.
VUVGroupMax: VUV photons max group velocity in cm/ns.
angle_bin_timing: Angle binning for timing parameterization in deg.
tau_fast: LAr fast scintillation constant in ns.
tau_slow: LAr slow scintillation constant in ns.
scint_ratio: LAr slow to fast scintillation ratio.
Distances_exp: Timing parameter.
Distances_landau: Timing parameter.
Expo_over_Landau_norm: Timing parameter.
Mpv: Timing parameter.
Norm_over_entries: Timing parameter.
Width: Timing parameter.
EF: Electric field in the TPC in KeV/cm.
fWion: Ionization work function for LAr in MeV.
fWph: Ion+excitation work function for LAr in MeV.
fRecombA: LArQL model parameter.
fRecombk: LArQL model parameter.
fLarqlChi0A: LArQL model parameter.
fLarqlChi0B: LArQL model parameter.
fLarqlChi0C: LArQL model parameter.
fLarqlChi0D: LArQL model parameter.
fLarqlAlpha: LArQL model parameter.
fLarqlBeta: LArQL model parameter.
```

The parameter file has to be in the same folder as the executable, so you will have to move it to the build folder.

Geometry file
---------------
It is required to provide a geometry file under the build directory containing the optical detector's IDs, their central positions and their heights and widths. The information should be ordered in colums as follows:
```
ID XPosition YPosition ZPosition Heigh Width
0  <xpos_0>  <ypos_0>  <zpos_0>    H.    W.
1  <xpos_1>  <ypos_1>  <zpos_1>    H.    W.   
.
.
.
```

The first raw is only for ilustration purposes and should not be included in the actual geometry file. The name of the geometry file has to be provided in the parameters file "params.json". At the moment, the simulation admits two different configurations: one with the detectors in the XZ plane, corresponding to a DUNE-VD-like geometry and another with the detectors in the YZ plane, as in the DUNE-HD geometry.


Root input files
---------------
It is also required a .root input file containing the information on the energy depositions of the events that are to be simulated. The file should contain an "event" tree with the following branches:



Running the code
---------------
Finally, the code is ready to be run as follows:
```
./main --file <input_file_name>
```

Examining the code
---------------
Now were will go through the code to show how it is structured and the most relevant bits. In the first place, we have the SemiAnalytical class, which is in charge of computing the number of photons that reach each photon detector. Its two most relevant functions are: 
```c++
void SemiAnalyticalModel::detectedDirectVisibilities(std::vector<double>& DetectedVisibilities, Point_t const& ScintPoint) const
```
This function computes the visibility of each optical detector from the scintillation point.
```c++
void SemiAnalyticalModel::detectedNumPhotons(std::vector<int>& DetectedNumPhotons, const std::vector<double>& OpDetVisibilities,const int NumPhotons) const
```
This computes the number of photons reaching each detector. 

The second relevant class is PropagationTimeModel. When an object of this class is initialized a vector of TF1 objects containing a set of generated timing Ladau+Exponential functions is generated (this usually takes ~40s). The most relevant function of this class is:
```c++
void PropagationTimeModel::getVUVTimes(std::vector<double>& arrivalTimes, const double distance, const size_t angle_bin)
```
This function computes the number arrival times from the semi-analytical model parameterization.  

Code workflow
---------------
Having described the most important clases of the code, we can roughly understand its workflow. The first relevant part is to initalize an object of the SemiAnalyticalModel and the PropagationTimeModel which will be used for the relevant calculations:

```c++
std::unique_ptr<SemiAnalyticalModel> semi;
semi = std::make_unique<SemiAnalyticalModel>(OpParams); //Initialize SemiAnalyticalModel object
std::unique_ptr<PropagationTimeModel> PropTime;
PropTime = std::make_unique<PropagationTimeModel>(OpParams); //Initialize PropagationTimeModel object
```
Now, for each event we will loop over energy depositions, defining for each one an object of the class EnergyDeposition, from which we will get the relevant information for the semi-analytical model calculations:
```c++
// Initialize the energy deposition object with its StartPoint and the EndPoint:
SemiAnalyticalModel::Point_t StartPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
SemiAnalyticalModel::Point_t EndPoint{hitX_start->at(nHit), hitY_start->at(nHit), hitZ_start->at(nHit)};
std::unique_ptr<EnergyDeposition> Edep;
Edep = std::make_unique<EnergyDeposition>(OpParams, edep->at(nHit), StartPoint, EndPoint, time_start->at(nHit), time_end->at(nHit) ,length->at(nHit));
```
With this information we can already compute the number of photons that will reach every optical channel:
```c++
semi->detectedDirectVisibilities(OpDetVisibilities, ScintPoint); //Compute visibility for each optical detector.
double nphot=Edep->LArQL(); //Compute number of photons generated with LArQL model.
semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot); //Compute the number of photons detected by each channel.
```
Finally, we can loop over each optical channel to sample the arrival times of the detected photons:
With this information we can already compute the number of photons that will reach every optical channel:
```c++
PropTime->propagationTime(transport_time, ScintPoint, channel);
for (size_t i = 0; i < n_detected; ++i) //Loop over detected photons
{
    int time;
    time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+PropTime->ScintTime() ); 
    ++photonHitCollection[channel].DetectedPhotons[time]; //Add an entry to [OpChannel,time]
}
```
After looping over all the hits, we will have an output root file containing an object:
```c++
std::vector<std::vector<double>> SimPhoton;
```
The first dimension refers to each one of the optical channel, whereas the second dimension stores each timetick (1ttick=1ns) at which a photon is detected.

