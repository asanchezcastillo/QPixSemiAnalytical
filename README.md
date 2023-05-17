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

| Parameter name                 |    Type                                     |  Description |             
|-----------------------------|---------------------------------------------|--------------|  
| DetectorPlaneXY | bool | Whether the detector plane is the XZ plane (true) or in the YZ plane (false) |
| geometry_file | string | Name of the geometry file |
| nOpDet | int | Number of optical channels |
| vuv_absorption_length | double | LAr VUV absorption length (in cm) |
| GH_params | std::vector< std::vector<\double\> > | GH parameters for each angle bin |
| GH_border | std::vector< std::vector<\double\> > | GH border corrections for each angle bin |
| DeltaAngle | int | Size of the angular bin for GH curves (in deg) |
| StepSize | int | Step size to generate the timing parameterization functions (in cm) |
| max_d | int | Maximum distance to generate the timing parameterization functions (in cm) |
| min_d | int | Minimum distance to generate the timing parameterization functions (in cm) |
| VUVGroupMean | double | Mean VUV group velocity in LAr [Journal of Instrumentation, 15(09):P09009–P09009, sep 2020] |
| VUVGroupMax | double | Mean VUV group velocity in LAr [G4 simulations] |
| angle_bin_timing | int | Size of the angular bin for timing parameterization functions (in deg) |
| tau_fast | int | Fast scintillation timing constant (in ns) [EPL 91 (2010) 6, 62002, Europhys.Lett. 91 (2010) 6, 62002] |
| tau_slow | int | Slow scintillation timing constant (in ns) [EPL 91 (2010) 6, 62002, Europhys.Lett. 91 (2010) 6, 62002] |
| scint_ratio | double | Fast to slow scintillation ratio [[J Chem Phys vol 91 (1989) 1469]] |
| Distances_landau | std::vector< std::vector<\double\> > | Vector containing the edges of bin distances for timimng parameterization (in cm) |
| Mpv1 | std::vector< std::vector<\double\> > | Vector containing MPV of the first Landau (in cm) |
| Norm_over_entries1 | std::vector< std::vector<\double\> > | Vector containing normalization of the first Landau (in cm) |
| Width1 | std::vector< std::vector<\double\> > | Vector containing width of the first Landau (in cm) |
| Mpv2 | std::vector< std::vector<\double\> > | Vector containing MPV of the second Landau (in cm) |
| Norm_over_entries2 | std::vector< std::vector<\double\> > | Vector containing normalization of the second Landau (in cm) |
| Width2 | std::vector< std::vector<\double\> > | Vector containing width of the second Landau (in cm) |
| fWion | double | Ionization work function for LAr (in MeV) [F. Marinho et al 2022 JINST 17 C07009] |
| fWph | double | Ion+excitation work function for LAr (in MeV) [F. Marinho et al 2022 JINST 17 C07009] |
| fRecombA | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fRecombk | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlChi0A | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlChi0B | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlChi0C | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlChi0D | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlAlpha | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| fLarqlBeta | double | LArQL model parameter [F. Marinho et al 2022 JINST 17 C07009] |
| EF | double | Electric field (in KeV/cm) |


The references included are for the parameter contained in this repository. This file has to be in the same folder as the executable, so you may have to move it to the build folder.

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
```c++
input_tree->SetBranchAddress("run",&run);
input_tree->SetBranchAddress("hit_start_x",&hitX_start);
input_tree->SetBranchAddress("hit_end_x",&hitX_end);
input_tree->SetBranchAddress("hit_start_y",&hitY_start);
input_tree->SetBranchAddress("hit_end_y",&hitY_end);
input_tree->SetBranchAddress("hit_start_z",&hitZ_start);
input_tree->SetBranchAddress("hit_end_z",&hitZ_end);
input_tree->SetBranchAddress("hit_energy_deposit",&edep);
input_tree->SetBranchAddress("hit_start_t",&time_start);
input_tree->SetBranchAddress("hit_end_t",&time_end);
input_tree->SetBranchAddress("hit_length",&length);
```

Running the code
---------------
Finally, the code is ready to be run as follows:
```
./OpticalSimulation -i <input_file_name> -o <output_file_name> -n <number_of_events>
```

Note: the input file name and the output file name are mandatory parameters whereas the number of events is not. If no n is given, the whole input file will be processed.

Examining the code
---------------
Now we will go through the code to show how it is structured and the most relevant bits. In the first place, we have the SemiAnalytical class, which is in charge of computing the number of photons that reach each photon detector. Its two most relevant functions are: 
```c++
void SemiAnalyticalModel::detectedDirectVisibilities(std::vector<double>& DetectedVisibilities, Point_t const& ScintPoint) const
```
This function computes the visibility of each optical detector from the scintillation point.
```c++
void SemiAnalyticalModel::detectedNumPhotons(std::vector<int>& DetectedNumPhotons, const std::vector<double>& OpDetVisibilities,const int NumPhotons) const
```
This function computes the number of photons reaching each detector. The second relevant class is PropagationTimeModel. When an object of this class is initialized a vector of TF1 objects containing a set of generated timing Ladau+Exponential functions is generated (this usually takes ~40s). The most relevant function of this class is:
```c++
void PropagationTimeModel::getVUVTimes(std::vector<double>& arrivalTimes, const double distance, const size_t angle_bin)
```
This function computes the number arrival times from the semi-analytical model parameterization.  

Code workflow
---------------
Having described the most important clases of the code, we can roughly understand its workflow. The first relevant part is to initalize an object of the SemiAnalyticalModel and the PropagationTimeModel which will be used for the relevant calculations:

```c++
std::unique_ptr<SemiAnalyticalModel> semi;
//Initialize SemiAnalyticalModel object
semi = std::make_unique<SemiAnalyticalModel>(OpParams); 
std::unique_ptr<PropagationTimeModel> PropTime;
//Initialize PropagationTimeModel object
PropTime = std::make_unique<PropagationTimeModel>(OpParams); 
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
//Compute visibility for each optical detector.
semi->detectedDirectVisibilities(OpDetVisibilities, ScintPoint); 
//Compute number of photons generated with LArQL model.
double nphot=Edep->LArQL(); 
//Compute the number of photons detected by each channel.
semi->detectedNumPhotons(DetectedNum, OpDetVisibilities, nphot); 
```
Finally, we can loop over each optical channel to sample the arrival times of the detected photons:
With this information we can already compute the number of photons that will reach every optical channel:
```c++
PropTime->propagationTime(transport_time, ScintPoint, channel);
//Loop over detected photons
for (size_t i = 0; i < n_detected; ++i) 
{
    int time;
    time =  static_cast<int>( ( (Edep->TimeStart() + Edep->TimeEnd())/2 ) + transport_time[i]+PropTime->ScintTime() ); 
    //Add an entry to [OpChannel,time]
    ++photonHitCollection[channel].DetectedPhotons[time]; 
}
```
After looping over all the hits, we will have an output root file containing an object:
```c++
std::vector<std::vector<double>> SimPhoton;
```
The first dimension refers to each one of the optical channel, whereas the second dimension stores each timetick (1ttick=1ns) at which a photon is detected.

Analysis example
---------------
An example on how to analyze the output of the simulation can be found in a Jupyter Notebook under the Analysis folder. This folder contains an output file (Output.root) obtained for a one APA Dune-HD-like detector containing 960 pixels (24 in the beam direction and 40 in the Y direction). The output of the analysis notebook is an animation representing the number of photons detected by each pixel for each time tick (ns).