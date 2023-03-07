#ifndef SIMPHOTONS_H
#define SIMPHOTONS_H


//  SimPhotons
//  - Class used to store the SemiAnalyticalModel output.

// March 2023 by A. Sánchez Castillo

#include <limits> 
#include <map>
#include <string>
#include <vector>


class SimPhotons {
public:

  SimPhotons() = default;

  /// Constructor: associated to optical detector channel `chan`, and empty.
  SimPhotons(int chan) : OpChannel(chan) {}

  int OpChannel; ///< Optical detector channel associated to this data.

  /// Number of photons detected at each given time: time tick -> photons.
  std::map<int, int> DetectedPhotons;

}; 

#endif