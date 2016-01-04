/*
 * Copyright (C) 2015 Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*** INTRO
 *
 * 2009 Livshitz Rudy Model of a ventricular guinea pig myocyte
 * Biophysical Journal, 2009
 *
 * main.cpp, v3.0
 *
 * Author: Francis Ortega (v1.0 - 3.0) (2011 -2015)
 *
 *** NOTES
 *
 * Example main which runs simple 5 beat pace simulation and outputs it to
 * 'voltage.dat'.
 *
 ***/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>

#include "include/LivRudy2009.hpp"

int main() {
  LivRudy2009 model;

  std::cout << "Starting simulation" << std::endl;
  double voltage;
  double stim = 0;
  double dataDt = 0.1; // Dt of data output
  double maxDt = model.getDt(); // Maximum dt for adaptive timestep
  double dt = model.getDt(); // Adaptive timestep
  double dVdt; // dVdt used to modify timestep
  double dVdtThresh = maxDt * 2; // If Vm changes less than this, set dt to max
  double v0 = model.getVm(); // Previous timestep voltage
  int steps = dt / dataDt;
  int bcl = 500; // ms
  int beats = 105;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  int protocolLength = beats * bcl / dt;
  std::vector<double> voltageData;
  voltageData.reserve(protocolLength);

  // Unitless to prevent rounding errors
  int bclCounter = bcl / dt;
  int stimCounter = stimLength / dt;

  // Each time increment is equivalent to dataDt
  for (int time = 0; time < protocolLength; time++) {
    dVdt = std::abs(model.getVm() - v0) / dataDt; // Calculate dVdt
    v0 = model.getVm();

    // Stimulation
    if (time % bclCounter <= stimCounter) {
      stim = -1 * stimAmp;
    }
    else
      stim = 0;

    // Adaptive dt calculation
    if (dVdt < dVdtThresh) {
      dt = maxDt;
      model.setDt(dt);
      steps = dataDt / dt;
    }
    else { // dVdt is > than dVdtThresh, so reduce dt
      steps = std::ceil(dVdt / dVdtThresh); // Round up to integer

      if (steps > maxDt / 0.001)
        steps = maxDt / 0.001; // Set min dt to 1000kHz

      dt = maxDt / steps;
      model.setDt(dt);
    }

    for (int i = 0; i < steps; i++) {
      model.iClamp(stim);
    }

    if (model.getStatus()) // If model did not crash, save voltage
      voltageData.push_back(model.getVm());
    else { // Model crash
      std::cout << "ERROR: Model crash" << std::endl;
      break;
    }
  }

  std::cout << "Simulation finished" << std::endl <<
      "Min element: " <<
      *std::min_element(voltageData.begin(), voltageData.end()) << std::endl <<
      "Max element: " <<
      *std::max_element(voltageData.begin(), voltageData.end()) << std::endl <<
      "Simulation length: " <<
      voltageData.size() << std::endl;

  std::vector<double> conditions(model.getConditions());
  std::cout << "Conditions at end of simulation: " << std::endl;
  for (auto it = conditions.begin(); it != conditions.end(); it++) {
    std::cout << *it << std::endl;
  }

  // Data output
  std::ofstream dataFile("voltage.dat");
  dataFile << std::setprecision(12);
  std::copy(voltageData.begin(), voltageData.end(),
            std::ostream_iterator<double>(dataFile,"\n"));
  dataFile.close();
}
