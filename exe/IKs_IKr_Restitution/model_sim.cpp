/*
 * Copyright (C) 2017 Weill Medical College of Cornell University
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
 * model_sim.cpp, v1.0
 *
 * Author: Francis Ortega (v1.1) (2017)
 *
 *** NOTES
 *
 * Restitution portrait simulation specifically for exploring the impact of
 * IKs/IKr ratio.
 *
 * Arguments: Start BCL, End BCL, Increment BCL, IKs conductance,
 * IKr conductance, data output file
 *
 * v1.1 - Added start, ending, and increment BCL as arguments
 *
 ***/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>

#include "../../include/LivRudy2009.hpp"
#include "../../include/APD_Calculator.hpp"

int main(int argc, char *argv[]) {
  // Arguments: Start BCL, End BCL, Increment BCL, IKs conductance,
  // IKr conductance, data output file
  if (argc != 7) {
    std::cout << "Error: invalid number of arguments: " << argc <<
        " instead of " << 6 << std::endl;
    return EXIT_FAILURE;
  }

  LivRudy2009 model;
  model.setGKs(atof(argv[4]));
  model.setGKr(atof(argv[5]));
  std::cout << "Starting restitution portrait simulation" << std::endl;
  double voltage;
  double stim = 0;
  double dataDt = 0.1; // Dt of data output
  double maxDt = model.getDt(); // Maximum dt for adaptive timestep
  double dt = model.getDt(); // Adaptive timestep
  double dVdt; // dVdt used to modify timestep
  double dVdtThresh = maxDt * 2; // If Vm changes less than this, set dt to max
  double v0 = model.getVm(); // Previous timestep voltage
  int steps = dt / dataDt;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  // Restitution portrait parameters
  int startBCL = atof(argv[1]); // Starting BCL (ms)
  int endBCL = atof(argv[2]); // Desired ending BCL (ms)
  int incrementBCL = atof(argv[3]); // BCL change (ms)
  int numBeatsPerBCL = 2000; // Number of stimulations for a BCL
  int numApdSave = 10; // Number of APDs to save for each BCL

  // Data structures for output
  std::vector<double> apdData;
  std::vector<double> bclData;

  int stimCounter = stimLength / dt; // Unitless to prevent rounding errors

  // APD calculation class initialization
  APD_Calculator apdCalc(model.getVm(), dt);

  bool done = false; // Flag for loop
  int bcl = startBCL;

  // Repeat protocol for every desired BCL, decrementing until end is reach
  while (!done) {
    // Initial BCL is 1000 beats to ensure steady-state of model
    int beats;
    if (bcl == startBCL)
      beats = 2000;
    else
      beats = numBeatsPerBCL;

    int beatCount = 0;
    int bclCounter = bcl / dt; // Unitless to prevent rounding errors
    int protocolLength = beats * bcl / dt;
    std::cout << "BCL: " << bcl << std::endl;
    // Each time increment is equivalent to dataDt
    for (int time = 0; time < protocolLength; time++) {
      dVdt = std::abs(model.getVm() - v0) / dataDt; // Calculate dVdt
      v0 = model.getVm();

      // Stimulation
      if (time % bclCounter <= stimCounter)
        stim = -1 * stimAmp;
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

      if (model.getStatus()) { // Check for model crash
        // Push voltage to APD calculator
        apdCalc.push_voltage(model.getVm());
      }
      else { // Model crash
        std::cout << "ERROR: Model crash" << std::endl;
        break;
      }
    }

    // Check for block, by checking if each stimulation produced an AP, if
    // not, skip BCL
    if (apdCalc.get_apd_length() != beats)
      break;
    // Save last x APDs along with corresponding BCL
    std::vector<double> tmp(apdCalc.get_apd(numApdSave));
    apdData.insert(apdData.end(), tmp.begin(), tmp.end());
    apdCalc.reset(); // Reset apd calculator data
    std::fill_n(std::back_inserter(bclData), numApdSave, bcl);

    if (bcl == endBCL) // Finished
      break;
    else if (startBCL > endBCL) // BCL intended to decrement
      bcl -= incrementBCL;
    else if (startBCL < endBCL) // BCL intended to increment
      bcl += incrementBCL;
  }

  std::cout << "Simulation finished" << std::endl;

  // APD and BCL data output
  std::ofstream dataFile(argv[6]);
  auto itA = apdData.begin();
  auto itB = bclData.begin();
  double ratio = model.getGKs() / model.getGKr();
  dataFile << "APD,BCL,Ratio" << std::endl;
  while (itA != apdData.end() || itB != bclData.end()) {
    dataFile << std::setprecision(8) << *itA << "," << *itB << ","
             << std::setprecision(2) << ratio << std::endl;
    ++itA;
    ++itB;
  }
  dataFile.close();
}
