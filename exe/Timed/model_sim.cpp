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
 * main.cpp, v4.0
 *
 * Author: Francis Ortega (v1.0 - 4.0) (2011 -2017)
 *
 *** NOTES
 *
 * Example main which runs simple 500 beat pace simulation and outputs it to
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
#include <chrono>

#include "../../include/LivRudy2009.hpp"
#include "../../include/APD_Calculator.hpp"

int main() {
  using namespace std::chrono;
  LivRudy2009 model;

  std::cout << "Starting simulation" << std::endl;
  double voltage;
  double stim = 0;
  double dataDt = 0.2; // Dt of data output
  double dataTime = 0;
  double maxDt = 0.1; // Adaptive timestep maximum dt - 10 kHz
  double minDt = 0.001; // Adaptive timestep minimum dt - 1000 kHz
  double dt = maxDt; // Dt changes each time step, initially set to max
  model.setDt(dt); // Set model dt

  // Changes in voltage and intracellular concentrations determine adaptive
  // changes in dt
  double dVdt;
  double dCaidt;
  double dNaidt;
  double dKidt;

  // If Vm or intracellular concentrations changes less than this, set dt to max
  double dXdtThresh = maxDt * 2;
  double maxDx; // Biggest change between voltage or ions

  double v0 = model.getVm(); // Previous timestep voltage
  double cai0 = model.getCai(); // Previous timestep intracellular calcium
  double nai0 = model.getNai(); // Previous timestep intracellular sodium
  double ki0 = model.getKi(); // Previous timestep intracellular potassium

  int steps = dataDt / dt;
  int bcl = 500; // ms
  int beats = 1000;
  int beatIdx = 0;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  int protocolLength = beats * bcl / dataDt;

  // Unitless to prevent rounding errors
  int bclCounter = std::round(bcl / dataDt);
  int stimCounter = std::round(stimLength / dataDt);

  std::vector<double> timePerBeat;

  high_resolution_clock::time_point t_start, t_end;
  duration<double, std::milli> t_span;

  // Each time increment is equivalent to dataDt
  for (int time = 0; time < protocolLength; time++) {
    // Voltage and intracellular concentrations
    // dXdt for each
    dVdt = std::abs(model.getVm() - v0) / dataDt; // Calculate dVdt
    dCaidt = std::abs(model.getCai() - cai0) / dataDt * 1e4; // Calculate dVdt
    dNaidt = std::abs(model.getNai() - nai0) / dataDt; // Calculate dVdt
    dKidt = std::abs(model.getKi() - ki0) / dataDt; // Calculate dVdt
    // Previous timestep values
    v0 = model.getVm();
    cai0 = model.getCai();
    nai0 = model.getNai();
    ki0 = model.getKi();

    // Stimulation and adaptive dt calculation
    if (time % bclCounter < stimCounter) { // Stimulation on
      if (time % bclCounter == 0) // Start of beat
        t_start = high_resolution_clock::now();
      stim = -1 * stimAmp;
      dt = minDt; // Set min dt to minimum Dt during stimulation
      steps = dataDt / dt;
      model.setDt(dt);
    }
    else { // Stimulation off
      stim = 0;
      // Adaptive dt calculation - Fast changes in volage or ions will result in
      // smaller dt
      // Calcium is scaled due to its relatively small concentration compared
      // to other ions

      // Find maximum change
      maxDx = std::max({dVdt, dCaidt, dNaidt, dKidt});

      // If changes are slow, set to maximum Dt
      if (maxDx < dXdtThresh) {
        dt = maxDt; // Set to maximum dt
        steps = dataDt / dt; // Calculate intregration steps
        model.setDt(dt); // Set model dt
      }
      // If changes are fast, lower dt
      else {
        steps = std::ceil(maxDx / dXdtThresh);
        if (steps > dataDt / minDt)
          steps = dataDt / minDt;

        dt = dataDt / steps; // Calculate dt based on desired integration steps
        model.setDt(dt); // Set model dt
      }
    }

    for (int i = 0; i < steps; i++) {
      model.iClamp(stim);
    }

    // End of beat, calculate computation time taken
    if (time % bclCounter == bclCounter - 1) {
      t_end = high_resolution_clock::now();
      t_span = t_end - t_start;
      timePerBeat.push_back(t_span.count());
    }

    if (model.getStatus()) { // If model did not crash
      // Increment data time
      dataTime += dataDt;
    }
    else { // Model crash
      std::cout << "ERROR: Model crash" << std::endl;
      break;
    }
  }

  // Total time in seconds
  double totalTime =
      std::accumulate(timePerBeat.begin(), timePerBeat.end(), 0.0) / 1000;
  // Average time in seconds
  double avgTime = totalTime / timePerBeat.size();

  std::cout << "Simulation finished" << std::endl;
  std::cout << "Max time per beat (ms): "
            << *std::max_element(timePerBeat.begin(), timePerBeat.end())
            << std::endl
            << "Min time per beat (ms): "
            << *std::min_element(timePerBeat.begin(), timePerBeat.end())
            << std::endl
            << "Average time per beat (s): "
            << avgTime << std::endl
            << "Total time (s) (" << timePerBeat.size() << " beats): "
            << totalTime << std::endl;
}
