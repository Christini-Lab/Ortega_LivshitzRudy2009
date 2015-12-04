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
  double dt = model.getDt();
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

  for (int time = 0; time < protocolLength; time++) {
    if (time % bclCounter <= stimCounter)
      stim = -1 * stimAmp;
    else
      stim = 0;

    voltage = model.iClamp(stim);

    // Voltage check
    if( isnan(voltage) || voltage < -150 || voltage > 150 ) {
      std::cout << "Error: voltage out of range" << std::endl;
      break;
    }
    voltageData.push_back(voltage);
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
