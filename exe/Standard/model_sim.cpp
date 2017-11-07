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

#include "../../include/LivRudy2009.hpp"
#include "../../include/APD_Calculator.hpp"

int main() {
  LivRudy2009 model;

  std::cout << "Starting simulation" << std::endl;
  double voltage;
  double stim = 0;
  double dataDt = 0.2; // Dt of data output
  double dataTime = 0;
  double dt = 0.001;
  model.setDt(dt);
  int steps = dataDt / dt;
  int bcl = 500; // ms
  int beats = 10;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  int protocolLength = beats * bcl / dataDt;
  std::vector< std::vector<double> >
      data(5, std::vector<double>(protocolLength));

  // Unitless to prevent rounding errors
  int bclCounter = bcl / dataDt;
  int stimCounter = stimLength / dataDt;

  // APD calculation class initialization
  APD_Calculator apdCalc(model.getVm(), dataDt);

  // Each time increment is equivalent to dataDt
  for (int time = 0; time < protocolLength; time++) {

    if (time % bclCounter < stimCounter)
      stim = -1 * stimAmp;
    else
      stim = 0;

    for (int i = 0; i < steps; i++) {
      model.iClamp(stim);
    }

    if (model.getStatus()) { // If model did not crash, save data
      data.at(0).at(time) = dataTime; // Time
      data.at(1).at(time) = model.getVm(); // Voltage
      data.at(2).at(time) = stim; // Intracellular Na
      data.at(3).at(time) = model.getKi(); // Intracellular K
      data.at(4).at(time) = model.getCai(); // Intracellular Ca

      // Push voltage to APD calculator
      apdCalc.push_voltage(model.getVm());

      // Increment data time
      dataTime += dataDt;
    }
    else { // Model crash
      std::cout << "ERROR: Model crash" << std::endl;
      break;
    }
  }

  std::cout << "Simulation finished" << std::endl
            << "Min voltage: "
            << *std::min_element(data.at(1).begin(), data.at(1).end())
            << std::endl
            << "Max voltage: "
            << *std::max_element(data.at(1).begin(), data.at(1).end())
            << std::endl
            << "Simulation length: "
            << data.at(0).size() << std::endl;

  std::vector<double> conditions(model.getConditions());
  std::cout << "Conditions at end of simulation: " << std::endl;
  for (auto it = conditions.begin(); it != conditions.end(); it++) {
    std::cout << *it << std::endl;
  }

  // Data output
  std::ofstream dataFile("data.dat");
  dataFile << "Time,Voltage,Nai,Ki,Cai" << std::endl;
  dataFile << std::setprecision(8);
  for (int idx = 0; idx < protocolLength; idx++)
    dataFile << data.at(0).at(idx) << ","
             << data.at(1).at(idx) << ","
             << data.at(2).at(idx) << ","
             << data.at(3).at(idx) << ","
             << data.at(4).at(idx) << std::endl;
  dataFile.close();

  // APD data output, seperate output since it is by beat
  std::vector<double> apdData(apdCalc.get_all_apd());
  dataFile.open("apd.dat");
  dataFile << std::setprecision(6);
  std::copy(apdData.begin(), apdData.end(),
            std::ostream_iterator<double>(dataFile,"\n"));
  dataFile.close();
}
