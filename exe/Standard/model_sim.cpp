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

#include "../../include/LivRudy2009.hpp"
#include "../../include/APD_Calculator.hpp"

int main() {
  LivRudy2009 model;

  std::cout << "Starting simulation" << std::endl;
  double voltage;
  double stim = 0;
  double dataDt = 0.1; // Dt of data output
  double dataTime = 0;
  double dt = 0.01; // Dt of model integration
  model.setDt(dt);
  int steps = std::round(dataDt / dt);
  int bcl = 500; // ms
  int beats = 50;
  int beatIdx = 0;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  int protocolLength = std::round(beats * bcl / dataDt);
  std::vector< std::vector<double> >
      data(6, std::vector<double>(protocolLength));

  // Diastolic voltage, current, and intracellular concentrations per beat
  std::vector< std::vector<double> > dia_data(5, std::vector<double>(beats));
  std::vector< std::vector<double> > sys_data(5, std::vector<double>(beats));

  // Unitless to prevent rounding errors
  int bclCounter = std::round(bcl / dataDt);
  int stimCounter = std::round(stimLength / dataDt);

  // APD calculation class initialization
  APD_Calculator apdCalc(model.getVm(), dataDt);

  // Each time increment is equivalent to dataDt
  for (int time = 0; time < protocolLength; time++) {
      // Stimulation
      // Before stimulus starts, gather diastolic data
      if (time % bclCounter == 0) {
        std::cout << "Beat: " << beatIdx + 1 << std::endl;
        dia_data.at(0).at(beatIdx) = beatIdx + 1;
        dia_data.at(1).at(beatIdx) = model.getVm();
        dia_data.at(2).at(beatIdx) = model.getNai();
        dia_data.at(3).at(beatIdx) = model.getKi();
        dia_data.at(4).at(beatIdx) = model.getCai();
        stim = -1 * stimAmp;
      }
      else if (time % bclCounter == stimCounter) {
        // At end of stimulus, gather systolic data
        sys_data.at(0).at(beatIdx) = beatIdx + 1;
        sys_data.at(1).at(beatIdx) = model.getVm();
        sys_data.at(2).at(beatIdx) = model.getNai();
        sys_data.at(3).at(beatIdx) = model.getKi();
        sys_data.at(4).at(beatIdx) = model.getCai();
        beatIdx++;
        stim = 0;
      }

    // Simulate model at desired dt, broken up in steps based on data Dt
    for (int i = 0; i < steps; i++) {
      model.iClamp(stim);
    }
    if (model.getStatus()) { // If model did not crash, save data
      data.at(0).at(time) = dataTime; // Time
      data.at(1).at(time) = model.getVm(); // Voltage
      data.at(2).at(time) = stim; // Injected current
      data.at(3).at(time) = model.getNai(); // Intracellular Na
      data.at(4).at(time) = model.getKi(); // Intracellular K
      data.at(5).at(time) = model.getCai(); // Intracellular Ca

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
  // Voltage and intracellular concentration traces
  std::ofstream dataFile("trace_data.dat");
  dataFile << "Time,Voltage,Current,Nai,Ki,Cai" << std::endl;
  dataFile << std::setprecision(8);
  for (int idx = 0; idx < protocolLength; idx++)
    dataFile << data.at(0).at(idx) << ","
             << data.at(1).at(idx) << ","
             << data.at(2).at(idx) << ","
             << data.at(3).at(idx) << ","
             << data.at(4).at(idx) << ","
             << data.at(5).at(idx) << std::endl;
  dataFile.close();

  // Diastolic concentration
  dataFile.open("dia_data.dat");
  dataFile << "Beat,Voltage,Nai,Ki,Cai" << std::endl;
  dataFile << std::setprecision(8);
  for (int idx = 0; idx < beats; idx++)
    dataFile << dia_data.at(0).at(idx) << ","
             << dia_data.at(1).at(idx) << ","
             << dia_data.at(2).at(idx) << ","
             << dia_data.at(3).at(idx) << ","
             << dia_data.at(4).at(idx) << std::endl;
  dataFile.close();

  // Systolic concentration
  dataFile.open("sys_data.dat");
  dataFile << "Beat,Voltage,Nai,Ki,Cai" << std::endl;
  dataFile << std::setprecision(8);
  for (int idx = 0; idx < beats; idx++)
    dataFile << sys_data.at(0).at(idx) << ","
             << sys_data.at(1).at(idx) << ","
             << sys_data.at(2).at(idx) << ","
             << sys_data.at(3).at(idx) << ","
             << sys_data.at(4).at(idx) << std::endl;
  dataFile.close();

  // APD data output, seperate output since it is by beat
  std::vector<double> apdData(apdCalc.get_all_apd());
  dataFile.open("apd.dat");
  dataFile << "APD" << std::endl;
  dataFile << std::setprecision(6);
  std::copy(apdData.begin(), apdData.end(),
            std::ostream_iterator<double>(dataFile,"\n"));
            dataFile.close();
}
