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

  // Devenyi & Ortega 2016 GPig modifications
  // model.setGNa(model.getGNa() * 0.967);
  // model.setGNab(model.getGNab() * 1.99);
  // model.setGCaL(model.getGCaL() * 0.525);
  // model.setGCaT(model.getGCaT() * 0.322);
  // model.setGCab(model.getGCab() * 1.559);
  // model.setGK1(model.getGK1() * 1.073);
  // model.setGKr(model.getGKr() * 1.782);
  // model.setGKs(model.getGKs() * 0.042);
  // model.setGKp(model.getGKp() * 0.059);
  // model.setGNaK(model.getGNaK() * 1.645);
  // model.setGNCX(model.getGNCX() * 1.096);
  // model.setGpCa(model.getGpCa() * 0.270);
  // model.setGserca(model.getGserca() * 2.14);

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
  int beats = 120 * 15; // 15 minutes
  int beatIdx = 0;
  int stimAmp = 40; // pA/pF
  int stimLength = 1; // ms

  int protocolLength = beats * bcl / dataDt;


  // Time, voltage, current, and intracellular concentrations for last 5 beats
  std::vector< std::vector<double> >
      data(6, std::vector<double>(5 * bcl / dataDt));
  int offset = protocolLength - data.at(0).size();


  // Diastolic voltage and intracellular concentrations per beat
  std::vector< std::vector<double> > dia_data(5, std::vector<double>(beats));
  std::vector< std::vector<double> > sys_data(5, std::vector<double>(beats));

  // Unitless to prevent rounding errors
  int bclCounter = bcl / dataDt;
  int stimCounter = stimLength / dataDt;

  // APD calculation class initialization
  APD_Calculator apdCalc(model.getVm(), dataDt);

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
      // Before stimulus starts, gather diastolic data
      if (time % bclCounter == 0) {
        std::cout << "Beat: " << beatIdx + 1 << std::endl;
        dia_data.at(0).at(beatIdx) = beatIdx + 1;
        dia_data.at(1).at(beatIdx) = model.getVm();
        dia_data.at(2).at(beatIdx) = model.getNai();
        dia_data.at(3).at(beatIdx) = model.getKi();
        dia_data.at(4).at(beatIdx) = model.getCai();
      }
      stim = -1 * stimAmp;
      dt = minDt; // Set min dt to minimum Dt during stimulation
      steps = dataDt / dt;
      model.setDt(dt);
    }
    else { // Stimulation off
      if (time % bclCounter == stimCounter) {
        // At end of stimulus, gather systolic data
        sys_data.at(0).at(beatIdx) = beatIdx + 1;
        sys_data.at(1).at(beatIdx) = model.getVm();
        sys_data.at(2).at(beatIdx) = model.getNai();
        sys_data.at(3).at(beatIdx) = model.getKi();
        sys_data.at(4).at(beatIdx) = model.getCai();
        beatIdx++;
      }
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

    if (model.getStatus()) { // If model did not crash, save data
      if (time >= offset) {
        data.at(0).at(time - offset) = dataTime; // Time
        data.at(1).at(time - offset) = model.getVm(); // Voltage
        data.at(2).at(time - offset) = stim; // Injected current
        data.at(3).at(time - offset) = model.getNai(); // Intrac Na
        data.at(4).at(time - offset) = model.getKi(); // Intra K
        data.at(5).at(time - offset) = model.getCai(); // Intra Ca
        // Increment data time
        dataTime += dataDt;
      }
      // Push voltage to APD calculator
      apdCalc.push_voltage(model.getVm());
    }
    else { // Model crash
      std::cout << "ERROR: Model crash" << std::endl;
      break;
    }
  }
  std::cout << "Simulation finished" << std::endl;
  std::cout << "Min voltage: "
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
  std::ofstream dataFile;

  // Trace data output
  dataFile.open("trace_data.dat");
  dataFile << "Time,Voltage,Current,Nai,Ki,Cai" << std::endl;
  dataFile << std::setprecision(8);
  for (int idx = 0; idx < data.at(0).size(); idx++)
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
  // Vector will be empty if no APs are detected
  std::vector<double> apdData(apdCalc.get_all_apd());
  dataFile.open("apd.dat");
  dataFile << "APD" << std::endl;
  dataFile << std::setprecision(6);
  std::copy(apdData.begin(), apdData.end(),
            std::ostream_iterator<double>(dataFile,"\n"));
  dataFile.close();
}
