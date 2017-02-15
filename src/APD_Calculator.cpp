#include "../include/APD_Calculator.hpp"
#include <iostream>
// Constructor
APD_Calculator::APD_Calculator(double vRm, double n_dt, int psw, double ut,
                               int repol) :
    vmRest(vRm), dt(n_dt), peakSearchWindow(psw), upstrokeThreshold(ut),
    repolarizationPercent(repol) {};

// Reset algorithm and apd data vector
void APD_Calculator::reset() {
  mode = APDMode::SEARCH;
  apd.clear();
}

// Return most recent action potential duration calculated
double APD_Calculator::get_apd() {
  return apd.back();
}

// Return most recent num action potentials duration calculated
std::vector<double> APD_Calculator::get_apd(int num) {
  std::vector<double> retval(apd.end() - num, apd.end());
  return retval;
}

// Push voltage in real-time to APD calculating algorithm
APD_Calculator::APDMode APD_Calculator::push_voltage(double voltage) {
  switch (mode) {
    case APDMode::SEARCH:
      // Lowest membrane voltage during search is assumed to be Vrm
      if (voltage <= vmRest) {
        vmRest = voltage;
      }
      // Search for upstroke threshold crossing denoting start of AP
      if (voltage >= upstrokeThreshold) {
        mode = APDMode::AMPLITUDE;
        steps = 0; // Used to determine AP duration
        peakVoltage = voltage; // Set starting peak voltage
      }
      break;

    case APDMode::AMPLITUDE:
      // Find peak of AP
      if (peakVoltage < voltage) {
        peakVoltage = voltage;
        peakTime = steps;
      }
      // Keep searching for peak until search window limit is reached
      else if ((steps - peakTime) * dt > peakSearchWindow) {
        double apAmp;
        apAmp = peakVoltage - vmRest;

        // Calculate downstroke threshold based on AP amp and repolarization %
        downstrokeThreshold =
            peakVoltage - (apAmp * (repolarizationPercent / 100.0));
        mode = APDMode::END;
      }
      break;

    case APDMode::END:
      // Search for downstroke threshold crossing denoting end of AP
      if (voltage <= downstrokeThreshold) {
        apd.push_back(steps * dt);
        mode = APDMode::SEARCH; // Return to search for next AP
      }
      break;
  }
  steps++; // Increment data counter
}
