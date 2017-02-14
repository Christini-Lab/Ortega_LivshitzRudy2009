/*** INTRO
 * Action potential duration calculator
 *
 * Calculate_APD.hpp, v1.0
 *
 * Author: Francis A. Ortega (v1.0) (2017)
 *
 *** NOTES
 *
 * Accepts a stream of membrane voltage values in real time, one at a time, to
 * calculate action potential duration if an AP is detected.
 *
 ***/

// Header Guard
#ifndef CALCULATE_APD_HPP
#define CALCULATE_APD_HPP

#include <vector>

class APD_Calculator {
 public:
  /*
    Constructor
    -- vRm: initial resting membrane potential value (mV)
    -- n_dt: time interval of data points (ms)
    -- psw: peak search window after a peak is found (ms)
    -- ut: upstroke threshold to denote start of AP (mV)
    -- repol: repolarization percentage for end of AP (%)
  */
  APD_Calculator(double vRm = -40, double n_dt = 0.1, int psw = 5, double ut = -40, int repol = 90);
  ~APD_Calculator() {};

  /*
    Mode of action potential duration algorithm
  */
  enum class APDMode {SEARCH, AMPLITUDE, END};

  /*
    Push membrane voltage to algorithm for real-time APD calculation
    -- voltage: membrane voltage at a single time point (mV)
  */
  APDMode push_voltage(double voltage);

  /*
    Reset algorithm and apd data vector
  */
  void reset();

  /*
    Return current mode of APD algorithm
  */
  APDMode get_mode() { return mode; };

  /*
    Return most recent action potential duration calculated
  */
  double get_apd();

  /*
    Set time increment of data
    -- new_dt: time increment (ms)
  */
  void set_dt(double new_dt) { dt = new_dt; };

  /*
    Set upstroke threshold
    -- new_ut: upstroke threshold (mV)
  */
  void set_upstrokeThreshold(double new_ut) { upstrokeThreshold = new_ut; };

  /*
    Set APD repolarization %
    -- new_repol: action potential duration repolarization percentage
  */
  void set_repolarizationPercent(double new_repol)
  { repolarizationPercent = new_repol; };


 private:
  APDMode mode = APDMode::SEARCH; // default start in SEARCH mode
  double dt; // (ms)
  std::vector<double> apd; // (ms)
  double voltage; // (mV)
  int peakSearchWindow; // (ms)
  double vmRest; // Resting membrane potential, used to calculate amp (mV)
  double upstrokeThreshold; // (mV)
  int repolarizationPercent; // (%)
  double peakVoltage; // (mV)
  int peakTime;
  double downstrokeThreshold; // (mV)
  double apAmplitude; // (mV)
  int steps; // Number of steps from start to end of action potential
};

#endif // CALCULATE_APD_HPP
