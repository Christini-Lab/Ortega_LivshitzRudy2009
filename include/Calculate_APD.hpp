/*** INTRO
 * Action potential duration calculator
 *
 * Calculate_APD.hpp
 *
 * Author: Francis A. Ortega (2017)
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

class CalculateAPD {
 public:
  CalculateAPD() {};
  ~CalculateAPD() {};

  /*
    Mode of action potential duration algorithm
  */
  enum class APDMode {SEARCH, AMPLITUDE, END};

  /*
    Push membrane voltage to algorithm
    -- voltage: membrane voltage at a single time point (mV)
  */
  APDMode push_voltage(double voltage);

  /*
    Set time increment of data
    -- dt: time increment (ms)
  */
  void set_dt(double dt);

  /*
    Reset algorithm to SEARCH state
  */
  void reset();

  /*
    Return current mode of APD algorithm
  */
  APDMode get_mode();

  /*
    Return most recent action potential duration calculated
  */
  double get_apd();

 private:
  APDMode mode = APDMode::SEARCH; // default start in SEARCH mode
  double dt; // (ms)
  double apd; // (ms)
  double voltage; // (mV)
  double upstrokeThreshold = -40.0; // (mV)
  int repolarizationPercent = 90; // (%)
  double peakVoltage; // (mV)
  double downstrokeThreshold; // (mV)
  double apAmplitude; // (mV)
  int steps; // Number of steps from start to end of action potential
};

#endif // CALCULATE_APD_HPP
