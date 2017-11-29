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
 * LivRudy2009.h, v4.0
 *
 * Author: Francis Ortega (v1.0-4.0)(2011 - 2017)
 *
 *** NOTES
 *
 * Model equations taken from Matlab model created by Dr. Eric Sobie based
 * on Livshitz 2009 paper, which was used on Devenyi and Ortega 2016 dynamic
 * clamp paper. More recently, the model has been revised to more closely match
 * original MATLAB model from Livshitz and Rudy. Optimizations and modifications
 * have been split off into different git branches. Extensive commenting has
 * been added that is missing from the original.
 *
 * v1.0 - Initial version
 * v1.1 - Added rate changer
 * v1.2 - Replaced calls to PowFast library with calls to RealTimeMath library
 * v2.0 - Made compatible with and added to Model Cell Library (MCL)
 * v2.1 - Added parameter functions for use with gScaling DynClamp
 * v3.0 - Pulled out of MCL to be used in GA fitting of DynClamp scalings
 * v4.0 - Split off all computational optimizations into its own git branch
 *      - Implemented original Livshitz Rudy calcium buffering equations
 *      - Reorganized and reviewed comments and equations
 *      - Model now completely matches original, now defunct, MATLAB code
 *
 ***/

// Header guard
#ifndef LIVRUDY2009_H
#define LIVRUDY2009_H

#include <vector>
#include <complex> // Calcium buffering involves use of complex numbers

class LivRudy2009 {

 public:
  LivRudy2009(void);
  ~LivRudy2009(void);

  /*
    Current clamp model
    Returns 1 if model passed crash tests, 0 if it crashed
    -- current: current injected into model (A)
  */
  int iClamp(double current);

  /*
    Voltage clamp model
     Returns 1 if model passed crash tests, 0 if it crashed
    -- voltage: clamped model voltage (mV)
  */
  int vClamp(double voltage);

  /*
    Reset model to initial conditions
  */
  void reset();

  /*** Set functions ***/
  /*
    Set model dt
    -- dt: (ms)
  */
  void setDt(double new_dt) { DT = new_dt; };

  /*
    Set conductances
  */
  void setGNa(double new_g) { GNa_ = new_g; };
  void setGNab(double new_g) { GNab = new_g; };
  void setGCaL(double new_g) { GCaL_ = new_g; };
  void setGCaT(double new_g) { GCaT_ = new_g; };
  void setGCab(double new_g) {GCab_ = new_g; };
  void setGK1(double new_g) { GK1_ = new_g; };
  void setGKr(double new_g) { GKr_ = new_g; };
  void setGKs(double new_g) { GKs_ = new_g; };
  void setGKp(double new_g) { GKp_ = new_g; };
  void setGNaK(double new_g) { INaK_ = new_g; };
  void setGNCX(double new_g) { GNCX_ = new_g; };
  void setGJrel(double new_g) { GJrel_ = new_g; };
  void setGpCa(double new_g) { IpCa_ = new_g; };
  void setGserca(double new_g) { Jserca_ = new_g; };

  /*
    Set concentrations
  */
  void setNai(double new_c) { Nai = new_c; };
  void setKi(double new_c) { Ki = new_c; };
  void setCai_t(double new_c) { Cai = new_c; };
  void setCaNSR(double new_c) { CaNSR = new_c; };
  void setCaJSR_t(double new_c) { CaJSR = new_c; };

  /*
    Set conditions
    -- conditions: array of model state values
  */
  void setConditions(std::vector<double> &conditions);

  /*** Get functions ***/
  /*
    Get model dt
    -- dt: (ms)
  */
  double getDt() const { return DT; };

  /*
    Check crash status
  */
  const int getStatus();

  /*
    Get conductances
  */
  double getGNa() { return GNa_; };
  double getGNab() { return GNab; };
  double getGCaL() { return GCaL_; };
  double getGCaT() { return GCaT_; };
  double getGCab() { return GCab_; };
  double getGK1() { return GK1_; };
  double getGKr() { return GKr_; };
  double getGKs() { return GKs_; };
  double getGKp() { return GKp_; };
  double getGNaK() { return INaK_; };
  double getGNCX() { return GNCX_; };
  double getGJrel() { return GJrel_; };
  double getGpCa() { return IpCa_; };
  double getGserca() { return Jserca_; };
  double getVm() { return V; };
  double getI() { return Iion; };

  /*
    Get conditions
    Returns double array of model state values
  */
  std::vector<double> getConditions();

  /*
    Get current values
  */
  double getITotal() { return Iion; };
  double getINa() { return INa; };
  double getINab() { return INab; };
  double getICaL() { return ICaL; };
  double getICaL_Na() { return ICaL_Na; };
  double getICaL_K() { return ICaL_K; };
  double getICab() { return ICab; };
  double getICaT() { return ICaT; };
  double getIpCa() { return IpCa; };
  double getIKr() { return IKr; };
  double getIKs() { return IKs; };
  double getIK1() { return IK1; };
  double getIKp() { return IKp; };
  double getINCX() { return INCX; };
  double getINaK() { return INaK; };
  double getI_Inject() { return I_Inject; };

  /*
    Get internal concentrations
  */
  double getNai() { return Nai; };
  double getKi() { return Ki; };
  double getCai() { return Cai; };
  double getCai_t() { return Cai; };
  double getCaJSR() { return CaJSR; };
  double getCaJSR_t() { return CaJSR; };
  double getCaNSR() { return CaNSR; };

 private:
  // Model solver
  void solve();
  // Calcium buffering approximation
  double calcium_buffer(double, double, double, double, double);

  double DT; // Model time-step (ms)

  // Physical constants
  const double F = 96485.0; // Faraday's constant (C/mol)
  const double R = 8314.0; // Gas constant (mJ/K)
  const double T = 310.0; // Absolute temperature at 37C (K)
  const double RTF = R * T / F;
  const double FRT = 1 / RTF;
  const double pi = 4.0 * atan(1.0);

  // Cell Geometry
  const double length_cell = 0.01; // Length of the cell (cm)
  const double radius = 0.0011; // Radius of the cell (cm)
  const double Vcell = 1000 * pi * radius * radius *
      length_cell; // Cell volume (uL)
  const double Ageo = 2 * pi * radius * radius + 2 * pi * radius *
      length_cell; // Geometric membrane area (cm^2)
  const double Acap = 2 * Ageo; // Capacitive membrane area (cm^2)
  double Vmyo = Vcell * 0.68; // Myoplasm volume (uL)
  double Vmito = Vcell * 0.24; // Mitochondria volume (uL)
  double VNSR = Vcell * 0.0552;; // NSR volume (uL)
  double VJSR = Vcell * 0.0048; // JSR volume (uL)

  // Cell Capacitance, implied 1 uF/cm^2
  // const double Cm = 1.0;

  // Voltage
  double V; // Membrane voltage (mV)
  double dVdt; // Change in voltage / Change in time (mV/ms)

  // Static extracellular ionic concentrations
  double Ko; // Extracellular K concentration (mM)
  double Nao; // Extracellular Na concentration (mM)
  double Cao; // Extracellular Ca concentration (mM)

  // Myoplasmic ionic concentrations
  double Nai; // Intracellular Na concentration (mM)
  double Ki; // Intracellular K concentration (mM)
  double Cai_t; // Total buffered and free intracellular Ca concentration (mM)
  double Cai; // Buffered myoplasmic Ca concentration (mM)
  double CaNSR; // NSR Ca concentation (mM)
  double CaJSR; // Buffered JSR Ca concentration (mM)
  double CaJSR_t; // Total JSR buffered and free Ca concentration (mM)

  // Myoplasmic ionic concentration changes
  double dNai; // Change in intracellular Na (mM)
  double dKi; // Change in intracellular K (mM)
  double dCai_t; // Change in total intracellular Ca (mM)

  // JRS Ca concentration changes
  double dCaJSR_t; // Change in total JSR Ca concentration (mM)
  double Jrel; // JSR Ca release to myoplasm (mM/ms)
  double GJrel_; // Jrel scaling parameter, nominal value = 1.0
  double Jrelinf; // Steady-state value of JSR Ca release
  double tau_rel; // JSR Ca release time constant
  double Krel_inf; // Half-saturation coefficient of steady-state JSR Ca release
  double Krel_tau; // Tau of the half-saturation coefficient
  double hrel; // JSR Ca release hill coefficient
  double alpha_rel; // JSR Ca release mplitude coefficient
  double beta_tau; // Maximal value of JSR Ca release time constant
  double dJrel; // Change in JSR Ca release to myoplasm (mM)

  // NSR Ca concentration changes
  double dCaNSR; // Change in NSR Ca concentration (mM)
  double Jserca; // Ca uptake from myoplasm to NSR due to SERCA (mM/ms)
  double Jserca_; // Jserca scaling factor, nominal = 1.0
  double Vserca; // Maximal current through SERCA (mM/ms)
  double Kmserca; // Half-saturation concentration of SERCA (mM)
  double CaNSR_max; // Maximal Ca concentration in NSR (mM)
  double Jtr; // NSR Ca transfer to JSR (mM/ms)
  double Jtr_tau; // Time constant of NSR Ca transfer to JSR (ms)

  // Buffers in cytosol
  double TRPNtot; // Maximal Ca buffered in troponin (mM)
  double KmTRPN; // Equilibrium constant of buffering for troponin (mM)
  double CMDNtot; // Maximal Ca buffered in calmodulin (mM)
  double KmCMDN; // Equilibrium connstant of buffering for calmodulin (mM)

  // Calcium buffering analytical computation parameters
  double alp0, alp1, alp2;
  double q, r;
  std::complex<double> qr, root_qr, cuberoot_rqr, t;

  // Buffers in JSR
  double CSQNtot; // Maximal Ca buffered in calsequestrin (mM)
  double KmCSQN; // Equilibrium connstant of buffering for calsequestrin (mM)

  // Na currents
  double ENa; // Na reversal potential (mV)

  // Fast Na channel current
  double INa; // Fast Na current (uA/uF)
  double GNa_; // Fast Na current conductance (mS/uF)
  double m; // Fast Na current activation gate
  double minf; // Fast Na current activation gate steady-state value
  double taum; // Fast Na current activation gate time constant (1/ms)
  double h; // Fast Na current inactivation gate
  double hinf; // Fast Na current inactivation gate steady-state value
  double tauh; // Fast Na current inactivation gate time constant (1/ms)
  double j; // Fast Na current slow inactivation gate
  double jinf; // Fast Na current slow inactivation gate steady-state value
  double tauj; // Fast Na current slow inactivation gate time constant (1/ms)
  double lambda_na; // Auxiliary function to remove singularities
  double am; // Fast Na current alpha-m rate constant (1/ms)
  double bm; // Fast Na current beta-m rate constant (1/ms)
  double ah; // Fast Na current alpha-h rate constant (1/ms)
  double bh; // Fast Na current beta-h rate constant (1/ms)
  double aj; // Fast Na current alpha-j rate constant (1/ms)
  double bj; // Fast Na current beta-j rate constant (1/ms)

  // Time-independent background Na current
  double INab; // Background Na current (uA/uF)
  double GNab; // Background Na conductance (mS/uF)

  // K currents
  double EK; // K reversal potential (mV)

  // Time-independent K current
  double IK1; // Time-independent K current (uA/uF)
  double GK1_; // Time-independent K current conductance (mS/uF)
  double xK1; // Time-independent K current steady-state value

  // Fast component of delayed rectifier K current
  double IKr; // Rapidly activating K current (uA/uF)
  double GKr_; // Rapidly activating K current conductance (mS/uF)
  double xKr; // Rapidly activating K current activation gate
  double xKrinf; // Rapidly activating K activation gate steady-state value
  double tauxKr; // Rapidly activating K activation gate time constant (1/ms)
  double RKr; // Rapidly activating K current inactivation gate

  // Slow component of delayed rectifier current
  double IKs; // Slowly activating K current (uA/uF)
  double EKs; // Slowly activating K current reversal potential (mV)
  double GKs_; // Slowly activating K current conductance (mS/uF)
  double xs1; // Slowly activating K current fast activation gate
  double tauxs1; // Fast activation gate time constant (1/ms)
  double xs2; // Slowly activating K current slow activation gate
  double tauxs2; // Slow activation gate time constant (1/ms)
  double xsinf; // Slowly activating K current activation steady-state value
  double pKNa; // Slowly activating K current Na to K permeability ratio

  // Plateau K current
  double IKp; // Plateau K current (uA/uF)
  double GKp_; // Plateau K current conductance (mS/uF)
  double Kp; // Plateau K current factor

  // Ca currents
  double ECa; // Ca reversal potential (mV)

  // L-type Ca channel current
  double ICaL; // Ca current through L-type Ca channel (uA/uF)
  double ICa_; // Ca maximal current through L-type Ca channel (uA/uF)
  double ICaL_K; // K current through L-type Ca channel (uA/uF)
  double ICaK_; // K maximal current through L-type Ca channel (uA/uF)
  double ICaL_Na; // Na current through L-type Ca channel (uA/uF)
  double ICaNa_; // Na maximal current through L-type Ca channel (uA/uF)
  double GCaL_; // L-type Ca current scaling factor, nominal = 1.0
  double d; // L-type Ca current V-dependent activation gate
  double dinf; // L-type V-dependent activation gate steady-state value
  double taud; // L-type V-dependent activation gate time connstant (1/ms)
  double f; // L-type Ca current V-dependent inactivation gate
  double finf; // L-type V-dependent inactivation gate steady-state value
  double tauf; // L-type V-dependent inactivation gate time constant (1/ms)
  double fCa; // L-type Ca current Ca-dependent inactivation gate
  double KmCa; // Ca-dependent inactivation gate half-saturation constant (mM)
  double PCa; // Ca membrane permeability (cm/s)
  double PCa_Na; // Na membrane permeability (cm/s)
  double PCa_K; // K membrane permeability (cm/s)
  double gamma_Cao; // L-type Ca activity coefficient of extracellular Ca
  double gamma_Cai; // L-type Ca activity coefficient of intracellular Ca
  double gamma_Nao; // L-type Ca aActivity coefficient of extracellular Na
  double gamma_Nai; // L-type Ca activity coefficient of intracellular Na
  double gamma_Ko; // L-type Ca activity coefficient of extracellular K
  double gamma_Ki; // L-type Ca activity coefficient of intracellular K

  // T-type Ca channel current
  double ICaT; // T-type Ca current (uA/uF)
  double GCaT_; // T-type Ca current conductance (mS/uF)
  double b; // T-type Ca current activation gate
  double binf; // T-type Ca current activation gate steady-state value
  double taub; // T-type Ca current activation gate time constant (1/ms)
  double g; // T-type Ca current inactivation gate
  double ginf; // T-type Ca current inactivation gate steady-state value
  double taug; // T-type Ca current inactivation gate time constant (1/ms)
  double lambda_g; // Auxiliary function to remove singularities

  // Time-independent background Ca current
  double ICab; // Time-independent background Ca current (uA/uF)
  double GCab_; // Time independent background Ca current conductance (mS/uF)

  // Pumps and transporters

  // Na-K Pump
  double INaK; // Na-K pump current (uA/uF)
  double INaK_; // Na-K pump maximal current (uA/uF)
  double fNaK; // Na-K pump voltage dependence parameter
  double sigma_NaK; // Na-K pump extracellular Na dependence factor
  double KmNa_NaK; // Na-K pump half-saturation concentration of Na (mM)
  double KmK_NaK; // Na-K pump half-saturation concentration of K (mM)

  // Na-Ca exchanger
  double INCX; // Na-Ca exchanger current (uA/uF)
  double GNCX_; // Na-Ca exchanger scaling parameter, nominal value = 1.0
  double kNCX; // Na-Ca exchanger scaling factor (uA/uF)
  double ksat; // Na-Ca exchanger half-saturation concentration (mM)
  double eta; // Position of energy barrier controlling voltage dependence

  // Sarcolemmal Ca pump
  double IpCa; // Sarcolemmal Ca pump current (uA/uF)
  double IpCa_; // Sarcolemmal Ca pump maximal current (uA/uF)
  double KmpCa; // Sarcolemmal Ca pump half-saturation concentration (mM)

  // Summated ionic currents
  double I_Inject; // Injected current (uA/uF)
  double NaIon, KIon, CaIon; // Currents flow per ion (uA/uF)
  double Iion; // Total membrane current (uA/uF)
};

#endif
