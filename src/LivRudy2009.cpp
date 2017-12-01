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
 * LivRudy2009.cpp, v4.0
 *
 * Author: Francis Ortega (v1.0 - 4.0) (2011 - 2017)
 *
 *** NOTES
 *
 * Notes in Header
 *
 ***/

#include <cmath>
#include <iostream>

#include "../include/LivRudy2009.hpp"

LivRudy2009::LivRudy2009(void) { // Model initialization
  DT = 0.1; // Model time-step (ms)

  reset(); // Set to initial conditions

  // Model constants

  // Static extracellular ionic concentrations
  // Original Livshitz Rudy concentrations
  // Ko = 4.5; // Extracellular K concentration (mM)
  // Nao = 140.0; // mM
  // Cao = 1.8; // mM
  // Concentrations used in experiments
  Ko = 5.4 ; // Extracellular K concentration (mM)
  Nao = 137.0; // Extracellular Na concentration (mM)
  Cao = 2.0; // Extracellular Ca concentration (mM)

  // JSR constants
  GJrel_ = 1.0; // Jrel scaling parameter, nominal value = 1.0
  alpha_rel = 0.125; // JSR Ca release amplitude coefficient
  Krel_inf = 1; // Half-saturation coefficient of steady-state JSR Ca release
  Krel_tau = 0.0123; // Tau of the half-saturation coefficient
  hrel = 9; // JSR Ca release hill coefficient
  beta_tau = 4.75; // Maximal value of JSR Ca release time constant

  // NSR constants
  Jserca_ = 1.0; // Jserca scaling factor, nominal = 1.0
  Vserca = 8.75e-3; // Maximal current through SERCA (mM/ms)
  Kmserca = 9.0e-4; // Half-saturation concentration of SERCA (mM)
  CaNSR_max = 15.0; // Maximal Ca concentration in NSR (mM)
  Jtr_tau = 120; // Time constant of NSR Ca transfer to JSR (ms)

  // Myoplasmic Ca buffering constants
  TRPNtot = 70e-3; // Maximal Ca buffered in troponin (mM)
  KmTRPN = 0.5e-3; // Equilibrium constant of buffering for troponin (mM)
  CMDNtot = 50e-3; // Maximal Ca buffered in calmodulin (mM)
  KmCMDN = 2.38e-3; // Equilibrium connstant of buffering for calmodulin (mM)

  // JSR Ca buffering constants
  CSQNtot = 10; // Maximal Ca buffered in calsequestrin (mM)
  KmCSQN = 0.8; // Equilibrium connstant of buffering for calsequestrin (mM)

  // Na current constants
  GNa_= 16; // Fast Na current conductance (mS/uF)
  GNab = 0.004; // Background Na current conductance (mS/uF)

  // K current constants
  GK1_ = 0.75; // Time-independent K current conductance (mS/uF)
  GKr_ = 0.02614; // Rapidly activating K current conductance (mS/uF)
  GKs_ = 0.433; // Slowly activating K current conductance (mS/uF)
  pKNa = 0.01833; // Slowly activating K current Na to K permeability ratio
  GKp_ = 5.52e-3; // Plateau K current conductance (mS/uF)

  // Ca current consants
  // L-type Ca current
  GCaL_ = 1.0; // L-type Ca current scaling factor, nominal = 1.0
  PCa = 5.4e-4; // Ca membrane permeability (cm/s)
  PCa_Na = 6.75e-7; // Na membrane permeability (cm/s)
  PCa_K = 1.93e-7; // K membrane permeability (cm/s)
  gamma_Cao = 0.341; // L-type Ca activity coefficient of extracellular Ca
  gamma_Cai = 1; // L-type Ca activity coefficient of intracellular Ca
  gamma_Nao = 0.75; // L-type Ca aActivity coefficient of extracellular Na
  gamma_Nai = 0.75; // L-type Ca activity coefficient of intracellular Na
  gamma_Ko = 0.75; // L-type Ca activity coefficient of extracellular K
  gamma_Ki = 0.75; // L-type Ca activity coefficient of intracellular K
  KmCa = 6e-4; // Ca-dependent inactivation gate half-saturation constant (mM)
  // T-type & background current
  GCaT_ = 0.05; // T-type Ca current conductance (mS/uF)
  GCab_ = 0.003016; // Time independent background Ca conductance (mS/uF)

  // Pumps and transporters
  // Na-K pump
  INaK_ = 2.25; // Na-K pump maximal current (uA/uF)
  KmNa_NaK = 10; // Na-K pump half-saturation concentration of Na (mM)
  KmK_NaK = 1.5; // Na-K pump half-saturation concentration of K (mM)
  // Na-Ca exchanger
  GNCX_ = 1.0; // Na-Ca exchanger scaling parameter, nominal value = 1.0
  kNCX = 0.00025; // Na-Ca exchanger scaling factor (uA/uF)
  ksat = 0.0001; // Na-Ca exchanger half-saturation concentration (mM)
  eta = 0.15; // Position of energy barrier controlling voltage dependence
  // Sarcolemmal Ca pump
  IpCa_ = 1.15; // Sarcolemmal Ca pump maximal current (uA/uF)
  KmpCa = 5e-4; // Sarcolemmal Ca pump half-saturation concentration (mM)
}

LivRudy2009::~LivRudy2009(void){
}


// Model Solver
void LivRudy2009::solve(){
  // Calculate free Ca through rapid buffer approximation
  Cai = calcium_buffer(Cai_t, TRPNtot, KmTRPN, CMDNtot, KmCMDN);
  CaJSR = calcium_buffer(CaJSR_t, CSQNtot, KmCSQN, 0, 0);

  // Reversal potentials
  ENa = RTF * log(Nao / Nai); // Na reversal potential (mV)
  EK = RTF * log(Ko / Ki); // K reversal potential (mV)
  ECa = 0.5 * RTF * log(Cao / Cai); // Ca reversal potential (mV)

  // Fast Na channel current
  // Inactivation gate
  // Auxiliary function to remove singularities
  lambda_na = 1.0 - 1.0 / (1.0 + exp(-(V + 40) / 0.024));
  // Alpha-h rate constant (1/ms)
  ah = lambda_na * 0.135 * exp(-(80.0 + V) / 6.8);
  // Beta-h rate constant (1/ms)
  bh = (1.0 - lambda_na) / (0.13 * (1.0 + exp((V + 10.66)/(-11.1)))) +
      lambda_na * (3.56 * exp(0.079 * V) + 3.1 * 1e5 * exp(0.35 * V));
  hinf = ah / (ah + bh); // Inactivation gate steady-state value
  tauh = 1 / (ah + bh); // Inactivation gate time constant (1/ms)
  // Slow inactivation gate
  // Alpha-j rate constant (1/ms)
  aj =  lambda_na *
      (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) *
      (V + 37.78) / (1 + exp(0.311 * (V + 79.23)));
  // Beta-j rate constant (1/ms)
  bj = (1 - lambda_na) *
      (0.3 * exp(-2.535e-7 * V) / (1 + exp(-0.1 * (V + 32)))) + lambda_na *
      (0.1212 * exp(-0.01052 * V) / (1 + exp(-0.1378 * (V + 40.14))));
  tauj = 1.0 / (aj + bj); // Slow inactivation gate time constant (1/ms)
  jinf = aj / (aj + bj); // Slow inactivation gate steady-state value
  // Activation gate
  // Fast Na current alpha-m rate constant (1/ms)
  if (V > -47.14 && V < -47.12) // if V = -47.13, divide by 0 error
    am = 3.2;
  else
    am = 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
  bm = 0.08 * exp(-V / 11.0); // Beta-m rate constant (1/ms)
  minf = am / (am + bm); // Activation gate steady-state value
  taum = 1.0 / (am + bm); // Activation gate time constant (1/ms)
  INa = GNa_ * m * m * m * h * j * (V - ENa); // Fast Na current (uA/uF)

  // Time-independent background Na current
  INab = GNab * (V - ENa); // Background Na current (uA/uF)

  // K currents

  // Time-independent K current
  // Steady-state value
  xK1 = 0.004 * (1.0 + exp(0.6987 * (V - EK + 11.724))) /
      (1.0 + exp(0.6168 * (V - EK + 4.872)));
  // Time-independent K current (uA/uF)
  IK1 = GK1_ * sqrt(Ko / 5.4) * (V - EK) / (1.0 + xK1);

  // Fast component of the delayed rectifier K current
  // Activation gate steady-state value
  xKrinf = 1.0 / (1.0 + exp(-(V + 21.5) / 7.5));
  // Activation gate time constant (1/ms)
  if (V > -14.21 && V < -14.19) // if V = -14.2, divide by 0 error
    tauxKr = 85.830287334611480;
  else
    tauxKr = (1.0 / (0.00138 * (V + 14.2) /
                     (1.0 - exp(-0.123 * (V + 14.2))) + 0.00061 *
                     (V + 38.9) / (exp(0.145 *(V + 38.9)) -1.0)));
  RKr = 1.0 / (exp((V + 9.0) / 22.4) + 1.0); // Inactivation gate
  // Rapidly activating K current (uA/uF)
  IKr = GKr_ * sqrt(Ko / 5.4) * xKr * RKr * (V - EK);

  // Slow component of the delayed rectifier K current
  // Reversal potential (mV)
  EKs = RTF * log((Ko + pKNa * Nao)/(Ki + pKNa * Nai));
  // Fast activation gate time constant (1/ms)
  if (V > -30.01 && V < -29.99) // if V = -30, divide by 0 error
    tauxs1 = 417.9441667499822;
  else
    tauxs1 = 10000.0 / (0.719 * (V + 30.0) / (1 - exp(-0.148 * (V + 30.0))) +
                        1.31 * (V + 30.0) / (exp(0.0687 * (V + 30.0)) - 1.0));
  tauxs2 = 4.0 * tauxs1; // Slow activation gate time constant (1/ms)
  // Activation steady-state value
  xsinf = 1.0 / (1.0 + exp(-(V - 1.5) / 16.7));
  // Slowly activating K current (uA/uF)
  IKs = GKs_ * (1.0 + 0.6 / (exp(1.4 * log(3.8e-5 / Cai)) + 1.0)) * xs1 *
      xs2 * (V - EKs); // pow() removed

  // Plateau K current
  Kp = 1.0 / (1.0 + exp((7.488 - V) / 5.98)); // Plateau K current factor
  IKp = GKp_ * Kp * (V - EK); // Plateau K current (uA/uF)

  // Ca currents

  // L-type Ca channel current
  // V-dependent activation gate
  // V-dependent activation gate steady-state value
  dinf = (1.0 / (1.0 + exp(-(V + 10) / 6.24))) *
      (1.0 / (1.0 + exp(-(V + 60) / 0.024)));
  // V-dependent activation gate time connstant (1/ms)
  if (V > -10.01 && V < -9.99) // if V = -10, divide by 0 error
    taud = 2.289374849326888;
  else
    taud =  1.0 / (1.0 + exp(-(V + 10) / 6.24)) *
        (1 - exp(-(V + 10) / 6.24))/(0.035 * (V + 10));
  // V-dependent inactivation gate
  // V-dependent inactivation gate steady-state value
  finf = 1.0 / (1.0 + exp((V + 32) / 8.0)) +
      (0.6) / (1.0 + exp((50 - V) / 20.0));
  // V-dependent inactivation gate time constant (1/ms)
  tauf = 1.0 /
      (0.0197 * exp(-(0.0337 * (V + 10)) *
                    (0.0337 * (V + 10))) + 0.02);
  // Ca maximal current through L-type Ca channel (uA/uF)
  ICa_ = PCa * 4.0 * F * FRT * V * (gamma_Cai * Cai * exp(2.0 * V *FRT) -
                                    gamma_Cao *Cao) / (exp(2.0 * V *FRT) - 1.0);
  // K maximal current through L-type Ca channel (uA/uF)
  ICaK_ = PCa_K * F * FRT * V * (gamma_Ki *Ki * exp(V * FRT) -
                                 gamma_Ko * Ko) / (exp(V * FRT) - 1.0);
  // Na maximal current through L-type Ca channel (uA/uF)
  ICaNa_ = PCa_Na * F * FRT * V * (gamma_Nai * Nai * exp(V * FRT) -
                                   gamma_Nao * Nao) / (exp(V * FRT) - 1.0) ;
  fCa = 1.0 / (Cai / KmCa + 1.0); // Ca-dependent inactivation gate
  // Ca current through L-type Ca channel (uA/uF)
  ICaL = GCaL_ * ICa_ * d * f * fCa;
  // K current through L-type Ca channel (uA/uF)
  ICaL_K = GCaL_ * ICaK_* d * f * fCa;
  // Na current through L-type Ca channel (uA/uF)
  ICaL_Na = GCaL_ * ICaNa_ * d * f * fCa;

  // T-type Ca current
  // Activation gate
  // Activation gate steady-state value
  binf = 1.0 / (1.0 + exp(-(V + 14.0) / 10.8));
  // Activation gate time constant (1/ms)
  taub = (3.7 + 6.1 / (1 + exp((V + 25.0) / 4.5)));
  // Inactivation gate
  // Auxiliary function to remove singularities
  lambda_g = 1.0 - 1.0 / (1.0 + exp(-V / 0.0024));
  // Inactivation gate steady-state value
  ginf = 1.0 / (1.0 + exp((V + 60.0) / 5.6));
  // Inactivation gate time constant (1/ms)
  taug = (lambda_g * (-0.875 * V + 12.0) + 12.0 * (1.0 - lambda_g));
  ICaT = GCaT_ * b*b * g * (V - ECa); // T-type Ca current (uA/uF)

  // Time-independent background Ca current
  ICab = GCab_ * (V - ECa); // Time-independent background Ca current (uA/uF)

  // Pumps and transporters

  // Na-K pump
  sigma_NaK = (exp(Nao / 67.3) - 1) / 7.0; // Extracellular Na dependence factor
  fNaK = 1.0/(1.0 + 0.1245 * exp(-0.1 * V * FRT) + 0.0365 * sigma_NaK *
              exp(-V * FRT)); //Voltage dependence parameter
  // Na-K pump current (uA/uF)
  INaK = INaK_ * fNaK * Ko / ((Ko + KmK_NaK) *
                              (1.0 + ((KmNa_NaK / Nai) *
                                      (KmNa_NaK / Nai)))); // pow() removed

  // Na-Ca exchanger
  INCX = GNCX_ * kNCX * exp((eta - 1.0) * V * FRT) *
      ((Nai * Nai * Nai) * Cao * exp(V * FRT) - (Nao * Nao * Nao) * Cai) /
      ( 1.0 + ksat * exp((eta - 1.0) * V * FRT) *
        ((Nai * Nai * Nai) * Cao * exp(V * FRT) +
         (Nao * Nao * Nao) * Cai)); // Na-Ca exchanger current (uA/uF)

  // Sarcolemmal calcium pump
  IpCa = IpCa_ * Cai / (Cai + KmpCa); // Sarcolemmal Ca pump current (uA/uF)

  // Intracellular Ca fluxes

  // JSR Ca compartment
  // Steady-state value of JSR Ca release - Added GJrel_ as a scaling factor and
  // removed pow()
  Jrelinf = GJrel_ * alpha_rel * beta_tau * ICaL /
      (exp(hrel * log(Krel_inf / CaJSR)) + 1);
  tau_rel = beta_tau / (Krel_tau / CaJSR + 1); // JSR Ca release time constant
  // Change in JSR release to myoplasm (mM)
  dJrel = - (Jrelinf + Jrel) / tau_rel;

  // NSR Ca compartment
  // Ca uptake from myoplasm to NSR due to SERCA (mM/ms)
  Jserca = Jserca_ * Vserca * (Cai / (Cai + Kmserca) - CaNSR / CaNSR_max);
  // NSR Ca transfer to JSR (mM/ms)
  Jtr = (CaNSR - CaJSR) / Jtr_tau;

  // Summated ionic currents
  NaIon = INa + INab + 3 * INCX + ICaL_Na + 3 * INaK;
  // Injected current assumed to be due to potassium flux
  KIon = I_Inject + IKr + IKs + IK1 + IKp + ICaL_K - 2 * INaK;
  CaIon = ICaL + ICab + IpCa - 2 * INCX + ICaT;
  Iion = NaIon + KIon + CaIon;

  // Units are in mM. Note that summated currents must be multipled by cell
  // capacitance to get the correct units. Since cell capacitance = 1 uF/cm^2,
  // it doesn't explicitly appear in the equation below.
  // Derivatives for ionic concentration
  dNai = -NaIon * Acap / (Vmyo * F);
  dKi = -KIon * Acap / (Vmyo * F);
  dCai_t = -Jserca * VNSR / Vmyo + Jrel * VJSR / Vmyo -
      CaIon * Acap / (2 * Vmyo * F);
  dCaJSR_t = Jtr - Jrel;
  dCaNSR = Jserca - Jtr * VJSR / VNSR;

  // Derivative for voltage
  dVdt = -(Iion);

  // Update voltage and ionic concentrations
  V += DT * dVdt; // Membrane voltage (mV)
  Nai += DT * dNai; // Intracellular Na concentration (mM)
  Ki += DT * dKi; // Intracellular K concentration (mM)
  // Total buffered and free intracellular Ca concentration (mM)
  Cai_t += DT * dCai_t;
  // Total buffered and free JSR Ca concentration (mM)
  CaJSR_t += DT * dCaJSR_t;
  CaNSR += DT * dCaNSR;
  Jrel += DT * dJrel;

  // Update gating variables - Euler Method
  // Fast Na channel gates
  h = (hinf - (hinf - h) * exp(-DT / tauh));
  j = (jinf - (jinf - j) * exp(-DT / tauj));
  m = (minf - (minf - m) * exp(-DT / taum));
  // Fast component of the delayed rectifier K channel gates
  xKr = (xKrinf - (xKrinf - xKr) * exp(-DT / tauxKr));
  // Slow component of the delayed rectifier K channel gates
  xs1 = (xsinf - (xsinf - xs1) * exp(-DT / tauxs1));
  xs2 = (xsinf - (xsinf - xs2) * exp(-DT / tauxs2));
  // L-type Ca channel gates
  d = (dinf - (dinf - d) * exp(-DT / taud));
  f = (finf - (finf - f) * exp(-DT / tauf));
  // T-type Ca channel gates
  b = (binf - (binf - b) * exp(-DT / taub));
  g = (ginf - (ginf - g) * exp(-DT / taug));
}

// Fast buffer implementation
// a - Maximal Ca buffered
// b - Equilibrium constant of buffering
// 1 - Buffer 1
// 2 - Buffer 2
double LivRudy2009::calcium_buffer(
    double ca_t, double a1, double b1, double a2, double b2) {
  alp2 = a1 + a2 + b1 + b2 - ca_t;
  alp1 = b1 * b2 - ca_t * (b1 + b2) + a1 * b2 + a2 * b1;
  alp0 = -b1 * b2 * ca_t;
  q = (3.0 * alp1 - (alp2 * alp2)) / 9.0;
  r = (9.0 * alp2 * alp1 - 27.0 * alp0 - 2.0 * (alp2 * alp2 * alp2)) / 54.0;
  qr = pow(q,3.0) + pow(r,2.0);
  root_qr = pow(qr, 0.5);
  cuberoot_rqr = pow(r + root_qr, 1.0/3.0);
  t = cuberoot_rqr - q / cuberoot_rqr;

  return abs(t - alp2/3.0);
}

// Voltage Clamp Function
int LivRudy2009::vClamp(double voltage) {
  if (Cai < 0 ||
      Nai < 0 ||
      Ki < 0 ||
      CaJSR < 0 ||
      CaNSR < 0 ||
      V < -200 ||
      V > 200)
    return 0;

  // Clamp model voltage
  V = voltage;

  // Run model solver
  solve(); // voltage free to change during this time period

  // Returns 0 if any of the following conditions are out of bounds
  return 1;
}

// Current Clamp Function
int LivRudy2009::iClamp(double current) {
  if (Cai < 0 ||
      Nai < 0 ||
      Ki < 0 ||
      CaJSR < 0 ||
      CaNSR < 0 ||
      V < -200 ||
      V > 200)
    return 0;

  // Inject current into model
  I_Inject = current;

  // Run model solver
  solve();

  // Returns 0 if any of the following conditions are out of bounds
  return 1;
}

// Crash status
const int LivRudy2009::getStatus() {
  return !(Cai < 0 ||
          Nai < 0 ||
          Ki < 0 ||
          CaJSR < 0 ||
          CaNSR < 0 ||
          V < -200 ||
          V > 200);
}

// Model Reset Function
void LivRudy2009::reset(){ // Reset to initial conditions
  // Initial conditions - 2Hz pacing
  V = -84.6964;
  Cai_t = 0.0218933;
  CaNSR = 2.44962;
  CaJSR_t = 8.17769;
  Nai = 13.919;
  Ki = 137.879;
  m = 0.0016256;
  h = 0.983832;
  j = 0.989807;
  d = 9.88131e-324;
  f = 0.998914;
  b = 0.00143546;
  g = 0.974959;
  xKr = 0.000252518;
  xs1 = 0.0223878;
  xs2 = 0.0726541;
  Jrel = 1.73625e-41;

  Cai = calcium_buffer(Cai_t, TRPNtot, KmTRPN, CMDNtot, KmCMDN);
  CaJSR = calcium_buffer(CaJSR_t, CSQNtot, KmCSQN, 0, 0);
}

// Condition functions
void LivRudy2009::setConditions(std::vector<double> &conditions) {
  if (conditions.size() != 17)
    std::cout << "Error: 17 conditions are required, " << conditions.size() <<
        " conditions were given." << std::endl;

  V = conditions.at(0);
  Cai_t = conditions.at(1);
  CaNSR = conditions.at(2);
  CaJSR_t = conditions.at(3);
  Nai = conditions.at(4);
  Ki = conditions.at(5);
  m = conditions.at(6);
  h = conditions.at(7);
  j = conditions.at(8);
  d = conditions.at(9);
  f = conditions.at(10);
  b = conditions.at(11);
  g = conditions.at(12);
  xKr = conditions.at(13);
  xs1 = conditions.at(14);
  xs2 = conditions.at(15);
  Jrel = conditions.at(16);
}

std::vector<double> LivRudy2009::getConditions() {
  std::vector<double> conditions;

  conditions.push_back(V);
  conditions.push_back(Cai_t);
  conditions.push_back(CaNSR);
  conditions.push_back(CaJSR_t);
  conditions.push_back(Nai);
  conditions.push_back(Ki);
  conditions.push_back(m);
  conditions.push_back(h);
  conditions.push_back(j);
  conditions.push_back(d);
  conditions.push_back(f);
  conditions.push_back(b);
  conditions.push_back(g);
  conditions.push_back(xKr);
  conditions.push_back(xs1);
  conditions.push_back(xs2);
  conditions.push_back(Jrel);

  return conditions;
}
