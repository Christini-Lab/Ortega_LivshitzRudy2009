/*
 * Copyright (C) 2016 Weill Medical College of Cornell University
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
 * LivRudy2009.cpp, v3.0
 *
 * Author: Francis Ortega (v1.0 - 3.0) (2011 - 2016)
 *
 *** NOTES
 *
 * Notes in Header
 *
 ***/

#include <cmath>
#include <iostream>

#include "../include/LivRudy2009.hpp"

#define fastEXP RTMath.fastEXP //define shortcut for fastEXP function
//#define fastEXP exp

LivRudy2009::LivRudy2009(void) { // Model initialization

  // Model parameters
  DT = 0.1;

  reset(); // Set to initial conditions

  // Constants
  // Physical Constants
  F = 96485; // Faraday's constant, C/mol
  R = 8314; // gas constant, mJ/K
  T = 273 + 37; // absolute temperature, K
  RTF = R * T / F;
  FRT = 1 / RTF;
  pi = 4.0 * atan(1.0); // Pi
  // Cell Geometry
  length_cell = 0.01; // Length of the cell (cm)
  radius = 0.0011; // Radius of the cell (cm)
  Vcell = 1000 * pi * radius * radius *
      length_cell; // 3.801e-5 uL Cell volume (uL)
  Ageo = 2 * pi * radius * radius + 2 * pi * radius *
      length_cell; // 7.671e-5 cm^2 Geometric membrane area (cm^2)
  Acap = 2 * Ageo; // 1.534e-4 cm^2 Capacitive membrane area (cm^2)
  Vmyo = Vcell * 0.68; // Myoplasm volume (uL)
  Vmito = Vcell * 0.24; // Mitochondria volume (uL)
  VNSR = Vcell * 0.0552; // NSR volume (uL)
  VJSR = Vcell * 0.0048; // JSR volume (uL)
  Vss = Vcell * 0.02;
  // Cell Capacitance
  Cm = 1.0;
  // Fixed ionic concentrations
  // Original concentrations
  /*
  Ko = 4.5 ; // uM
  Nao = 140 ; // uM
  Cao = 1.8 ; // uM
  */
  // Concentration used in dynamic clamp experiments
  Ko = 5.4 ; // uM
  Nao = 137; // uM
  Cao = 2.0 ; // uM
  // Na current constants
  GNa_= 16; // mS/cm^2
  GNab = 0.004;
  //double GNaL_= 6.5e-3;
  // Ca current constants
  PCa = 5.4e-4; // cm/s
  PCa_Na = 6.75e-7; // cm/s
  PCa_K = 1.93e-7; // cm/s
  PCab = 1.995084e-7; // cm/s
  gamma_Cao = 0.341; // dimensionless
  gamma_Cai = 1; // dimensionless
  gamma_Nao = 0.75; // dimensionless
  gamma_Nai = 0.75; // dimensionless
  gamma_Ko = 0.75; // dimensionless
  gamma_Ki = 0.75; // dimensionless
  GCaL_ = 1.0;
  //const double hLca = 1; // dimensionless, Hill coefficient
  KmCa = 6e-4; // Half saturation constant, mM
  // T-type & background currents
  GCaT_ = 0.05;
  GCab_ = 0.003016;
  // K Currents
  GK1_ = 0.75;
  GKr_ = 0.02614;
  GKs_ = 0.433;
  pKNa = 0.01833; // relative permeability of IKs, Na to K
  GKp_ = 5.52e-3;
  INaK_ = 2.25; // Max. current through Na-K pump (uA/uF)
  KmNa_NaK = 10; // Half-saturation concentration of NaK pump (mM)
  KmK_NaK = 1.5; // Half-saturation concentration of NaK pump (mM)
  ksat = 0.0001;
  eta = 0.15;
  alpha_rel = 0.125;
  Krel_inf = 1;
  hrel = 9;
  beta_tau = 4.75;
  Krel_tau = 0.0123;
  // Pumps and Transporters
  IpCa_ = 1.15; // Max. Ca current through sarcolemmal Ca pump (uA/uF)
  KmpCa = 5e-4; // Half-saturation concentration of sarcolemmal Ca pump (mM)
  GJrel_ = 1.0; // Jrel scaling factor
  Jserca_ = 1.0; // Jserca scaling factor
  Vserca = 8.75e-3; // mM/ms
  Kmserca = 9.0e-4; // mM
  CaNSR_max = 15.0;
  tau_transfer = 120;
  kNCX = 0.00025;
  GNCX_ = 1.0; // Na-Ca exchanger scaling factor
  // Buffers in cytosol
  TRPNtot = 70e-3;
  KmTRPN = 0.5e-3;
  CMDNtot = 50e-3;
  KmCMDN = 2.38e-3;
  // Buffers in JSR
  CSQNtot = 10;
  KmCSQN = 0.8;
}

LivRudy2009::~LivRudy2009(void){
}


// Model Solver
void LivRudy2009::solve(){
  // Reversel Potentials
  ENa = RTF * log(Nao / Nai);
  EK = RTF * log(Ko / Ki);
  EKs = RTF * log((Ko + pKNa * Nao)/(Ki + pKNa * Nai));
  ECa = 0.5 * RTF * log(Cao / Cai);

  // Na currents
  // H-gate
  lambda_na = 1 - 1 / (1 + exp(-(V + 40) / 0.024));
  ah = lambda_na * 0.135 * exp(-(80 + V) / 6.8);
  bh = (1 - lambda_na) / (0.13 * (1 + exp((V + 10.66)/(-11.1)))) +
      lambda_na * (3.56 * exp(0.079 * V) + 3.1 * 1e5 * exp(0.35 * V));
  hinf = ah / (ah + bh); // hinf
  hinf= 1 / (ah + bh); // tauh
  // J-gate
  aj =  lambda_na *
      (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) *
      (V + 37.78) / (1 + exp(0.311 * (V + 79.23)));
  bj = (1 - lambda_na) *
      (0.3 * exp(-2.535e-7 * V) / (1 + exp(-0.1 * (V + 32)))) + lambda_na *
      (0.1212 * exp(-0.01052 * V) / (1 + exp(-0.1378 * (V + 40.14))));
  tauj = 1 / (aj + bj); // tauj
  jinf = aj / (aj + bj); // jinf
  // M-gate
  if (V > -47.14 && V < -47.12) // if V = -47.13, divide by 0 error
    am = 3.199985789461998;
  else
    am = 0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
  bm = 0.08 * exp(-V / 11.0);
  minf = am / (am + bm); // minf
  taum = 1 / (am + bm); // taum
  INa = GNa_ * m * m * m * h * j * (V - ENa);
  INab = GNab * (V - ENa);

  // L-type Ca current
  // D-gate
  dinf_0 = 1 / (1 + exp(-(V + 10) / 6.24));
  dinf_1 = 1 / (1 + exp(-(V + 60) / 0.024));
  dinf = dinf_0 * dinf_1; // dinf
  if (V > -10.01 && V < -9.99)// if V = -10, divide by 0 error
    taud = 2.289374849326888; // taud
  else
    taud =  1 / (1 + exp(-(V + 10) / 6.24)) *
        (1 - exp(-(V + 10) / 6.24))/(0.035 * (V + 10)); // taud
  // F-gate
  finf = 1 / (1 + exp((V + 32) / 8.0)) +
      (0.6) / (1 + exp((50 - V) / 20.0)); // finf
  tauf = 1 /
      (0.0197 * exp(-(0.0337 * (V + 10)) *
                    (0.0337 * (V + 10))) + 0.02); // tauf
  ICa_ = PCa * 4 * F * FRT * V * (gamma_Cai * Cai * fastEXP(2 * V *FRT) -
                                  gamma_Cao *Cao) / (fastEXP(2 * V *FRT) - 1);
  ICaK_ = PCa_K * F * FRT * V * (gamma_Ki *Ki * fastEXP(V * FRT) -
                                 gamma_Ko * Ko) / (fastEXP(V * FRT) - 1);
  ICaNa_ = PCa_Na * F * FRT * V * (gamma_Nai * Nai * fastEXP(V * FRT) -
                                   gamma_Nao * Nao) / (fastEXP(V * FRT) - 1) ;
  fCa = 1 / (Cai / KmCa + 1);
  ICaL = GCaL_ * ICa_ * d * f * fCa;
  ICaL_K = GCaL_ * ICaK_* d * f * fCa;
  ICaL_Na = GCaL_ * ICaNa_ * d * f * fCa;

  // Background calcium current
  ICab = GCab_ * (V - ECa);

  // Sarcolemmal calcium pump
  IpCa = IpCa_ * Cai / (Cai + KmpCa);

  // T-type Ca current
  // B-gate
  binf = 1 / (1 + exp(-(V + 14.0) / 10.8)); // binf
  taub = (3.7 + 6.1 / (1 + exp((V + 25.0) / 4.5))); // taub

  // G-gate
  lambda_g = 1 - 1 / (1 + exp(-V / 0.0024));
  ginf = 1 / (1 + exp((V + 60.0) / 5.6)); // ginf
  taug = (lambda_g * (-0.875 * V+12.0) + 12.0 *
                 (1 - lambda_g)); // taug
  ICaT = GCaT_ * b*b * g * (V - ECa);

  // K currents
  // Time independent K current
  xK1 = 0.004 * (1 + fastEXP(0.6987 * (V - EK + 11.724))) /
      (1 + fastEXP(0.6168 * (V - EK + 4.872)));
  IK1 = GK1_ * sqrt(Ko / 5.4) * (V - EK) / (1 + xK1);

  // Fast component of the delayed rectifier K current
  // IKr
  xKrinf = 1 / (1 + exp(-(V + 21.5) / 7.5)); // xKrinf
  if (V > -14.21 && V < -14.19) // if V = -14.2, divide by 0 error
    tauxKr = 85.830287334611480; // tauxKr
  tauxKr = (1 / (0.00138 * (V + 14.2) /
                      (1 - exp(-0.123 * (V + 14.2))) + 0.00061 *
                      (V + 38.9) / (exp(0.145 *(V + 38.9)) -1))); // tauxKr
  if (V > -30.01 && V < -29.99) // if V = -30, divide by 0 error
    tauxs1 = 417.9441667499822;
  else
    tauxs1 = 10000 / (0.719 * (V + 30) / (1 - exp(-0.148 * (V + 30))) +
                       1.31 * (V + 30) / (exp(0.0687 * (V + 30)) - 1));
  tauxs2 = 4 * tauxs1; // tauxs2
  xsinf = 1 / (1 + exp(-(V - 1.5) / 16.7)); // xsinf
  RKr = 1 / (fastEXP((V + 9) / 22.4) + 1);
  IKr = GKr_ * sqrt(Ko / 5.4) * xKr * RKr * (V - EK);

  // Fast component of the delayed rectifier K current
  //IKs = GKs_ * (1 + 0.6/(pow(3.8e-5/Cai,1.4)+1)) * xs1 * xs2 * (V - EKs);
  IKs = GKs_ * (1 + 0.6 / (fastEXP(1.4 * log(3.8e-5 / Cai)) + 1)) * xs1 *
      xs2 * (V - EKs); // pow() removed

  // Plateau K current
  Kp = 1 / (1 + fastEXP((7.488 - V) / 5.98));
  IKp = GKp_ * Kp * (V - EK);

  // Pumps and transporters
  // Na-K pump
  sigma_NaK = (fastEXP(Nao / 67.3) - 1) / 7.0;
  fNaK = 1/(1 + 0.1245 * fastEXP(-0.1 * V * FRT) + 0.0365 * sigma_NaK *
            fastEXP(-V * FRT));
  //INaK = INaK_ * fNaK * Ko / ( (Ko + KmK_NaK) * pow( 1 +
  //((KmNa_NaK/Nai)*(KmNa_NaK/Nai)),2) );
  INaK = INaK_ * fNaK * Ko / ((Ko + KmK_NaK) *
                              (1 + ((KmNa_NaK / Nai) *
                                    (KmNa_NaK / Nai)))); // pow() removed

  // Na-Ca exchanger
  INCX = GNCX_ * kNCX * fastEXP((eta - 1) * V * FRT) *
      ((Nai * Nai * Nai) * Cao * fastEXP(V * FRT) - (Nao * Nao * Nao) * Cai) /
      ( 1 +ksat * fastEXP((eta - 1) * V * FRT) *
        ((Nai * Nai * Nai) * Cao * fastEXP(V * FRT) +
         (Nao * Nao * Nao) * Cai));

  // Intracellular Ca fluxes
  // SR Ca release, uptake, and leak
  //Jrelinf = alpha_rel * beta_tau * ICaL / (pow((Krel_inf/CaJSR),hrel) + 1);
  // Added GJrel_ as a scaling factor
  Jrelinf = GJrel_ * alpha_rel * beta_tau * ICaL /
      (fastEXP(hrel * log(Krel_inf / CaJSR)) + 1);
  tau_rel = beta_tau / (Krel_tau / CaJSR + 1);
  dJreldt = - (Jrelinf + Jrel) / tau_rel;

  Jserca = Jserca_ * Vserca * (Cai / (Cai + Kmserca) - CaNSR / CaNSR_max);

  Jtr = (CaNSR - CaJSR) / tau_transfer;

  // Buffering factors for rapid buffering approximation
  BJSR = 1.0 / (1 + CSQNtot * KmCSQN / ((KmCSQN + CaJSR) * (KmCSQN + CaJSR)));
  Bi = 1.0 / (1 + (CMDNtot * KmCMDN / ((Cai + KmCMDN) * (Cai + KmCMDN))) +
              (TRPNtot * KmTRPN / ((Cai + KmTRPN) * (Cai + KmTRPN))));

  // Total Current
  Iion = INa + INab + ICaL + ICaL_Na + ICaL_K + ICab + ICaT + IpCa + IKr +
      IKs + IK1 + IKp + INCX + INaK + I_Inject;

  // Derivatives for ionic concentration
  dNai = -(INa + INab + ICaL_Na + 3 * INCX + 3 * INaK) * Acap / (Vmyo * F);
  // Injected current assumed to be due to potassium flux
  dKi = -(I_Inject + IKr + IKs + IK1 + IKp + ICaL_K - 2 * INaK) * Acap /
      (Vmyo * F);
  dCai = Bi * (-Jserca * VNSR / Vmyo + Jrel * VJSR / Vmyo -
               (ICaL + ICaT + ICab + IpCa - 2 * INCX) *
               Acap / (2 * Vmyo * F));
  dCaJSR = BJSR * (Jtr - Jrel);
  dCaNSR = Jserca - Jtr * VJSR / VNSR;

  // Derivative for voltage
  dVdt = -(Iion);

  // Update voltage and ionic concentrations
  V += DT * dVdt;
  Nai += DT * dNai;
  Ki += DT * dKi;
  Cai += DT * dCai;
  CaJSR += DT * dCaJSR;
  CaNSR += DT * dCaNSR;
  Jrel += DT * dJreldt;

  // Update gating variables - Euler Method
  h = (hinf - (hinf - h) *
       fastEXP(-DT / tauh)); // Rush-Larsen approximation used for m gate
  j = (jinf - (jinf - j) *
       fastEXP(-DT / tauj)); // Rush-Larsen approximation used for m gate
  m = (minf - (minf - m) *
       fastEXP(-DT / taum)); // Rush-Larsen approximation used for m gate
  d = (dinf - (dinf - d) *
       fastEXP(-DT / taud));
  f = (finf - (finf - f) *
       fastEXP(-DT / tauf));
  b = (binf - (binf - b) *
       fastEXP(-DT / taub));
  g = (ginf - (ginf - g) *
       fastEXP(-DT / taug));
  xKr = (xKrinf - (xKrinf - xKr) *
       fastEXP(-DT / tauxKr));
  xs1 = (xsinf - (xsinf - xs1) *
       fastEXP(-DT / tauxs1));
  xs2 = (xsinf - (xsinf - xs2) *
       fastEXP(-DT / tauxs2));
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
  // Initial conditions at 0 beats
  /*
  V = -84.7;
  Cai = 0.0822e-3;
  CaNSR = 1.25;
  CaJSR = 1.25;
  Nai = 9.71;
  Ki = 142.82;
  m = 2.46e-4;
  h = 0.99869;
  j = 0.99887;
  d = 1e-4;
  f = 0.983;
  b = 1e-4;
  g = 0.983;
  xKr = 0.229;
  xs1 = 1e-4;
  xs2 = 1e-4;
  Jrel = 1e-4;
  */

  // Initial conditions after 700 beats (Steady State)
  /*
  V = -88.9248;
  Cai = 1.1366e-04;
  CaNSR = 1.7519;
  CaJSR = 1.6213;
  Nai = 14.1538;
  Ki = 136.5910;
  m = 8.0083e-04;
  h = 0.9931;
  j = 0.9957;
  d = 9.8446e-26;
  f = 0.9998;
  b = 9.6988e-04;
  g = 0.9942;
  xKr = 1.2471e-04;
  xs1 = 0.0049;
  xs2 = 0.0243;
  Jrel = 2.4967e-18;
  */

  // Initial conditions after conc change - 800 beats
  // Used in dynamic clamp experiments
  V = -8.185872e+01;
  Cai = 1.797384e-04;
  CaNSR = 2.463960e+00;
  CaJSR = 1.524945e+00;
  Nai = 1.357382e+01;
  Ki = 1.239044e+02;
  m = 2.601169e-03;
  h = 9.697101e-01;
  j = 9.806867e-01;
  d = 5.928788e-322;
  f = 9.981991e-01;
  b = 1.866354e-03;
  g = 9.650771e-01;
  xKr = 4.291283e-04;
  xs1 = 2.954278e-02;
  xs2 = 8.283927e-02;
  Jrel = 1.101473e-38;
}

// Condition functions
void LivRudy2009::setConditions(std::vector<double> &conditions) {
  if (conditions.size() != 17)
    std::cout << "Error: 17 conditions are required, " << conditions.size() <<
        " conditions were given." << std::endl;

  V = conditions.at(0);
  Cai = conditions.at(1);
  CaNSR = conditions.at(2);
  CaJSR = conditions.at(3);
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
  conditions.push_back(Cai);
  conditions.push_back(CaNSR);
  conditions.push_back(CaJSR);
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
