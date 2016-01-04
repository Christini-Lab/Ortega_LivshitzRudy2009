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
 * LivRudy2009.cpp, v3.0
 *
 * Author: Francis Ortega (v1.0 - 3.0) (2011 - 2015)
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

  // Gating Variable Lookup Table
  lkup = new double[20000][20];
  Vx = 0; // Voltage placeholder for lookup table
  V_min = -1000;

  // Lookup Table Initialization
  for (z = 0; z < 20000; z++) {
    Vx = V_min + 0.1 * z;

    // Voltage
    lkup[z][0] = V_min + 0.1 * z;

    // H-gate
    lambda_na = 1 - 1 / (1 + exp(-(Vx + 40) / 0.024));
    ah = lambda_na * 0.135 * exp(-(80 + Vx) / 6.8);
    bh = (1 - lambda_na) / (0.13 * (1 + exp((Vx + 10.66)/(-11.1)))) +
        lambda_na * (3.56 * exp(0.079 * Vx) + 3.1 * 1e5 * exp(0.35 * Vx));
    lkup[z][1] = ah / (ah + bh); // hinf
    lkup[z][2] = 1 / (ah + bh); // tauh

    // J-gate
    aj =  lambda_na *
        (-1.2714e5 * exp(0.2444 * Vx) - 3.474e-5 * exp(-0.04391 * Vx)) *
        (Vx + 37.78) / (1 + exp(0.311 * (Vx + 79.23)));
    bj = (1 - lambda_na) *
        (0.3 * exp(-2.535e-7 * Vx) / (1 + exp(-0.1 * (Vx + 32)))) + lambda_na *
        (0.1212 * exp(-0.01052 * Vx) / (1 + exp(-0.1378 * (Vx + 40.14))));
    lkup[z][3] = 1 / (aj + bj); // tauj
    lkup[z][4] = aj / (aj + bj); // jinf

    // M-gate
    if ( Vx > -47.14 && Vx < -47.12) // if V = -47.13, divide by 0 error
      am = 3.41333;
    else
      am = 0.32 * (Vx + 47.13) / (1 - exp(-0.1 * (Vx + 47.13)));
    bm = 0.08 * exp(-Vx / 11.0);
    lkup[z][5] = am / (am + bm); // minf
    lkup[z][6] = 1 / (am + bm); // taum

    // D-gate
    dinf_0 = 1 / (1 + exp(-(Vx + 10) / 6.24));
    dinf_1 = 1 / (1 + exp(-(Vx + 60) / 0.024));
    lkup[z][7] = dinf_0 * dinf_1; // dinf
    if ( Vx > -10.01 && Vx < -9.99)// if V = -10, divide by 0 error
      lkup[z][8] = 2.30655; // taud
    else
      lkup[z][8] =  1 / (1 + exp(-(Vx + 10) / 6.24)) *
          (1 - exp(-(Vx + 10) / 6.24))/(0.035 * (Vx + 10)); // taud

    // F-gate
    lkup[z][9] = 1 / (1 + exp((Vx + 32) / 8.0)) +
        (0.6) / (1 + exp((50 - Vx) / 20.0)); // finf
    lkup[z][10] = 1 /
        (0.0197 * exp(-(0.0337 * (Vx + 10)) *
                      (0.0337 * (Vx + 10))) + 0.02); // tauf

    // B-gate
    lkup[z][11] = 1 / (1 + exp(-(Vx + 14.0) / 10.8)); // binf
    lkup[z][12] = (3.7 + 6.1 / (1 + exp((Vx + 25.0) / 4.5))); // taub

    // G-gate
    lambda_g = 1 - 1 / (1 + exp(-Vx / 0.0024));
    lkup[z][13] = 1 / (1 + exp((Vx + 60.0) / 5.6)); // ginf
    lkup[z][14] = (lambda_g * (-0.875 * Vx+12.0) + 12.0 *
                   (1 - lambda_g)); // taug

    // IKr
    if ( Vx > -30.01 && Vx < -29.99 ) // if V = -30, divide by 0 error
      tau_xs1 = 411.501;
    else
      tau_xs1 = 10000 / (0.719 * (Vx + 30) / (1 - exp(-0.148 * (Vx + 30))) +
                         1.31 * (Vx + 30) / (exp(0.0687 * (Vx + 30)) - 1));
    lkup[z][15] = 1 / (1 + exp(-(Vx + 21.5) / 7.5)); // xKrinf
    if ( Vx > -14.21 && Vx < -14.19 ) // if V = -14.2, divide by 0 error
      lkup[z][16] = 87.1735; // tauxKr
    lkup[z][16] = (1 / (0.00138 * (Vx + 14.2) /
                        (1 - exp(-0.123 * (Vx + 14.2))) + 0.00061 *
                        (Vx + 38.9) / (exp(0.145 *(Vx + 38.9)) -1))); // tauxKr
    lkup[z][17] = tau_xs1; // tau_xs1
    lkup[z][18] = 4 * tau_xs1; // tau_xs2
    lkup[z][19] = 1 / (1 + exp(-(Vx - 1.5) / 16.7)); // xsinf
  }
}

LivRudy2009::~LivRudy2009(void){
  delete[] lkup;
}


// Model Solver
void LivRudy2009::solve(){
  // Reversel Potentials
  ENa = RTF * log(Nao / Nai);
  EK = RTF * log(Ko / Ki);
  EKs = RTF * log((Ko + pKNa * Nao)/(Ki + pKNa * Nai));
  ECa = 0.5 * RTF * log(Cao / Cai);

  // Na currents
  INa = GNa_ * m * m * m * h * j * (V - ENa);
  INab = GNab * (V - ENa);

  // L-type Ca current
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
  ICaT = GCaT_ * b*b * g * (V - ECa);

  // K currents
  // Time independent K current
  xK1 = 0.004 * (1 + fastEXP(0.6987 * (V - EK + 11.724))) /
      (1 + fastEXP(0.6168 * (V - EK + 4.872)));
  IK1 = GK1_ * sqrt(Ko / 5.4) * (V - EK) / (1 + xK1);

  // Fast component of the delayed rectifier K current
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
  Jrelinf = alpha_rel * beta_tau * ICaL /
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

  // Set gating variables using lookup table
  ilow = fabs((V - V_min) / 0.1);
  linext = -(-V + lkup[ilow + 1][0]) / 0.1;
  hinf = (lkup[ilow + 1][1] - lkup[ilow][1]) * linext + lkup[ilow + 1][1];
  tauh = (lkup[ilow + 1][2] - lkup[ilow][2]) * linext + lkup[ilow + 1][2];
  tauj = (lkup[ilow + 1][3] - lkup[ilow][3]) * linext + lkup[ilow + 1][3];
  jinf = (lkup[ilow + 1][4] - lkup[ilow][4]) * linext + lkup[ilow + 1][4];
  minf = (lkup[ilow + 1][5] - lkup[ilow][5]) * linext + lkup[ilow + 1][5];
  taum = (lkup[ilow + 1][6] - lkup[ilow][6]) * linext + lkup[ilow + 1][6];
  dinf = (lkup[ilow + 1][7] - lkup[ilow][7]) * linext + lkup[ilow + 1][7];
  taud = (lkup[ilow + 1][8] - lkup[ilow][8]) * linext + lkup[ilow + 1][8];
  finf = (lkup[ilow + 1][9] - lkup[ilow][9]) * linext + lkup[ilow + 1][9];
  tauf = (lkup[ilow + 1][10] - lkup[ilow][10]) * linext + lkup[ilow + 1][10];
  binf = (lkup[ilow + 1][11] - lkup[ilow][11]) * linext + lkup[ilow + 1][11];
  taub = (lkup[ilow + 1][12] - lkup[ilow][12]) * linext + lkup[ilow + 1][12];
  ginf = (lkup[ilow + 1][13] - lkup[ilow][13]) * linext + lkup[ilow + 1][13];
  taug = (lkup[ilow + 1][14] - lkup[ilow][14]) * linext + lkup[ilow + 1][14];
  xKrinf = (lkup[ilow + 1][15] - lkup[ilow][15]) * linext + lkup[ilow + 1][15];
  tauxKr = (lkup[ilow + 1][16] - lkup[ilow][16]) * linext + lkup[ilow + 1][16];
  tauxs1 = (lkup[ilow + 1][17] - lkup[ilow][17]) * linext + lkup[ilow + 1][17];
  tauxs2 = (lkup[ilow + 1][18] - lkup[ilow][18]) * linext + lkup[ilow + 1][18];
  xsinf = (lkup[ilow + 1][19] - lkup[ilow][19]) * linext + lkup[ilow + 1][19];

  // Update gating variables - Euler Method
  h = (hinf - (hinf - h) *
       fastEXP(-DT / tauh)); // Rush-Larsen approximation used for m gate
  j = (jinf - (jinf - j) *
       fastEXP(-DT / tauj)); // Rush-Larsen approximation used for m gate
  m = (minf - (minf - m) *
       fastEXP(-DT / taum)); // Rush-Larsen approximation used for m gate
  d += DT * ((dinf - d)/taud);
  f += DT * ((finf - f)/tauf);
  b += DT * ((binf - b) / taub);
  g += DT * ((ginf - g) / taug);
  xKr += DT * ((xKrinf - xKr)/tauxKr);
  xs1 += DT * ((xsinf - xs1) / tauxs1);
  xs2 += DT * ((xsinf - xs2) / tauxs2);
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
