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

#define exp RTMath.fastEXP

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
  // Ko = 4.5 ; // mM
  // Nao = 140 ; // mM
  // Cao = 1.8 ; // mM
  // Concentration used in dynamic clamp experiments
  Ko = 5.4 ; // mM
  Nao = 137; // mM
  Cao = 2.0 ; // mM
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

  // Gating Variable Lookup Table
  V_min = -200;
  V_step = 0.01;
  int s = V_min * -2 / V_step;
  lkup = new double[s][20];
  Vx = 0; // Voltage placeholder for lookup table

  // Lookup Table Initialization
  for (z = 0; z < s; z++) {
    Vx = V_min + V_step * z;

    // Voltage
    lkup[z][0] = V_min + V_step * z;

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
      am = 3.199985789461998;
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
      lkup[z][8] = 2.289374849326888; // taud
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
      tau_xs1 = 417.9441667499822;
    else
      tau_xs1 = 10000 / (0.719 * (Vx + 30) / (1 - exp(-0.148 * (Vx + 30))) +
                         1.31 * (Vx + 30) / (exp(0.0687 * (Vx + 30)) - 1));
    lkup[z][15] = 1 / (1 + exp(-(Vx + 21.5) / 7.5)); // xKrinf
    if ( Vx > -14.21 && Vx < -14.19 ) // if V = -14.2, divide by 0 error
      lkup[z][16] = 85.830287334611480; // tauxKr
    else
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
  // Buffering
  Cai = calcium_buffer(Cai_t, TRPNtot, KmTRPN, CMDNtot, KmCMDN);
  CaJSR = calcium_buffer(CaJSR_t, CSQNtot, KmCSQN, 0, 0);

  // Reversel Potentials
  ENa = RTF * log(Nao / Nai);
  EK = RTF * log(Ko / Ki);
  EKs = RTF * log((Ko + pKNa * Nao)/(Ki + pKNa * Nai));
  ECa = 0.5 * RTF * log(Cao / Cai);

  // Na currents
  // H-gate
  INa = GNa_ * m * m * m * h * j * (V - ENa);
  INab = GNab * (V - ENa);

  // L-type Ca current
  ICa_ = PCa * 4.0 * F * FRT * V * (gamma_Cai * Cai * exp(2.0 * V *FRT) -
                                  gamma_Cao *Cao) / (exp(2.0 * V *FRT) - 1.0);
  ICaK_ = PCa_K * F * FRT * V * (gamma_Ki *Ki * exp(V * FRT) -
                                 gamma_Ko * Ko) / (exp(V * FRT) - 1.0);
  ICaNa_ = PCa_Na * F * FRT * V * (gamma_Nai * Nai * exp(V * FRT) -
                                   gamma_Nao * Nao) / (exp(V * FRT) - 1.0) ;
  fCa = 1.0 / (Cai / KmCa + 1.0);
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
  xK1 = 0.004 * (1.0 + exp(0.6987 * (V - EK + 11.724))) /
      (1.0 + exp(0.6168 * (V - EK + 4.872)));
  IK1 = GK1_ * sqrt(Ko / 5.4) * (V - EK) / (1.0 + xK1);

  // Fast component of the delayed rectifier K current
  RKr = 1.0 / (exp((V + 9.0) / 22.4) + 1.0);
  IKr = GKr_ * sqrt(Ko / 5.4) * xKr * RKr * (V - EK);

  // Fast component of the delayed rectifier K current
  //IKs = GKs_ * (1 + 0.6/(pow(3.8e-5/Cai,1.4)+1)) * xs1 * xs2 * (V - EKs);
  IKs = GKs_ * (1.0 + 0.6 / (exp(1.4 * log(3.8e-5 / Cai)) + 1.0)) * xs1 *
      xs2 * (V - EKs); // pow() removed

  // Plateau K current
  Kp = 1.0 / (1.0 + exp((7.488 - V) / 5.98));
  IKp = GKp_ * Kp * (V - EK);

  // Pumps and transporters
  // Na-K pump
  sigma_NaK = (exp(Nao / 67.3) - 1) / 7.0;
  fNaK = 1.0/(1.0 + 0.1245 * exp(-0.1 * V * FRT) + 0.0365 * sigma_NaK *
            exp(-V * FRT));
  //INaK = INaK_ * fNaK * Ko / ( (Ko + KmK_NaK) * pow( 1 +
  //((KmNa_NaK/Nai)*(KmNa_NaK/Nai)),2) );
  INaK = INaK_ * fNaK * Ko / ((Ko + KmK_NaK) *
                              (1.0 + ((KmNa_NaK / Nai) *
                                    (KmNa_NaK / Nai)))); // pow() removed

  // Na-Ca exchanger
  INCX = GNCX_ * kNCX * exp((eta - 1.0) * V * FRT) *
      ((Nai * Nai * Nai) * Cao * exp(V * FRT) - (Nao * Nao * Nao) * Cai) /
      ( 1.0 + ksat * exp((eta - 1.0) * V * FRT) *
        ((Nai * Nai * Nai) * Cao * exp(V * FRT) +
         (Nao * Nao * Nao) * Cai));

  // Intracellular Ca fluxes
  // SR Ca release, uptake, and leak
  //Jrelinf = alpha_rel * beta_tau * ICaL / (pow((Krel_inf/CaJSR),hrel) + 1);
  // Added GJrel_ as a scaling factor
  Jrelinf = GJrel_ * alpha_rel * beta_tau * ICaL /
      (exp(hrel * log(Krel_inf / CaJSR)) + 1);
  tau_rel = beta_tau / (Krel_tau / CaJSR + 1);
  dJreldt = - (Jrelinf + Jrel) / tau_rel;

  Jserca = Jserca_ * Vserca * (Cai / (Cai + Kmserca) - CaNSR / CaNSR_max);

  Jtr = (CaNSR - CaJSR) / tau_transfer;

  // Total Current
  NaIon = INa + INab + 3 * INCX + ICaL_Na + 3 * INaK;
  KIon = I_Inject + IKr + IKs + IK1 + IKp + ICaL_K - 2 * INaK;
  CaIon = ICaL + ICab + IpCa - 2 * INCX + ICaT;
  Iion = NaIon + KIon + CaIon;

  // Derivatives for ionic concentration
  dNai = -NaIon * Acap / (Vmyo * F);
  // Injected current assumed to be due to potassium flux
  dKi = -KIon * Acap / (Vmyo * F);
  dCai_t = -Jserca * VNSR / Vmyo + Jrel * VJSR / Vmyo -
               CaIon * Acap / (2 * Vmyo * F);
  dCaJSR_t = Jtr - Jrel;
  dCaNSR = Jserca - Jtr * VJSR / VNSR;

  // Derivative for voltage
  dVdt = -(Iion);

  // Update voltage and ionic concentrations
  V += DT * dVdt;
  Nai += DT * dNai;
  Ki += DT * dKi;
  Cai_t += DT * dCai_t;
  CaJSR_t += DT * dCaJSR_t;
  CaNSR += DT * dCaNSR;
  Jrel += DT * dJreldt;

  // Set gating variables using lookup table
  ilow = fabs((V - V_min) / V_step);
  linext = -(-V + lkup[ilow + 1][0]) / V_step;
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
  h = (hinf - (hinf - h) * exp(-DT / tauh));
  j = (jinf - (jinf - j) * exp(-DT / tauj));
  m = (minf - (minf - m) * exp(-DT / taum));
  d = (dinf - (dinf - d) * exp(-DT / taud));
  f = (finf - (finf - f) * exp(-DT / tauf));
  b = (binf - (binf - b) * exp(-DT / taub));
  g = (ginf - (ginf - g) * exp(-DT / taug));
  xKr = (xKrinf - (xKrinf - xKr) * exp(-DT / tauxKr));
  xs1 = (xsinf - (xsinf - xs1) * exp(-DT / tauxs1));
  xs2 = (xsinf - (xsinf - xs2) * exp(-DT / tauxs2));
}

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
  V = -84.6684;
  Cai_t = 0.0231045;
  CaNSR = 2.64824;
  CaJSR_t = 9.04718;
  Nai = 14.2605;
  Ki = 137.491;
  m = 0.00163319;
  h = 3.66776;
  j = 0.989748;
  d = 9.88131e-324;
  f = 0.998927;
  b = 0.00143909;
  g = 0.975127;
  xKr = 0.000252203;
  xs1 = 0.0227544;
  xs2 = 0.0746654;
  Jrel = 2.34309e-39;

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
