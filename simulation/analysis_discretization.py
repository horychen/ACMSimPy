#!/usr/bin/env python
"""
=======================================================================
  Future Work Item 4: Discretization Analysis
=======================================================================
Analyze the relationship between CL_TS, CLBW, and Nyquist limit.

Key question: Is CLBW=320 Hz at CL_TS=100μs safe from discretization artifacts?
"""
import numpy as np

CL_TS = 1e-4  # 100 μs sampling period
f_s = 1 / CL_TS  # 10 kHz sampling frequency
f_nyquist = f_s / 2  # 5 kHz Nyquist frequency

# Motor parameters (small_L)
R = 0.035      # Ω
L = 0.036e-3   # H
npp = 22
KE = 0.0125    # Wb

tau_LR = L / R  # L/R time constant

print('='*60)
print('  DISCRETIZATION ANALYSIS')
print('='*60)

print(f'\nSampling:')
print(f'  CL_TS = {CL_TS*1e6:.0f} μs')
print(f'  f_s = {f_s:.0f} Hz')
print(f'  f_Nyquist = {f_nyquist:.0f} Hz')

print(f'\nPlant (small_L):')
print(f'  L/R time constant = {tau_LR*1e3:.2f} ms → f_plant = {1/(2*np.pi*tau_LR):.0f} Hz')
print(f'  Plant pole at s = -R/L = {-R/L:.0f} rad/s = {R/L/(2*np.pi):.0f} Hz')

print(f'\n--- CLBW Analysis ---')
clbws = [100, 200, 320, 500, 800, 1000, 2000, 3000]
for clbw in clbws:
    # Current regulator Kp = CLBW * 2π * L
    omega_cl = clbw * 2 * np.pi
    Kp = omega_cl * L
    Ki = Kp * (R / L)  # typical Ki
    
    # Ratio to Nyquist
    ratio = clbw / f_nyquist
    
    # Discrete stability: z-domain pole
    # For first-order plant with PI: closed-loop pole ≈ exp(-omega_cl * CL_TS)
    z_pole = np.exp(-omega_cl * CL_TS)
    
    # Phase margin estimate: PM ≈ 90° - atan(omega_cl / (R/L))
    pm_deg = 90 - np.degrees(np.arctan(omega_cl / (R / L)))
    
    # Rule of thumb: CLBW should be < f_s/10 for safe discretization
    safe = '✅' if clbw < f_s / 10 else ('⚠' if clbw < f_s / 5 else '❌')
    
    print(f'  CLBW={clbw:>5d} Hz | f_cl/f_Nyq={ratio:>5.1%} | z_pole={z_pole:.4f} | PM≈{pm_deg:.0f}° | {safe}')

print(f'\n--- Rule of Thumb ---')
print(f'  f_s/10 = {f_s/10:.0f} Hz (conservative limit)')
print(f'  f_s/5  = {f_s/5:.0f} Hz (aggressive limit)')
print(f'  f_s/2  = {f_s/2:.0f} Hz (Nyquist, absolute limit)')
print(f'\n  → CLBW=320 Hz is at {320/f_s*100:.1f}% of f_s = well within safe zone ✅')

# For big_L
print(f'\n--- big_L comparison ---')
L_big = 0.1035
R_big = 1.97
tau_big = L_big / R_big
f_plant_big = 1 / (2 * np.pi * tau_big)
print(f'  L/R time constant = {tau_big*1e3:.1f} ms → f_plant = {f_plant_big:.1f} Hz')
print(f'  Max safe CLBW (10x plant pole) = {f_plant_big*10:.0f} Hz')
print(f'  Current preset CLBW = 200 Hz → {200/f_plant_big:.1f}x plant pole')
print(f'  This is {200/(f_s/10)*100:.0f}% of conservative limit → safe from discretization')
print(f'  But 200 Hz = {200*2*np.pi*tau_big:.1f}x L/R bandwidth → may be too aggressive for plant!')

# Speed loop analysis
print(f'\n--- Speed Loop Discretization ---')
VL_TS = CL_TS * 5  # VL_EXE_PER_CL_EXE = 5
f_vl = 1 / VL_TS
print(f'  VL_TS = {VL_TS*1e3:.1f} ms → f_vl = {f_vl:.0f} Hz')
print(f'  VL Nyquist = {f_vl/2:.0f} Hz')
print(f'  Typical VLBW = 40-120 Hz → {120/f_vl*100:.0f}% of f_vl = safe')

# Electrical frequency at operating speed
print(f'\n--- Operating Frequencies ---')
for name, ke, pp, rpm, cl in [('small_L', 0.0125, 22, 200, 320),
                                ('servo', 0.1, 4, 500, 500),
                                ('big_L', 0.0745, 24, 150, 80)]:
    f_elec = pp * rpm / 60
    samples_per_cycle = f_s / f_elec
    print(f'  {name:>7s}: f_elec={f_elec:>6.0f} Hz, {samples_per_cycle:>5.0f} samples/cycle, CLBW/f_elec={cl/f_elec:.1f}')
