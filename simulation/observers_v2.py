#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pure EMF-PLL observer: no flux integration, directly estimate back-EMF
and use a PLL to lock onto the EMF angle.

Key difference from PLLFluxObserver:
- No flux integration at all (avoids R mismatch drift)
- Current model: di/dt = -R/L*i + u/L - e/L
- EMF = measured_derivative - model_prediction (with known R, L)
- PLL tracks the EMF angle

Also: Enhanced SMO with adaptive gain.
"""
import numpy as np


class EMFDirectPLL:
    """
    Direct EMF estimation from current derivative, plus PLL angle tracking.
    
    Estimated EMF: e = u - R*i - L*di/dt
    Since di/dt is noisy, we use a current observer approach:
      di_est/dt = -R/L * i_est + u/L + e_est/L
      e_est is driven by current error through a proportional gain.
    
    Then PLL tracks angle of e_est.
    """
    def __init__(self, R, L, KE, dt, emf_gain=None, pll_bw=100.0, emf_lpf_fc=500.0):
        self.R = R
        self.L = L
        self.KE = KE
        self.dt = dt
        
        # EMF observer gain (proportional to how fast we want EMF tracking)
        # Larger = faster but noisier
        self.k_emf = emf_gain if emf_gain is not None else R / L * 2  # 2x plant bandwidth
        
        # Current observer state
        self.i_est = np.zeros(2)
        
        # EMF estimate (low-pass filtered)
        self.emf = np.zeros(2)
        tau = 1.0 / (2 * np.pi * emf_lpf_fc)
        self.lpf_a = dt / (tau + dt)
        
        # PLL
        self.theta_pll = 0.0
        self.omega_pll = 0.0
        zeta = 0.707
        wn = pll_bw * 2 * np.pi
        self.Kp_pll = 2 * zeta * wn
        self.Ki_pll = wn ** 2
        
    def step(self, i_alpha, i_beta, u_alpha, u_beta):
        dt, R, L = self.dt, self.R, self.L
        
        # Current estimation error
        err_a = i_alpha - self.i_est[0]
        err_b = i_beta - self.i_est[1]
        
        # EMF estimate from current error (proportional observer)
        emf_raw_a = self.k_emf * L * err_a  # e ≈ k * L * (i - i_est)
        emf_raw_b = self.k_emf * L * err_b
        
        # Update current observer: di/dt = -R/L*i + u/L - e/L
        # We use e_raw as our EMF injection
        self.i_est[0] += dt * (-R/L * self.i_est[0] + u_alpha/L + emf_raw_a/L)
        self.i_est[1] += dt * (-R/L * self.i_est[1] + u_beta /L + emf_raw_b/L)
        
        # Low-pass filter EMF
        self.emf[0] += self.lpf_a * (emf_raw_a - self.emf[0])
        self.emf[1] += self.lpf_a * (emf_raw_b - self.emf[1])
        
        # EMF = KE*omega*[-sin(theta), cos(theta)]
        # theta = atan2(-e_alpha, e_beta)
        theta_raw = np.arctan2(-self.emf[0], self.emf[1])
        
        # PLL
        phase_err = np.sin(theta_raw - self.theta_pll)
        self.omega_pll += self.Ki_pll * phase_err * dt
        omega_out = self.omega_pll + self.Kp_pll * phase_err
        self.theta_pll += omega_out * dt
        self.theta_pll = (self.theta_pll + np.pi) % (2 * np.pi) - np.pi
        
        return self.theta_pll, omega_out, self.emf.copy()


class AdaptiveSMO:
    """
    Enhanced SMO with:
    - Sigmoid function instead of signum (smoother)
    - Adaptive gain based on estimated speed
    - Compensation for phase delay from LPF
    """
    def __init__(self, R, L, KE, npp, dt, base_gain=None, lpf_fc=200.0, pll_bw=80.0):
        self.R = R
        self.L = L
        self.KE = KE
        self.npp = npp
        self.dt = dt
        
        self.base_gain = base_gain if base_gain is not None else 2.0  # V
        self.i_est = np.zeros(2)
        self.emf_filt = np.zeros(2)
        
        # Adaptive LPF
        self.lpf_fc = lpf_fc
        tau = 1.0 / (2 * np.pi * lpf_fc)
        self.lpf_a = dt / (tau + dt)
        
        # PLL
        self.theta_pll = 0.0
        self.omega_pll = 0.0
        zeta = 0.707
        wn = pll_bw * 2 * np.pi
        self.Kp_pll = 2 * zeta * wn
        self.Ki_pll = wn ** 2
    
    def _sigmoid(self, x, delta=0.5):
        return x / (abs(x) + delta)
    
    def step(self, i_alpha, i_beta, u_alpha, u_beta):
        dt, R, L = self.dt, self.R, self.L
        
        err_a = i_alpha - self.i_est[0]
        err_b = i_beta - self.i_est[1]
        
        # Adaptive gain: scale with estimated speed (higher speed = larger EMF)
        est_speed = abs(self.omega_pll)
        est_emf = self.KE * est_speed
        k = max(self.base_gain, est_emf * 1.5)
        
        z_a = k * self._sigmoid(err_a)
        z_b = k * self._sigmoid(err_b)
        
        # Current observer
        self.i_est[0] += dt * (-R/L * self.i_est[0] + u_alpha/L + z_a/L)
        self.i_est[1] += dt * (-R/L * self.i_est[1] + u_beta /L + z_b/L)
        
        # LPF
        self.emf_filt[0] += self.lpf_a * (z_a - self.emf_filt[0])
        self.emf_filt[1] += self.lpf_a * (z_b - self.emf_filt[1])
        
        # Phase compensation: LPF introduces delay ≈ atan(ω/ωc)
        theta_raw = np.arctan2(-self.emf_filt[0], self.emf_filt[1])
        
        # Compensate LPF phase lag
        if abs(self.omega_pll) > 1:
            phase_comp = np.arctan(self.omega_pll / (2 * np.pi * self.lpf_fc))
            theta_raw += phase_comp
        
        # PLL
        phase_err = np.sin(theta_raw - self.theta_pll)
        self.omega_pll += self.Ki_pll * phase_err * dt
        omega_out = self.omega_pll + self.Kp_pll * phase_err
        self.theta_pll += omega_out * dt
        self.theta_pll = (self.theta_pll + np.pi) % (2 * np.pi) - np.pi
        
        return self.theta_pll, omega_out, self.emf_filt.copy()
