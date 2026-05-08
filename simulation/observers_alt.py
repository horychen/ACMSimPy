#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PLL-based and Sliding Mode Observer (SMO) for sensorless PMSM control.
Drop-in replacements for the Active Flux observer in demo_sensorless_active_flux.py.
"""
import numpy as np


class PLLFluxObserver:
    """
    PLL-based flux observer for sensorless PMSM.
    
    Instead of feeding the AF angle directly into Park transform (creating
    positive feedback), this observer:
    1. Integrates stator flux from voltage model (same as AF)
    2. Extracts angle via atan2 (same as AF)  
    3. Uses a PLL to track the angle smoothly (NEW - this is the key difference)
    
    The PLL acts as a 2nd-order low-pass filter on the angle, providing
    inherent noise rejection and avoiding the direct angle→Park→voltage→angle
    positive feedback loop that destabilizes AF on high-L motors.
    """
    def __init__(self, R_obs, L, KE, dt, pll_bw=100.0, af_Kp=500.0, af_Ki=5000.0):
        self.R = R_obs
        self.L = L
        self.KE = KE
        self.dt = dt
        
        # AF flux integrator state
        self.psi_s = np.array([KE, 0.0])
        self.corr_int = np.zeros(2)
        self.af_Kp = af_Kp
        self.af_Ki = af_Ki
        
        # PLL state: [theta_pll, omega_pll]
        self.theta_pll = 0.0
        self.omega_pll = 0.0
        
        # PLL gains (2nd-order, damping ratio = 0.707)
        zeta_pll = 0.707
        wn = pll_bw * 2 * np.pi  # natural frequency
        self.Kp_pll = 2 * zeta_pll * wn
        self.Ki_pll = wn ** 2
        
    def step(self, i_alpha, i_beta, u_alpha, u_beta):
        """Run one step. Returns (theta_est, omega_est_elec_rad_s)."""
        dt = self.dt
        L, KE = self.L, self.KE
        
        # 1. Active flux from voltage model (same as original AF)
        af_alpha = self.psi_s[0] - L * i_alpha
        af_beta  = self.psi_s[1] - L * i_beta
        af_amp = np.sqrt(af_alpha**2 + af_beta**2)
        
        # PI amplitude correction
        psi_err = KE - af_amp
        self.corr_int += self.af_Ki * psi_err * dt
        
        if af_amp > 1e-10:
            u_a, u_b = af_alpha / af_amp, af_beta / af_amp
        else:
            u_a, u_b = 1.0, 0.0
            
        corr_a = (self.af_Kp * psi_err + self.corr_int[0]) * u_a
        corr_b = (self.af_Kp * psi_err + self.corr_int[1]) * u_b
        
        self.psi_s[0] += dt * (u_alpha - self.R * i_alpha + corr_a)
        self.psi_s[1] += dt * (u_beta  - self.R * i_beta  + corr_b)
        
        # Recompute active flux after integration
        af_alpha = self.psi_s[0] - L * i_alpha
        af_beta  = self.psi_s[1] - L * i_beta
        
        # 2. Raw angle from atan2
        theta_raw = np.arctan2(af_beta, af_alpha)
        
        # 3. PLL tracking: phase detector → loop filter → VCO
        # Phase error = sin(theta_raw - theta_pll) ≈ theta_raw - theta_pll for small errors
        phase_err = np.sin(theta_raw - self.theta_pll) * np.cos(theta_raw - self.theta_pll)
        # Simplified: use sin for large-signal stability
        # phase_err = np.sin(theta_raw - self.theta_pll)
        
        # PI loop filter
        self.omega_pll += self.Ki_pll * phase_err * dt
        omega_out = self.omega_pll + self.Kp_pll * phase_err
        
        # VCO (integrator)
        self.theta_pll += omega_out * dt
        self.theta_pll = (self.theta_pll + np.pi) % (2 * np.pi) - np.pi
        
        return self.theta_pll, omega_out, (af_alpha, af_beta)


class SlidingModeObserver:
    """
    Sliding Mode Observer (SMO) for sensorless PMSM.
    
    Estimates back-EMF directly from current estimation error using
    a signum function, then extracts angle via PLL.
    
    Model: di/dt = -R/L * i + u/L - e/L
    where e = [eα, eβ] = KE*ω*[-sin(θ), cos(θ)] is back-EMF
    
    The SMO estimates current, and uses the sliding mode term
    (signum of current error) to estimate back-EMF.
    """
    def __init__(self, R, L, KE, dt, smo_gain=None, lpf_fc=None, pll_bw=50.0):
        self.R = R
        self.L = L
        self.KE = KE
        self.dt = dt
        
        # SMO gain: must be > max|e| = KE * omega_max
        # Default: 5x rated back-EMF for robustness
        self.k_smo = smo_gain if smo_gain is not None else KE * 500 * 2 * np.pi / 60 * 5
        
        # Estimated currents
        self.i_est = np.zeros(2)  # [i_alpha_est, i_beta_est]
        
        # Low-pass filtered back-EMF estimates
        self.emf_filt = np.zeros(2)  # [e_alpha_filt, e_beta_filt]
        lpf_fc = lpf_fc if lpf_fc is not None else 200.0  # Hz
        tau = 1.0 / (2 * np.pi * lpf_fc)
        self.lpf_alpha = dt / (tau + dt)
        
        # PLL for angle extraction from filtered EMF
        self.theta_pll = 0.0
        self.omega_pll = 0.0
        zeta_pll = 0.707
        wn = pll_bw * 2 * np.pi
        self.Kp_pll = 2 * zeta_pll * wn
        self.Ki_pll = wn ** 2
        
    def _sat(self, x, delta=0.1):
        """Saturation function (smooth approximation of signum)."""
        if abs(x) < delta:
            return x / delta
        return 1.0 if x > 0 else -1.0
    
    def step(self, i_alpha, i_beta, u_alpha, u_beta):
        """Run one step. Returns (theta_est, omega_est_elec_rad_s)."""
        dt = self.dt
        R, L, k = self.R, self.L, self.k_smo
        
        # Current estimation error
        err_a = i_alpha - self.i_est[0]
        err_b = i_beta  - self.i_est[1]
        
        # Sliding mode switching term (back-EMF estimate)
        z_a = k * self._sat(err_a)
        z_b = k * self._sat(err_b)
        
        # Current observer: di_est/dt = -R/L * i_est + u/L + z/L
        self.i_est[0] += dt * (-R/L * self.i_est[0] + u_alpha/L + z_a/L)
        self.i_est[1] += dt * (-R/L * self.i_est[1] + u_beta /L + z_b/L)
        
        # Low-pass filter the switching signal to get smooth back-EMF
        self.emf_filt[0] += self.lpf_alpha * (z_a - self.emf_filt[0])
        self.emf_filt[1] += self.lpf_alpha * (z_b - self.emf_filt[1])
        
        # Raw angle from filtered EMF: e = KE*ω*[-sin(θ), cos(θ)]
        # So θ = atan2(-eα, eβ)
        theta_raw = np.arctan2(-self.emf_filt[0], self.emf_filt[1])
        
        # PLL tracking
        phase_err = np.sin(theta_raw - self.theta_pll)
        self.omega_pll += self.Ki_pll * phase_err * dt
        omega_out = self.omega_pll + self.Kp_pll * phase_err
        self.theta_pll += omega_out * dt
        self.theta_pll = (self.theta_pll + np.pi) % (2 * np.pi) - np.pi
        
        return self.theta_pll, omega_out, self.emf_filt.copy()
