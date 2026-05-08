#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test the framework's built-in ESO by enabling it and reading xS outputs."""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, os, sys, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from eval_staged_tuning import get_motor_preset

def F(msg): sys.stdout.write(msg+'\n'); sys.stdout.flush()
def wrap_pi(a): return (a + np.pi) % (2*np.pi) - np.pi

F('=== Testing Framework ESO ===')
motor_name = 'servo'
p = get_motor_preset(motor_name)
d = copy.deepcopy(p['d'])
R=d['init_R']; L=d['init_Lq']; n_pp=d['init_npp']; KE=d['init_KE']
CL_TS=d['CL_TS']; Js=d['init_Js']; VL_TS=CL_TS*d['VL_EXE_PER_CL_EXE']

cKp,cKi = get_coeffs_dc_motor_current_regulator(R,L,p['clbw'])
sKp,sKi = get_coeffs_dc_motor_SPEED_regulator(Js,n_pp,KE,p['zeta'],cKp/L)
d['CL_SERIES_KP']=cKp; d['CL_SERIES_KI']=cKi; d['VL_SERIES_KP']=sKp; d['VL_SERIES_KI']=sKi

# Case 1: ESO disabled (Stage 0 encoder)
F('\n--- Case 1: ESO disabled (encoder only) ---')
CTRL = The_Motor_Controller(CL_TS=CL_TS,VL_TS=VL_TS,init_npp=n_pp,
    init_IN=d['init_IN'],init_R=R,init_Ld=d['init_Ld'],
    init_Lq=L,init_KE=KE,init_Rreq=0.0,init_Js=Js,DC_BUS_VOLTAGE=d['DC_BUS_VOLTAGE'])
CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
CTRL.bool_apply_sweeping_frequency_excitation = False
CTRL.bool_overwrite_speed_commands = True
CTRL.bool_zero_id_control = True
CTRL.bool_apply_speed_closed_loop_control = True
CTRL.index_separate_speed_estimation = 0  # ENCODER
ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
ki_f = d.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False',10)
vl_lim = d.get('VL_LIMIT_OVERLOAD_FACTOR',10.0)
Vm=d['DC_BUS_VOLTAGE']/1.732; Im=vl_lim*1.414*d['init_IN']
reg_id = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
reg_iq = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
reg_spd = The_PID_Regulator(sKp,sKp*sKi,0,0,Im,Im,VL_TS)

CTRL.cmd_rpm = p['cmd_rpm']; ACM.TLoad = 0
mt, wd = ACMSimPyIncremental(t0=0,TIME=d['TIME_SLICE']*d['NUMBER_OF_SLICES'],ACM=ACM,CTRL=CTRL,
                              reg_id=reg_id,reg_iq=reg_iq,reg_speed=reg_spd)
F(f'  wd[0] (theta_d): [{wd[0][:5]}...]')
F(f'  wd[1] (speed rpm): [{wd[1][:5]}...{wd[1][-5:]}]')
F(f'  wd[15] (xS[0] theta): [{wd[15][:5]}...]')
F(f'  wd[16] (xS[1] speed): [{wd[16][:5]}...]')
F(f'  wd[27] (Tem): [{wd[27][:5]}...]')

# Check ESO gains
F(f'  ESO gains: ell1={CTRL.ell1:.1f}, ell2={CTRL.ell2:.1f}, ell3={CTRL.ell3:.1f}, ell4={CTRL.ell4:.1f}')
F(f'  xS state: {CTRL.xS}')

# Case 2: ESO enabled (index=1), encoder angle used for FOC (bool_use_sensorless_theta=0)
F('\n--- Case 2: ESO enabled, encoder angle for FOC ---')
d2 = copy.deepcopy(p['d']); d2['CL_SERIES_KP']=cKp; d2['CL_SERIES_KI']=cKi; d2['VL_SERIES_KP']=sKp; d2['VL_SERIES_KI']=sKi
CTRL2 = The_Motor_Controller(CL_TS=CL_TS,VL_TS=VL_TS,init_npp=n_pp,
    init_IN=d2['init_IN'],init_R=R,init_Ld=d2['init_Ld'],
    init_Lq=L,init_KE=KE,init_Rreq=0.0,init_Js=Js,DC_BUS_VOLTAGE=d2['DC_BUS_VOLTAGE'])
CTRL2.bool_apply_decoupling_voltages_to_current_regulation = False
CTRL2.bool_apply_sweeping_frequency_excitation = False
CTRL2.bool_overwrite_speed_commands = True
CTRL2.bool_zero_id_control = True
CTRL2.bool_apply_speed_closed_loop_control = True
CTRL2.index_separate_speed_estimation = 1  # ESO ENABLED
CTRL2.bool_use_sensorless_theta = 0  # Still use encoder angle for FOC
CTRL2.use_disturbance_feedforward_rejection = 0
ACM2 = The_AC_Machine(CTRL2, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
reg_id2 = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
reg_iq2 = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
reg_spd2 = The_PID_Regulator(sKp,sKp*sKi,0,0,Im,Im,VL_TS)

F(f'  ESO gains: ell1={CTRL2.ell1:.1f}, ell2={CTRL2.ell2:.1f}, ell3={CTRL2.ell3:.1f}, ell4={CTRL2.ell4:.1f}')
F(f'  xS initial: {CTRL2.xS}')

CTRL2.cmd_rpm = p['cmd_rpm']; ACM2.TLoad = 0
mt2,wd2 = ACMSimPyIncremental(t0=0,TIME=d2['TIME_SLICE']*d2['NUMBER_OF_SLICES'],ACM=ACM2,CTRL=CTRL2,
                               reg_id=reg_id2,reg_iq=reg_iq2,reg_speed=reg_spd2)

F(f'  After sim:')
F(f'  wd2[0] (theta_d): [{wd2[0][:5]}...{wd2[0][-5:]}]')
F(f'  wd2[1] (speed rpm): [{wd2[1][:5]}...{wd2[1][-5:]}]')
F(f'  wd2[15] (xS[0] theta): [{wd2[15][:5]}...{wd2[15][-5:]}]')
F(f'  wd2[16] (xS[1] speed rpm): [{wd2[16][:5]}...{wd2[16][-5:]}]')
F(f'  wd2[17] (xS[2] dist): [{wd2[17][:5]}...{wd2[17][-5:]}]')
F(f'  xS final: {CTRL2.xS}')
F(f'  omega_r_elec (from xS[1]): {CTRL2.omega_r_elec:.2f} rad/s')
F(f'  speed_observer_output_error: {CTRL2.speed_observer_output_error:.6f} rad')

# ESO OE = theta_d(encoder) - xS[0]
theta_enc = wrap_pi(wd2[0])
theta_eso = wd2[15]  # xS[0] is already in [-pi, pi] from framework wrap
eso_oe = wrap_pi(theta_enc - theta_eso)
ss = mt2 > 0.4
oe_rms = np.degrees(np.sqrt(np.mean(eso_oe[ss]**2)))
oe_max = np.degrees(np.max(np.abs(eso_oe[ss])))
speed_err = wd2[16] - wd2[1]  # both in rpm
spd_rms = np.sqrt(np.mean(speed_err[ss]**2))
F(f'\n  Framework ESO results:')
F(f'  ESO OE RMS: {oe_rms:.4f} deg, Max: {oe_max:.4f} deg')
F(f'  Speed err RMS: {spd_rms:.2f} rpm')

# Plot
out_dir = os.path.join(os.path.dirname(__file__), 'docs', 'sensorless_report_assets')
os.makedirs(out_dir, exist_ok=True)
fig,axes=plt.subplots(3,1,figsize=(16,12),sharex=True)
fig.suptitle(f'{motor_name} | Framework ESO | ell1={CTRL2.ell1:.0f}', fontsize=13, fontweight='bold')
axes[0].plot(mt2,wd2[1],'b',lw=0.8,label='true speed [rpm]')
axes[0].plot(mt2,wd2[16],'r--',lw=0.6,label=f'ESO speed (err={spd_rms:.2f}rpm)')
axes[0].legend(fontsize=9); axes[0].set_ylabel('Speed [rpm]'); axes[0].grid(True,alpha=0.3)
axes[1].plot(mt2,np.degrees(eso_oe),'b',lw=0.6,label=f'ESO OE RMS={oe_rms:.4f}deg')
axes[1].axhline(0,color='k',lw=0.5); axes[1].set_ylabel('ESO OE [deg]')
axes[1].legend(fontsize=9); axes[1].grid(True,alpha=0.3)
axes[2].plot(mt2,speed_err,'r',lw=0.6,label=f'Speed err RMS={spd_rms:.2f}rpm')
axes[2].set_ylabel('Speed err [rpm]'); axes[2].set_xlabel('Time [s]')
axes[2].legend(fontsize=9); axes[2].grid(True,alpha=0.3)
plt.tight_layout()
fn=os.path.join(out_dir,'oe_framework_eso_servo.png')
plt.savefig(fn,dpi=150); plt.close()
F(f'  Saved {os.path.basename(fn)}')

F('\nDone!')
