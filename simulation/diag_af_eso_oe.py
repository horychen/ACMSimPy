#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AF OE + ESO OE — Using Framework ESO (v6)
Run simulation with framework ESO enabled (index_separate_speed_estimation=1),
read ESO outputs from watch_data, compute AF OE in post-processing.
ESO bandwidth controlled by omega_ob parameter.
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, os, sys, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from eval_staged_tuning import get_motor_preset

def F(msg): sys.stdout.write(msg+'\n'); sys.stdout.flush()
def wrap_pi(a): return (a + np.pi) % (2*np.pi) - np.pi


def run_with_eso(d, CLBW_Hz, zeta, af_Kp, af_Ki, omega_ob, cmd_rpm, load_step):
    """Run sim with framework ESO at given observer BW, return AF+ESO OE."""
    dd = copy.deepcopy(d)
    R=dd['init_R']; L=dd['init_Lq']; n_pp=dd['init_npp']; KE=dd['init_KE']
    CL_TS=dd['CL_TS']; Js=dd['init_Js']; VL_TS=CL_TS*dd['VL_EXE_PER_CL_EXE']

    cKp,cKi = get_coeffs_dc_motor_current_regulator(R,L,CLBW_Hz)
    sKp,sKi = get_coeffs_dc_motor_SPEED_regulator(Js,n_pp,KE,zeta,cKp/L)
    dd['CL_SERIES_KP']=cKp; dd['CL_SERIES_KI']=cKi; dd['VL_SERIES_KP']=sKp; dd['VL_SERIES_KI']=sKi

    CTRL = The_Motor_Controller(CL_TS=CL_TS,VL_TS=VL_TS,init_npp=n_pp,
        init_IN=dd['init_IN'],init_R=R,init_Ld=dd['init_Ld'],
        init_Lq=L,init_KE=KE,init_Rreq=0.0,init_Js=Js,DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'])
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = True
    CTRL.bool_apply_speed_closed_loop_control = True
    
    # Enable framework ESO but keep encoder angle for FOC
    CTRL.index_separate_speed_estimation = 1
    CTRL.bool_use_sensorless_theta = 0  # encoder angle for Park
    CTRL.use_disturbance_feedforward_rejection = 0
    
    # Set 3rd-order ESO gains (matching framework formula)
    CTRL.ell1 = 3 * omega_ob
    CTRL.ell2 = 3 * omega_ob**2
    CTRL.ell3 = omega_ob**3 * Js / n_pp
    CTRL.ell4 = 0.0

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
    ki_f = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False',10)
    vl_lim = dd.get('VL_LIMIT_OVERLOAD_FACTOR',10.0)
    Vm=dd['DC_BUS_VOLTAGE']/1.732; Im=vl_lim*1.414*dd['init_IN']
    reg_id = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
    reg_iq = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
    reg_spd = The_PID_Regulator(sKp,sKp*sKi,0,0,Im,Im,VL_TS)

    N_slices = dd['NUMBER_OF_SLICES']
    all_t=[]; all_th=[]; all_w=[]; all_ia=[]; all_ib=[]; all_ua=[]; all_ub=[]
    all_eso_th=[]; all_eso_w=[]; all_eso_d=[]
    for sl in range(N_slices):
        t0=sl*dd['TIME_SLICE']; tm=t0+dd['TIME_SLICE']*0.5
        if tm<0.3: CTRL.cmd_rpm=cmd_rpm*min(tm/0.3,1); ACM.TLoad=0
        elif tm<0.6: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<0.8: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=load_step
        elif tm<1.0: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<1.3: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0
        else: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0
        mt,wd = ACMSimPyIncremental(t0=t0,TIME=dd['TIME_SLICE'],ACM=ACM,CTRL=CTRL,
                                     reg_id=reg_id,reg_iq=reg_iq,reg_speed=reg_spd)
        for k in range(len(mt)):
            all_t.append(mt[k])
            all_th.append(wrap_pi(wd[0][k]))  # encoder theta [-pi,pi]
            all_w.append(wd[1][k])             # true speed [rpm]
            all_ia.append(wd[6][k]); all_ib.append(wd[7][k])
            all_ua.append(wd[28][k]); all_ub.append(wd[29][k])
            all_eso_th.append(wd[15][k])       # xS[0] ESO theta
            all_eso_w.append(wd[16][k])        # xS[1] ESO speed [rpm]
            all_eso_d.append(wd[17][k])        # xS[2] ESO dist

    t=np.array(all_t); th=np.array(all_th); w=np.array(all_w)
    ia=np.array(all_ia); ib=np.array(all_ib); ua=np.array(all_ua); ub=np.array(all_ub)
    eso_th=np.array(all_eso_th); eso_w=np.array(all_eso_w); eso_d=np.array(all_eso_d)

    # AF post-processing
    af_oe=np.zeros(len(t)); af_ae=np.zeros(len(t))
    psi_s=np.array([KE,0.0]); corr_int=np.zeros(2)
    for i in range(len(t)):
        a=psi_s[0]-L*ia[i]; b=psi_s[1]-L*ib[i]; amp=np.sqrt(a*a+b*b)
        err=KE-amp; corr_int[0]+=af_Ki*err*CL_TS; corr_int[1]+=af_Ki*err*CL_TS
        if amp>1e-10: un_a,un_b=a/amp,b/amp
        else: un_a,un_b=1.0,0.0
        ca=(af_Kp*err+corr_int[0])*un_a; cb=(af_Kp*err+corr_int[1])*un_b
        psi_s[0]+=CL_TS*(ua[i]-R*ia[i]+ca); psi_s[1]+=CL_TS*(ub[i]-R*ib[i]+cb)
        a2=psi_s[0]-L*ia[i]; b2=psi_s[1]-L*ib[i]
        af_oe[i]=err; af_ae[i]=wrap_pi(th[i]-np.arctan2(b2,a2))

    # ESO OE
    eso_oe = wrap_pi(th - eso_th)

    return {
        't':t, 'th':th, 'w':w, 'KE':KE, 'L':L, 'n_pp':n_pp,
        'af_oe':af_oe, 'af_ae':af_ae,
        'eso_oe':eso_oe, 'eso_w':eso_w, 'eso_d':eso_d,
    }


# ======================== MAIN ========================
if __name__ == '__main__':
    out_dir = os.path.join(os.path.dirname(__file__), 'docs', 'sensorless_report_assets')
    os.makedirs(out_dir, exist_ok=True)

    # omega_ob sweep (rad/s) — framework default is 100
    obs = [50, 100, 200, 500]

    for motor_name in ['servo', 'small_L']:
        p = get_motor_preset(motor_name)
        KE=p['d']['init_KE']; L=p['d']['init_Lq']; n_pp=p['d']['init_npp']; Js=p['d']['init_Js']
        F(f'\n{"="*75}')
        F(f'  {motor_name}: KE={KE*1e3:.1f}mWb L={L*1e3:.2f}mH npp={n_pp} Js={Js}')
        F(f'{"="*75}')
        F(f'  {"wo":>5s} {"bw_Hz":>6s} | {"AF_OE%":>7s} {"AF_ang":>7s} | {"ESO_OE":>8s} {"spd_err":>8s} | ok')
        F(f'  {"-"*62}')

        for wo in obs:
            bw_hz = wo/(2*np.pi)
            t0w=_time.time()
            res = run_with_eso(p['d'], p['clbw'], p['zeta'], p['af_kp'], p['af_ki'],
                               omega_ob=wo, cmd_rpm=p['cmd_rpm'], load_step=p['load_step'])
            dt_wall=_time.time()-t0w

            ss = res['t'] > 0.4
            af_pct = np.sqrt(np.mean(res['af_oe'][ss]**2))/KE*100
            af_ang = np.degrees(np.sqrt(np.mean(res['af_ae'][ss]**2)))
            eso_oe_deg = np.degrees(np.sqrt(np.mean(res['eso_oe'][ss]**2)))
            spd_err = np.sqrt(np.mean((res['eso_w'][ss]-res['w'][ss])**2))
            ok='Y' if eso_oe_deg<1 and spd_err<50 else ('~' if eso_oe_deg<5 else 'N')
            F(f'  {wo:>4d} {bw_hz:>6.1f} | {af_pct:>6.2f}% {af_ang:>6.2f}d | {eso_oe_deg:>7.4f}d {spd_err:>7.2f}r | {ok} ({dt_wall:.1f}s)')

            # Plot for selected BWs
            if wo in [100, 200]:
                fig,axes=plt.subplots(4,1,figsize=(16,14),sharex=True)
                fig.suptitle(f'{motor_name} | AF OE + ESO OE | wo={wo} ({bw_hz:.0f}Hz)', fontsize=13, fontweight='bold')
                axes[0].plot(res['t'],res['w'],'b',lw=0.8,label='true')
                axes[0].plot(res['t'],res['eso_w'],'r--',lw=0.6,alpha=0.7,label=f'ESO (err={spd_err:.2f}rpm)')
                axes[0].set_ylabel('Speed [rpm]'); axes[0].legend(fontsize=9)
                axes[0].set_title('Speed: True vs ESO'); axes[0].grid(True,alpha=0.3)
                axes[1].plot(res['t'],res['af_oe']*1e3,'r',lw=0.6,label=f'AF OE={af_pct:.2f}%KE')
                axes[1].axhline(0,color='k',lw=0.5); axes[1].set_ylabel('AF OE [mWb]')
                axes[1].legend(fontsize=9); axes[1].set_title('AF OE = KE - |psi_AF|'); axes[1].grid(True,alpha=0.3)
                axes[2].plot(res['t'],np.degrees(res['eso_oe']),'b',lw=0.6,label=f'ESO OE={eso_oe_deg:.4f}deg')
                axes[2].axhline(0,color='k',lw=0.5); axes[2].set_ylabel('ESO OE [deg]')
                axes[2].legend(fontsize=9); axes[2].set_title('ESO OE = theta_enc - theta_ESO'); axes[2].grid(True,alpha=0.3)
                axes[3].plot(res['t'],res['eso_w']-res['w'],'r',lw=0.6,label=f'Spd err={spd_err:.2f}rpm')
                ax3r=axes[3].twinx(); ax3r.plot(res['t'],res['eso_d'],'g',lw=0.6,alpha=0.7,label='Dist est')
                axes[3].set_ylabel('Spd err [rpm]'); ax3r.set_ylabel('Disturbance',color='g')
                axes[3].legend(loc='upper left',fontsize=9); ax3r.legend(loc='upper right',fontsize=9)
                axes[3].set_xlabel('Time [s]'); axes[3].grid(True,alpha=0.3)
                plt.tight_layout()
                fn=os.path.join(out_dir,f'oe_af_eso_{motor_name}_wo{wo}.png')
                plt.savefig(fn,dpi=150); plt.close()
                F(f'       -> {os.path.basename(fn)}')

    F('\nAll done!')
