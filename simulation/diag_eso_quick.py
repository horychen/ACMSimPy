#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Minimal ESO OE test v5
Key fixes:
1. theta_enc from wd[0] is [0,2pi), convert to [-pi,pi] for ESO consistency
2. ESO angle wrap matches theta_enc wrap
3. Both AF and ESO operate on encoder data, NO feedback to controller
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, os, sys, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from eval_staged_tuning import get_motor_preset

def F(msg): sys.stdout.write(msg+'\n'); sys.stdout.flush()

def wrap_pi(a):
    """Wrap angle to [-pi, pi]"""
    return (a + np.pi) % (2*np.pi) - np.pi

F('=== Starting ESO OE test v5 ===')

for motor_name in ['servo', 'small_L']:
    p = get_motor_preset(motor_name)
    d = copy.deepcopy(p['d'])
    R=d['init_R']; L=d['init_Lq']; n_pp=d['init_npp']; KE=d['init_KE']
    CL_TS=d['CL_TS']; Js=d['init_Js']; VL_TS=CL_TS*d['VL_EXE_PER_CL_EXE']
    cmd_rpm=p['cmd_rpm']; load_step=p['load_step']
    F(f'\n{"="*70}')
    F(f'Motor: {motor_name}, KE={KE*1e3:.1f}mWb, L={L*1e3:.2f}mH, CL_TS={CL_TS*1e6:.0f}us, npp={n_pp}, Js={Js}')

    cKp,cKi = get_coeffs_dc_motor_current_regulator(R,L,p['clbw'])
    sKp,sKi = get_coeffs_dc_motor_SPEED_regulator(Js,n_pp,KE,p['zeta'],cKp/L)
    d['CL_SERIES_KP']=cKp; d['CL_SERIES_KI']=cKi; d['VL_SERIES_KP']=sKp; d['VL_SERIES_KI']=sKi

    CTRL = The_Motor_Controller(CL_TS=CL_TS,VL_TS=VL_TS,init_npp=n_pp,
        init_IN=d['init_IN'],init_R=R,init_Ld=d['init_Ld'],
        init_Lq=L,init_KE=KE,init_Rreq=0.0,init_Js=Js,DC_BUS_VOLTAGE=d['DC_BUS_VOLTAGE'])
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = True
    CTRL.bool_apply_speed_closed_loop_control = True
    CTRL.index_separate_speed_estimation = 0  # ENCODER
    CTRL.use_disturbance_feedforward_rejection = 0

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
    ki_f = d.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False',10)
    vl_lim = d.get('VL_LIMIT_OVERLOAD_FACTOR',10.0)
    Vm=d['DC_BUS_VOLTAGE']/1.732; Im=vl_lim*1.414*d['init_IN']
    reg_id = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
    reg_iq = The_PID_Regulator(cKp,cKp*cKi*ki_f,0,0,Vm,Vm,CL_TS)
    reg_spd = The_PID_Regulator(sKp,sKp*sKi,0,0,Im,Im,VL_TS)

    F('  Sim...')
    t0w = _time.time()
    all_t=[]; all_th=[]; all_w=[]; all_ia=[]; all_ib=[]; all_ua=[]; all_ub=[]; all_Tem=[]
    N_slices = d['NUMBER_OF_SLICES']
    for sl in range(N_slices):
        t0=sl*d['TIME_SLICE']; tm=t0+d['TIME_SLICE']*0.5
        if tm<0.3: CTRL.cmd_rpm=cmd_rpm*min(tm/0.3,1); ACM.TLoad=0
        elif tm<0.6: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<0.8: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=load_step
        elif tm<1.0: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<1.3: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0
        else: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0
        mt,wd = ACMSimPyIncremental(t0=t0,TIME=d['TIME_SLICE'],ACM=ACM,CTRL=CTRL,
                                     reg_id=reg_id,reg_iq=reg_iq,reg_speed=reg_spd)
        for k in range(len(mt)):
            all_t.append(mt[k])
            # wd[0] is divmod(ACM.theta_d, 2*pi)[1] -> [0, 2pi), convert to [-pi,pi]
            all_th.append(wrap_pi(wd[0][k]))
            all_w.append(wd[1][k])  # wd[1] = omega_r_mech in rpm
            all_ia.append(wd[6][k]); all_ib.append(wd[7][k])
            all_ua.append(wd[28][k]); all_ub.append(wd[29][k])
            all_Tem.append(wd[27][k])  # wd[27] = CTRL.Tem

    t=np.array(all_t); th=np.array(all_th); w=np.array(all_w)
    ia=np.array(all_ia); ib=np.array(all_ib); ua=np.array(all_ua); ub=np.array(all_ub)
    Tem=np.array(all_Tem)
    F(f'  Sim done: {len(t)} steps, {_time.time()-t0w:.1f}s')
    F(f'  theta range: [{th.min():.3f}, {th.max():.3f}] rad')
    F(f'  speed range: [{w.min():.1f}, {w.max():.1f}] rpm')
    F(f'  Tem range: [{Tem.min():.4f}, {Tem.max():.4f}] Nm')

    # ---- AF ----
    t0w=_time.time()
    af_oe=np.zeros(len(t)); af_ae=np.zeros(len(t))
    psi_s=np.array([KE,0.0]); corr_int=np.zeros(2)
    af_Kp=p['af_kp']; af_Ki=p['af_ki']
    for i in range(len(t)):
        a=psi_s[0]-L*ia[i]; b=psi_s[1]-L*ib[i]; amp=np.sqrt(a*a+b*b)
        err=KE-amp; corr_int[0]+=af_Ki*err*CL_TS; corr_int[1]+=af_Ki*err*CL_TS
        if amp>1e-10: un_a,un_b=a/amp,b/amp
        else: un_a,un_b=1.0,0.0
        ca=(af_Kp*err+corr_int[0])*un_a; cb=(af_Kp*err+corr_int[1])*un_b
        psi_s[0]+=CL_TS*(ua[i]-R*ia[i]+ca); psi_s[1]+=CL_TS*(ub[i]-R*ib[i]+cb)
        a2=psi_s[0]-L*ia[i]; b2=psi_s[1]-L*ib[i]
        af_oe[i]=err; af_ae[i]=wrap_pi(th[i]-np.arctan2(b2,a2))
    F(f'  AF: {_time.time()-t0w:.1f}s')

    ss = t > 0.4
    af_pct = np.sqrt(np.mean(af_oe[ss]**2))/KE*100
    af_ang = np.degrees(np.sqrt(np.mean(af_ae[ss]**2)))
    F(f'  AF OE: {af_pct:.2f}% KE, angle: {af_ang:.2f} deg')

    # ---- ESO sweep ----
    out_dir = os.path.join(os.path.dirname(__file__), 'docs', 'sensorless_report_assets')
    os.makedirs(out_dir, exist_ok=True)

    F(f'\n  {"BW":>5s} | {"ESO_OE":>10s} {"spd_err":>10s} | ok')
    F(f'  {"-"*42}')

    for bw in [5, 10, 20, 50, 100]:
        t0w=_time.time()
        wo=2*np.pi*bw
        ell1=4*wo; ell2=6*wo**2; ell3=4*wo**3; ell4=wo**4
        x=np.zeros(4); x[0]=th[0]  # init from encoder
        eso_oe=np.zeros(len(t)); eso_w=np.zeros(len(t)); eso_d=np.zeros(len(t))
        
        for i in range(len(t)):
            # Output error (angle diff with wrap)
            oe = wrap_pi(th[i]-x[0])
            eso_oe[i]=oe
            eso_w[i]=x[1]/(2*np.pi*n_pp)*60  # elec rad/s -> rpm
            eso_d[i]=x[2]
            
            # RK4
            def f(xx):
                d=wrap_pi(th[i]-xx[0])
                return np.array([
                    ell1*d+xx[1],
                    ell2*d+(Tem[i]+xx[2])*n_pp/Js,
                    ell3*d+xx[3],
                    ell4*d
                ])
            dt=CL_TS
            k1=f(x)*dt; k2=f(x+k1*0.5)*dt; k3=f(x+k2*0.5)*dt; k4=f(x+k3)*dt
            x=x+(k1+2*k2+2*k3+k4)/6.0
            x[0]=wrap_pi(x[0])
        
        oe_rms=np.degrees(np.sqrt(np.mean(eso_oe[ss]**2)))
        spd_e=np.sqrt(np.mean((eso_w[ss]-w[ss])**2))
        ok='Y' if oe_rms<1 and spd_e<50 else ('~' if oe_rms<5 else 'N')
        F(f'  {bw:>4d}Hz | {oe_rms:>10.4f}d {spd_e:>10.2f}r | {ok} ({_time.time()-t0w:.1f}s)')

        if bw in [20,50,100]:
            fig,axes=plt.subplots(4,1,figsize=(16,14),sharex=True)
            fig.suptitle(f'{motor_name} | AF OE + ESO OE | ESO BW={bw}Hz', fontsize=13, fontweight='bold')
            axes[0].plot(t,w,'b',lw=0.8,label='true')
            axes[0].plot(t,eso_w,'r--',lw=0.6,alpha=0.7,label=f'ESO (err={spd_e:.2f}rpm)')
            axes[0].set_ylabel('Speed [rpm]'); axes[0].legend(fontsize=9)
            axes[0].set_title('Speed: True vs ESO'); axes[0].grid(True,alpha=0.3)
            axes[1].plot(t,af_oe*1e3,'r',lw=0.6,label=f'AF OE={af_pct:.2f}%KE')
            axes[1].axhline(0,color='k',lw=0.5); axes[1].set_ylabel('AF OE [mWb]')
            axes[1].legend(fontsize=9); axes[1].set_title('AF OE = KE - |psi_AF|'); axes[1].grid(True,alpha=0.3)
            axes[2].plot(t,np.degrees(eso_oe),'b',lw=0.6,label=f'ESO OE RMS={oe_rms:.4f}deg')
            axes[2].axhline(0,color='k',lw=0.5); axes[2].set_ylabel('ESO OE [deg]')
            axes[2].legend(fontsize=9); axes[2].set_title('ESO OE = theta_enc - theta_ESO'); axes[2].grid(True,alpha=0.3)
            axes[3].plot(t,eso_w-w,'r',lw=0.6,label=f'Spd err RMS={spd_e:.2f}rpm')
            ax3r=axes[3].twinx(); ax3r.plot(t,eso_d,'g',lw=0.6,alpha=0.7,label='Dist est')
            axes[3].set_ylabel('Spd err [rpm]'); ax3r.set_ylabel('Disturbance',color='g')
            axes[3].legend(loc='upper left',fontsize=9); ax3r.legend(loc='upper right',fontsize=9)
            axes[3].set_xlabel('Time [s]'); axes[3].grid(True,alpha=0.3)
            plt.tight_layout()
            fn=os.path.join(out_dir,f'oe_af_eso_{motor_name}_bw{bw}.png')
            plt.savefig(fn,dpi=150); plt.close()
            F(f'         -> {os.path.basename(fn)}')

F('\nAll done!')
