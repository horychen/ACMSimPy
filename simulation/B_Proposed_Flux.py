import numpy as np
from tutorials_ep9_flux_estimator import  RK4_ObserverSolver_CJH_Style, DYNAMICS_FluxEstimator



def rhf_CmdErrFdkCorInFrameRho_Dynamics(x, CTRL, fe_htz):
    
    fx = np.zeros(NS_GLOBAL)
    
    fe_htz.rotor_flux_error[0] = ( CTRL.cmd_psi[0] - (x[0]-CTRL.Lq * CTRL.iab[0]) )
    fe_htz.rotor_flux_error[1] = ( CTRL.cmd_psi[1] - (x[1]-CTRL.Lq * CTRL.iab[1]) )

    fe_htz.emf_stator[0] = CTRL.uab[0] - CTRL.R * CTRL.iab[0] + fe_htz.OFFSET_VOLTAGE_ALPHA + fe_htz.VM_PROPOSED_PI_CORRECTION_GAIN_P * fe_htz.rotor_flux_error[0] + x[2]
    fe_htz.emf_stator[1] = CTRL.uab[0] - CTRL.R * CTRL.iab[1] + fe_htz.OFFSET_VOLTAGE_BETA  + fe_htz.VM_PROPOSED_PI_CORRECTION_GAIN_P * fe_htz.rotor_flux_error[1] + x[3]
    fx[0] = fe_htz.emf_stator[0]
    fx[1] = fe_htz.emf_stator[1]
    fx[2] = fe_htz.VM_PROPOSED_PI_CORRECTION_GAIN_I * fe_htz.rotor_flux_error[0]
    fx[3] = fe_htz.VM_PROPOSED_PI_CORRECTION_GAIN_I * fe_htz.rotor_flux_error[1]
    return fx


def VM_ProposedCmdErrFdkCorInFrameRho(fe_htz, CTRL):
    #/* Proposed VM based Flux Command Error fe_htzedback Correction in Controller Frame (xRho), implemented in AB frame + ODE4 */
    #// stator flux and integral states update
    
    RK4_ObserverSolver_CJH_Style(rhf_CmdErrFdkCorInFrameRho_Dynamics, fe_htz.xFlux, CTRL.CL_TS, CTRL)
    #// Unpack x
    fe_htz.psi_1[0]                           = fe_htz.xFlux[0]
    fe_htz.psi_1[1]                           = fe_htz.xFlux[1]
    fe_htz.correction_integral_term[0]        = fe_htz.xFlux[2]
    fe_htz.correction_integral_term[1]        = fe_htz.xFlux[3]
    fe_htz.u_offset[0] = fe_htz.correction_integral_term[0]
    fe_htz.u_offset[1] = fe_htz.correction_integral_term[1]
   
    #// rotor flux updates

    fe_htz.psi_2[0] = fe_htz.psi_1[0] - CTRL.Lq * CTRL.iab_curr[0]
    fe_htz.psi_2[1] = fe_htz.psi_1[1] - CTRL.Lq * CTRL.iab_curr[1]

    fe_htz.theta_d = np.arctan2(fe_htz.psi_2[1], fe_htz.psi_2[0])
    fe_htz.cosT = np.cos(fe_htz.theta_d)
    fe_htz.sinT = np.sin(fe_htz.theta_d)

