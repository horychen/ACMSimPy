@njit(nogil=True)
def DYNAMICS_MACHINE(t, x, ACM, CLARKE_TRANS_TORQUE_GAIN=1.5):
    fx = np.zeros(ACM.NS)

    # omega_d_mech = x[0]
    # omega_r_mech = x[1]
    KA    = x[2]
    iD    = x[3]
    iQ    = x[4]
    # ACM.theta_d = x[0]*ACM.npp
    # ACM.omega_r = x[1]*ACM.npp
    ACM.omega_slip = ACM.Rreq * iQ / KA
    ACM.omega   = x[1]*ACM.npp + ACM.omega_slip

    # 电磁子系统 (KA, iD, iQ as x[2], x[3], x[4])
    if ACM.Rreq > 0:
        # s KA
        fx[2] = ACM.Rreq*iD - ACM.Rreq / (ACM.Ld - ACM.Lq) * KA # [Apply Park Transorm to (31b)]
        # s iD
        fx[3] = (ACM.udq[0] - ACM.R*iD + ACM.omega*ACM.Lq*iQ - fx[2]) / ACM.Lq # (6a)
    else: 
            # note fx[3] * ACM.Lq = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ - fx[2]
            #  =>  fx[3] * ACM.Lq = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ - (ACM.Ld - ACM.Lq) * fx[3] - 0.0
            #  =>  fx[3] * ACM.Ld = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ
            #  =>  s iD
        # s iD
        fx[3] = (ACM.udq[0] - ACM.R*iD + ACM.omega*ACM.Lq*iQ) / ACM.Ld
        # s KA
        fx[2] = (ACM.Ld - ACM.Lq) * fx[3] + 0.0
    # s iQ
    fx[4] = (ACM.udq[1] - ACM.R*iQ - ACM.omega*ACM.Lq*iD - ACM.omega*ACM.KActive) / ACM.Lq

    # 机械子系统 (theta_d_mech, omega_mech as x[0], x[1])
    ACM.Tem = CLARKE_TRANS_TORQUE_GAIN * ACM.npp * KA * iQ # 电磁转矩计算
    fx[0] = x[1] + ACM.omega_slip / ACM.npp # mech. angular rotor position (accumulated)
    fx[1] = (ACM.Tem - ACM.TLoad) / ACM.Js  # mech. angular rotor speed

    return fx
