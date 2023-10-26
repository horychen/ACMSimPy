import numpy as np
def collect_data(watch_data, watch_index, CTRL, ACM, reg_id, reg_iq, reg_speed, fe_htz):
	watch_data[0][watch_index] = fe_htz.psi_2[0]
	watch_data[1][watch_index] = fe_htz.psi_2[1]
	watch_data[2][watch_index] = fe_htz.u_offset[0]
	watch_data[3][watch_index] = fe_htz.u_offset[1]
	watch_data[4][watch_index] = CTRL.cmd_rpm
	watch_data[5][watch_index] = ACM.omega_r_mech/(2*np.pi)*60
	watch_data[6][watch_index] = CTRL.cmd_idq[0]
	watch_data[7][watch_index] = CTRL.idq[0]
	watch_data[8][watch_index] = CTRL.cmd_idq[1]
	watch_data[9][watch_index] = CTRL.idq[1]
	watch_data[10][watch_index] = CTRL.cmd_udq[0]
	watch_data[11][watch_index] = CTRL.cmd_udq[1]
	watch_data[12][watch_index] = CTRL.cmd_uab[0]
	watch_data[13][watch_index] = CTRL.cmd_uab[1]
	watch_data[14][watch_index] = CTRL.theta_d
	watch_index += 1
	return watch_index