import numpy as np
def collect_data(watch_data, watch_index, CTRL, ACM, reg_id, reg_iq, reg_speed, fe_htz):
	watch_data[0][watch_index] = fe_htz.psi_2[1]
	watch_data[1][watch_index] = fe_htz.psi_A[1]
	watch_data[2][watch_index] = fe_htz.u_offset[0]
	watch_data[3][watch_index] = fe_htz.u_offset[1]
	watch_data[4][watch_index] = CTRL.cmd_rpm
	watch_data[5][watch_index] = ACM.omega_r_mech/(2*np.pi)*60
	watch_data[6][watch_index] = CTRL.omega_r_elec/(2*np.pi*ACM.npp)*60
	watch_data[7][watch_index] = CTRL.cmd_idq[0]
	watch_data[8][watch_index] = CTRL.idq[0]
	watch_data[9][watch_index] = CTRL.cmd_idq[1]
	watch_data[10][watch_index] = CTRL.idq[1]
	watch_data[11][watch_index] = CTRL.cmd_udq[0]
	watch_data[12][watch_index] = CTRL.cmd_udq[1]
	watch_data[13][watch_index] = CTRL.cmd_uab[0]
	watch_data[14][watch_index] = CTRL.cmd_uab[1]
	watch_data[15][watch_index] = CTRL.omega_syn
	watch_data[16][watch_index] = CTRL.Tem
	watch_data[17][watch_index] = ACM.Tem
	watch_data[18][watch_index] = CTRL.theta_tilde
	watch_data[19][watch_index] = fe_htz.theta_d
	watch_data[20][watch_index] = CTRL.theta_d
	watch_data[21][watch_index] = CTRL.theta_tilde
	watch_data[22][watch_index] = CTRL.window_counter
	watch_data[23][watch_index] = CTRL.iQ_avg_curr
	watch_data[24][watch_index] = CTRL.iQ_avg_prev
	watch_data[25][watch_index] = CTRL.iQ_sum_curr
	watch_data[26][watch_index] = CTRL.iQ_sum_prev
	watch_data[27][watch_index] = fe_htz.psi_2[0]
	watch_data[28][watch_index] = fe_htz.psi_A[0]
	watch_index += 1
	return watch_index