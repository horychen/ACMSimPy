# %%
import multiprocessing
from acmsimpy import *
from datetime import time



# %%
""" DEFINE PARALLEL RUN """
def RunSim_Parallel(end=7, TIME=5.5, MACHINE_TS=1e-4, CL_TS=1e-4):

    # output variables
    output_control_times = np.arange(0, TIME, CL_TS)
    list_speed = [np.zeros_like(output_control_times) for i in range(end-1)]

    # parameter analysis
    omega_ob_values = np.arange(1, end, 1)

    # parallel loop
    for parallel_index in prange(1, int(len(omega_ob_values))):

        # get parameter
        omega_ob = omega_ob_values[parallel_index]

        # Main script receives the parameter
        controL_times, id, iq, ia, ib, speed = ACMSimPy(
            omega_ob=omega_ob, TIME=TIME, MACHINE_TS=MACHINE_TS, CL_TS=CL_TS)

        # output variables
        # if parallel_index == 1:
        #     output_control_times = controL_times
        list_speed[parallel_index] = speed

    return output_control_times, list_speed

RunSim_Parallel(end=2); # Compile


# %%
%%time
""" Series RUN """
f = njit(nogil=True, parallel=False)(RunSim_Parallel)
controL_times, list_speed = f()


# %%
%%time
""" PARALLEL RUN """
f = njit(nogil=True, parallel=True)(RunSim_Parallel)
controL_times, list_speed = f()
f.parallel_diagnostics(level=4)


# %%
""" PLOT """
plt.figure(figsize=(12, 6))
for i in range(len(list_speed)):
    print('\t', min(list_speed[i]))
    plt.plot(controL_times, list_speed[i], label=str(i))
plt.legend(loc='upper left')



