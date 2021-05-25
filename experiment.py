# %%
from acmsimpy import *
from datetime import time

# Compile
ACMSimPy(omega_ob=1, TIME=0.1, MACHINE_TS=1e-4, CL_TS=1e-4)


# %%
%%time
@njit(nogil=True, parallel=False)  # THIS IS ACTUALLY FASTER
def RUN(end, TIME=5.5, MACHINE_TS=1e-4, CL_TS=1e-4):
    control_times = np.arange(0, TIME, CL_TS)
    list_speed = [np.zeros_like(control_times) for i in range(end-1)]
    parallel_index = 0

    for omega_ob in prange(1, int(end)):
        print(omega_ob, ' | ')
        controL_times, id, iq, ia, ib, speed = ACMSimPy(
            omega_ob=omega_ob, TIME=TIME, MACHINE_TS=MACHINE_TS, CL_TS=CL_TS)
        list_speed[parallel_index] = speed
        parallel_index += 1
    return controL_times, list_speed

end = 12+1
controL_times, list_speed = RUN(end=end, TIME=5.5, MACHINE_TS=1e-4)
""" PLOT """
plt.figure(figsize=(12,6))
for i in range(len(list_speed)):
    print('\t', min(list_speed[i]))
    plt.plot(controL_times, list_speed[i], label=str(i))
plt.legend(loc='upper left');




























# %%
%%time
@njit(nogil=True, parallel=True) # THIS IS NOT WORKING AS PARALLEL!!!
def RUN(end=13):
    for omega_ob in prange(1, int(end)):
        # ACMSimPy(omega_ob=omega_ob, TIME=5.5, MACHINE_TS=1e-4, CL_TS=1e-4)
        ACMSimPy(omega_ob, 5.5, 1e-4, 1e-4)
        # print('TESA')
RUN()
RUN.parallel_diagnostics(level=4)




# %%
%%time
if False:
    # Numba has parallel feature, so ThreadPoolExecutor is not needed
    from concurrent.futures import ThreadPoolExecutor, as_completed
    with ThreadPoolExecutor(12) as ex:  # THIS IS NOT WORKING AS PARALLEL!!!
        if True:
            the_threads = [ex.submit(ACMSimPy, omega_ob) for omega_ob in range(1, 7)]
            for f in as_completed(the_threads):
                print(f.result())
        else:
            for f in ex.map(ACMSimPy, np.arange(1, 7, 1)):
                print(f.result())
    print('Done')




# %% https://www.youtube.com/watch?v=fKl2JW_qrso
import multiprocessing
print(multiprocessing.cpu_count())
print('Conclusion: numba.jitclass is likely to not support multiprocessing.')

processes = []
for omega_ob in range(1,7):
    p = multiprocessing.Process(target=ACMSimPy, args=[omega_ob])
    p.start()
    processes.append(p)

for p in processes:
    p.join()
