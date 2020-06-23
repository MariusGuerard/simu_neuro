############
# PACKAGES #
############

# General Packages.
import time
import shutil

import popu


###################################################
# TMP before deleting/libraryzing 'parameter.py'. #
###################################################

DT = 0.2
BIN_SIZE = 250
EXP_SR = [1, 1] * 2

OUT_DIR = "../output/"

multithread_flag = False

##############
# SIMULATION #
##############

if multithread_flag == False:

    ### Serial simulation.
    # Network Initialization.
    pop = popu.Popu()
    pop.connectGaussian()

    # Simulation.
    ### replace time.time by decorator.
    ts1 = time.time()
    # s1 = pop.simul(stim_list=EXP_SR, bin_size_t=BIN_SIZE_T, dt=DT,
    #                write_output=True, save_fig=True,
    #                name_simul = "dt_" + str(DT), out_path=OUT_DIR)
    s1 = pop.simul(stim_list=EXP_SR, bin_size=BIN_SIZE, dt=DT,
                   write_output=True, save_fig=True,
                   name_simul = "dt_" + str(DT), out_path=OUT_DIR)

    ts2 = time.time()
    print("time = {}".format(ts2 - ts1))

else:
    ### Parallel computing.
    ### Needs to be executed with 'python init.py' from the terminal.
    # Device Initialization
    num_popu = 20
    NUM_CORE = num_popu
    D1 = device.Device(num_popu)

    # Simulation.
    exp_name = "test_mp"
    for i in range(num_popu):
        # D1.exp[exp_name + str(i)] = exp_mp
        D1.exp[exp_name + str(i)] = exp_sr
    print("Parallel simulation of {} ms".format(length_simu * num_popu))
    ts1 = time.time()
    s = D1.simul_mp(NUM_CORE, tsimu)
    ts2 = time.time()

    print("time = {}".format(ts2 - ts1))


#shutil.rmtree('../output/')
#########
# DRAFT #
#########
