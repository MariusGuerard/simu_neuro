### Definition of the class 'Device' that can contain different 'Popu'.

# Multithreading.
import multiprocessing as mp
import copy


import helper_function as hf
import popu



class Device:

    ##################
    # INITIALIZATION #
    ##################

    def __init__(self, Npop = 0):

        self.Npop = Npop

        # 'pop' contains the neuron populations.
        self.pop = []

        # 'exp' contains the experience description for each populations.
        self.exp = {}

        # For now, if user specify Npop at the initialization, we fill the
        # Device with 'Npop' hard copy of 'pop0'.
        # Notably useful for parallel computing.
        if Npop != 0:
            # pop0 = popu.Popu(N, ri = ratio_inhib)
            pop0 = popu.Popu()
            # pop.build_square(1,1,2,2)
            # pop0.connectGaussian(ratio_connec)
            pop0.connectGaussian()
            self.pop.append(pop0)
            if Npop > 1:
                for p_i in range(Npop - 1):
                    self.pop.append(copy.copy(pop0))


    def simul(self):
        for i in range(self.Npop):
            pop_i = self.pop[i]
            name_i = list(self.exp.keys())[i]
            exp_i = self.exp[name_i]
            pop_i.simul(tsimu, exp_i, ['on'],
                        name_simul = name_i)



    def mp_genpara(self, tsimu):
        """
        Create the variable called by pool.map for multiprocessing and store
        it into 'self.sim_arg'.

        Note: def simul(self, Ntimesteps, stim_list, outVar = ['on'],
        time_flag = True, write_output = True, name_simul = "tmp", save_fig = True):

        To set time_flag = [False] * (Npop - 1) allows to display only
        the time of one of the process.
        """
        Npop = self.Npop
        tsimu_l = [tsimu] * Npop
        stim_l = list(self.exp.values())
        outVar_l = [['on']] * Npop
        time_flag_l = [True] + [False] * (Npop - 1)
        write_output_l = [True] * Npop
        name_simul_l = list(self.exp.keys())
        save_fig_l = [False] * Npop
        self.sim_arg = list(zip(*(self.pop, tsimu_l, stim_l, outVar_l, time_flag_l,
                             write_output_l, name_simul_l, save_fig_l)))


    def simul_mp(self, num_core, tsimu_):
        """
        ! Must be used outside of ipython if matplotlib is used
        (with e.g. python init.py from terminal)!
        """

        self.mp_genpara(tsimu_)
        pool = mp.Pool(num_core)
        result_list = pool.map(hf.worker, [simu for simu in self.sim_arg])
        pool.close()
        pool.join()
        return result_list

        # pool = ProcessingPool(nodes=Npop)
        # pool.map(pop_i, [1,2,3], [1,1,1])
