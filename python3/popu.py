### POPULATION

import os
import shutil

import pandas as pd
import numpy as np
from random import gauss, randint

import neuron
import helper_function as hf
import simu_function as sf


# Default Number of neurons in the population.
N = 100

# Nsel is the number of neurons in which we input a fixed current at each step.
# (should be in the experiment module).
rin_sel = 0.1
Nsel = int(rin_sel * N)

# Ratio of inhibitory neurons in the population.
ratio_inhib = 0.2

# Each neurons is in average connected to 'ratio_connec' * N other neurons.
ratio_connec = 0.2

# Ratio of neurons in which we input a random current.
rin = 0.1

# Distribution of the random current.
param_noise = {'mean_noise': 70, 'sigma_noise': 30}

# Choice of the neuron used by default('izhi' or 'lif').
neuron_default = 'izhi'


### Definition of the class 'Popu' representing a Neurons' population.
class Popu:


    ##################
    # INITIALIZATION #
    ##################

    def __init__(self, ri=ratio_inhib, Nneuron=N, model=neuron_default):

        ### Total number of neurons.
        self.Nneuron = Nneuron

        # Number of inhibitory and excitatory neurons defined by 'ratio_inhib'
        # Note: Most commonly, 'ratio_inhib' = 0.2 in the litterature.
        self.Ninhib = int(Nneuron * ri)
        self.Nexcit = Nneuron - self.Ninhib

        ### 'n' is a list containing the populations' neurons.
        self.n = []

        ### If the model of Population is defined as 'izhi' (For Izhikevich),
        # create Nexcit random neurons as regular spiking excitatory neurons
        # and Ninhib random neurons as fast spiking inhibitory neurons.
        n_permut = np.random.permutation(list(range(Nneuron)))
        n_inhib = n_permut[:self.Ninhib]
        self.n_inhib = n_inhib

        if model == 'izhi':
            for i in range(Nneuron):
                if i in n_inhib:
                    self.n.append(neuron.Izhikevich(i, typ = -1, spike_typ = 'fs'))
                else:
                    self.n.append(neuron.Izhikevich(i))

        else:
            raise NameError("This model has not been implemented yet.")


        ### 'on' is a np.array representing which neuron are on (=1) or off (=0).
        self.on = np.zeros(Nneuron)

        ### 'connection' is a np.array representing the connections between neurons.
        # connection[i,j] = 1 if neuron i is pre_synaptic of neuron j,
        # 0 otherwise.
        self.connection = np.zeros((Nneuron, Nneuron))


    def connectGaussian(self, meanConnectionRate = 0.1,
                        sigmaConnectionRate = 0.01):
       """
       Connect neurons of this population with each others, every random being
       connected to 'meanConnectionRate' other neurons on average
       (+- 'sigmaConnectionRate'). Note: Izhi chose mean = 0.1
       """
       Nneuron = self.Nneuron
       meanConnection = Nneuron * meanConnectionRate
       sigmaConnection = Nneuron * sigmaConnectionRate
       for i in range(Nneuron):
           Nconnection_neuron = int(gauss(meanConnection, sigmaConnection))
           for j in range(Nconnection_neuron):
               r = randint(0, Nneuron - 1)
               if r != i:
                   self.connect_i_j(i, r)


    def connect_i_j(self, i, j):
        """
        Connect neuron 'i' of the population with neuron 'j'.
        """
        self.n[i].connectTo(self.n[j])
        self.connection[i, j] = 1

    ##############
    # SIMULATION #
    ##############



    def update_neurons(self, stim, dt, mean_noise, sigma_noise, Iother):
        """
        Update all the neurons.
        Note: This function has been updated: Iother is a vector of size N,
        indicating how much external current should be added to each neurons
        at this time step (used for ORs for example).
        """
        Nneuron = self.Nneuron

        # Baseline stimulation.
        # Nin = int(Nneuron * rin)
        # neuron_random = np.random.permutation(range(Nneuron))
        # neuron_in = neuron_random[:Nin]
        neuron_in = hf.permut(Nneuron, rin)


        for ni in range(Nneuron):
            # Intensity noise received by the neurons.
            Iext = 0
            Inoise = gauss(mean_noise, sigma_noise)
            if stim == 0:
                """
                [n0 Experiment] We are putting random current only on n0.
                """
                if ni == 0: Iext = Inoise
            if stim == 1:
                """
                [Standard Experiment] Random current on a random subset of
                neurons of size Nin only. Note: If rin is chosen well, this
                should simulate a spontaneous burst in the population.
                """
                if ni in neuron_in: Iext = Inoise

            if stim == 2:
                """
                [Narrow Experiment 1] Random current on the 50 last neurons
                in addition of the random subset.
                """
                if ni in neuron_in: Iext = Inoise
                if (ni > Nneuron - Nsel):
                    Iext =  Inoise

            if stim == 3:
                """
                [Narrow Experiment 2] Random current on the 50 first neurons
                in addition of the random subset.
                """
                if ni in neuron_in: Iext = Inoise
                if (ni <= Nsel): Iext = Inoise

            if stim == 4:
                """
                [Intermediate Experiment 1] The 50 first neurons receive 20%
                of the random current and the the 50 last neurons receive 80%
                of the random current. (in addition of the random subset).
                it describes the case of broad receptors recevivng odorant A,
                and where the 50 first neurons have 20% affinity with receptor
                A and the 50 last neurons 80% affinity with receptor A.
                """
                if ni in neuron_in: Iext = Inoise
                if (ni <= Nsel): Iext = gauss(0.2 * mean_noise, sigma_noise)
                if (ni > Nneuron - Nsel):
                    Iext =  gauss(0.8 * mean_noise, sigma_noise)

            if stim == 5:
                """
                [Intermediate Experiment 2] The 50 first neurons receive 80%
                of the random current and the the 50 last neurons receive 20%
                of the random current. (in addition of the random subset).
                it describes the case of broad receptors recevivng odorant B,
                and where the 50 first neurons have 80% affinity with receptor
                B and the 50 last neurons 20% affinity with receptor B.

                """
                if ni in neuron_in: Iext = Inoise
                if (ni <= Nsel): Iext = gauss(0.8 * mean_noise, sigma_noise)
                if (ni > Nneuron - Nsel):
                    Iext =  gauss(0.2 * mean_noise, sigma_noise)

            if stim == 6:
                """
                [Broad Experiment 1] The 50 first neurons receive 45%
                of the random current and the the 50 last neurons receive 55%
                of the random current. (in addition of the random subset).
                it describes the case of broad receptors recevivng odorant A,
                and where the 50 first neurons have 45% affinity with receptor
                A and the 50 last neurons 55% affinity with receptor A.
                """
                if ni in neuron_in: Iext = Inoise
                if (ni <= Nsel): Iext = gauss(0.45 * mean_noise, sigma_noise)
                if (ni > Nneuron - Nsel):
                    Iext =  gauss(0.55 * mean_noise, sigma_noise)


            if stim == 7:
                """
                [Broad Experiment 2] The 50 first neurons receive 55%
                of the random current and the the 50 last neurons receive 45%
                of the random current. (in addition of the random subset).
                it describes the case of broad receptors recevivng odorant B,
                and where the 50 first neurons have 55% affinity with receptor
                B and the 50 last neurons 45% affinity with receptor B.

                """
                if ni in neuron_in: Iext = Inoise
                if (ni <= Nsel): Iext = gauss(0.55 * mean_noise, sigma_noise)
                if (ni > Nneuron - Nsel):
                    Iext =  gauss(0.45 * mean_noise, sigma_noise)

            if type(stim) == str:
                """
                [OR experiment] Stimulation of the neurons thanks to ORs
                (in addition of the random subset).

                """
                Iext = Iother[ni]
                if ni in neuron_in: Iext += Inoise

            self.n[ni].update(dt, Iext)
            self.on[ni] = self.n[ni].on


    def update_synapses(self, dt):
        """
        Update all the synapses.
        """
        for ni in range(self.Nneuron):
            # Reset all the synapses to 0.
            self.n[ni].update_output_synapses(dt)
            # if the neuron with index 'ni' has fired in the last update
            # fire its synapses.
            if self.on[ni]:
                for syn in self.n[ni].outSyn:
                    syn.fire()

    def update(self, stim, dt, mean_noise,
               sigma_noise, Iother = np.zeros(0)):
        """
        Update all the neurons of the network in an asynchronous way.
        # First loop: update the states of the neurons.
        # Second loop: update the states of the synapses.
        """
        self.update_neurons(stim, dt, mean_noise, sigma_noise, Iother)
        self.update_synapses(dt)




    ###############
    # SIMULATIONS #
    ###############



    def write_output(self, name_simul, stim_list, out_dic, dt, out_path):

        """
        ! Careful, if a directory is already named 'name_simul',
        it will be deleted.
        Add: np.savetxt instead of np.tofile (keep the structure but might be
        heavier.
        """
        # dirPath = out_dir + name_simul
        dirPath = out_path + name_simul
        if name_simul == "tmp":
            create_or_replace_directory(dirPath)
        else:
            os.makedirs(dirPath)
        for key in list(out_dic.keys()):
            filePath = os.path.join(dirPath, key + ".csv")
            out_dic[key].tofile(filePath, sep = ",")
            # np.savetxt(filePath, out_dic[key].astype(int), fmt = "%d",
            #            delimiter = ',')
            # var_to_save = pd.DataFrame(out_dic[key])
            # var_to_save.to_csv(filePath, index = False)
        # Write connection (e.g. use to plot graph.)
        filePath = os.path.join(dirPath, "connection.csv")
        self.connection.tofile(filePath, sep = ",")
        # Write inhibitors neurons (If not the last one).
        filePath = os.path.join(dirPath, "inhibitors.csv")
        self.n_inhib.tofile(filePath, sep = ",")
        # Write position (e.g. use to plot graph).
        # Write some network and simulation parameters.
        filePath = os.path.join(dirPath, "parameter.csv")
        param_df = pd.DataFrame()
        param_df['dt'] = [dt]
        # param_df['length_simu'] = length_simu
        # param_df['Ntimestep'] = tsimu
        # param_df['Nneuron'] = [self.Nneuron]
        # param_df['Ninhib'] = [self.Ninhib]
        # param_df['Nexc'] = [self.Nexcit]
        # param_df['connectivity'] = [ratio_connec]
        # if synapse_flag == 'Renaud':
        #     param_df['Synapse_e' + synapse_flag] = W_RR_e
        #     param_df['Synapse_i' + synapse_flag] = W_RR_i
        # param_df.to_csv(filePath, sep = ",", index = False)
        ### The list of stimulations and the length of each stimulation.
        # Convert stim_list into a string unfortunately...
        filePath = os.path.join(dirPath, "experiment.csv")
        exp_df = pd.DataFrame()
        exp_df['experiment'] = stim_list
        #exp_df['stim_bin_size'] = [bin_size] * len(stim_list)
        exp_df.to_csv(filePath, sep = ",", index = False)


        ### Copy the scripts that gave this analysis.
        # Working directory.
        dirWork = os.getcwd()
        scriptFolderName = os.path.join(dirPath, 'history_scripts')
        shutil.copytree(dirWork, scriptFolderName)


        # output_vec.astype('int16').tofile("test.bin")
        # output_vec.astype('int16').tofile("test.txt", sep = ",")
        return param_df, dirPath

    def write_ORexp(self, ExpOR, dirPath):
        ### Write the parameters of the OR experiment
        ### (without the stimSeq already written in write_output).
        filePath = os.path.join(dirPath, "ORparameters.csv")
        tmp = ExpOR.param
        del tmp["stimSeq"]
        del tmp["ORs"]
        tmpdf = pd.DataFrame(tmp, index = [0])
        Nodor = ExpOR.ORs.Nodor
        tmpdf["Nodor"] = Nodor
        tmpdf["NOr"] = ExpOR.ORs.NOr
        tmpdf.to_csv(filePath, index = False)
        ### Write the Concentration of odorants (Nodor x timesteps).
        filePath = os.path.join(dirPath, "cOd.csv")
        ExpOR.cOd.tofile(filePath, sep = ",")
        ### Plot the concentration of odorants.
        plt.close()
        cOdPlot = ExpOR.cOd.transpose()
        t_sample = 3000
        for comp in range(Nodor):
            plot_values(cOdPlot[:t_sample, comp], label = comp,
                        show_flag = False)
        plt.legend()
        filePath = os.path.join(dirPath, "cOd.png")
        #plt.savefig(dirPath + "cOd.png")
        plt.savefig(filePath)





    def simul(self, stim_list, bin_size, dt, out_path=None, outVar=['on'],
              time_flag=True, write_output=True, name_simul="tmp",
              save_fig=True, write_OR=True, fig_window=2500, exp_or=None):
        """
        Simulate the population with the 'stim_list' experiment.
        Note: Create a warning if the length of experiment is not the same as
        the length of the simulation.
        """
        n_time_steps = int(len(stim_list) * bin_size)
        length_simu = n_time_steps * dt

        print("dt = {}".format(dt))
        print("Serial simulation of {} ms.".format(length_simu))
        print("Corresponding to {} timesteps".format(n_time_steps))

        # Warning if the length of the simulation is not the same as the
        # experiment.
        if len(stim_list) * bin_size < n_time_steps:
            print ("Warning [Popu.simul]: simu longer than the experiment")
        if len(stim_list) * bin_size > n_time_steps:
            print ("Warning [Popu.simul]: simu shorter than the experiment")

        # save_fig only if write_output and save_fig was initially True.
        save_fig *= write_output

        out_dic = dict()

        for key in outVar:
            out_dic[key] = np.zeros((n_time_steps, self.Nneuron))

        for t in range(n_time_steps):
            if time_flag:
                if t % 50 == 0:
                    print("t = " + str(t))
                # print "t = " + str(t)


            ### Fill all the variables specify by user.
            if 'on' in outVar:
                out_dic['on'][t, :] = self.on
            if 'Vm' in outVar:
                V_t = [self.n[i].Vm for i in range(self.Nneuron)]
                out_dic['Vm'][t, :] = np.round_(V_t, 2)
            if 'Iin' in outVar:
                I_t = [self.n[i].Iin for i in range(self.Nneuron)]
                out_dic['Iin'][t, :] = np.round_(I_t, 2)

            ### Simulate the network following the desired experiment.
            self.exp_gen(t, stim_list, dt=dt, exp_or=exp_or, bin_size=bin_size)

        # Note: It might be faster/lighter to store in binary but
        # from the test I made it was not significant.
        if write_output:
            print("Writing Output")
            # ! Careful, if a directory is already named 'name_simul',
            # it will be deleted.
            if exp_or == None:
                param_df, dirPath = self.write_output(name_simul,
                                                      stim_list, out_dic, dt,
                                                      out_path=out_path)
            else:
                param_df, dirPath = self.write_output(exp_or.name,
                                                      stim_list, out_dic, dt,
                                                      out_path=out_path)
                #if write_OR: exp_or.save(dirPath)




            if save_fig:
                dirPath += "/scatterplots/"
                os.makedirs(dirPath)
                scat_name = os.path.join(dirPath, "scatter")
                # print scat_name
                # divide_and_scatter(out_dic['on'], save_path = scat_name,
                #                    new_length = fig_window)
                sf.divide_and_scatter(out_dic['on'], dt=dt,
                                      save_path=scat_name,
                                      new_length=fig_window,
                                      inhib_list=self.n_inhib)
            return out_dic
        else:
            return out_dic

    ### EXPERIMENTS.



    ### Add-on: continuously varying (among one bin) current
    ### compatible with defExpOR.py.

    ### Possibility to have experiment length different from simulation length
    ### had been depreciated.

    def exp_gen(self, t, stim_list, dt, bin_size, exp_or=None):
        """
        Simulate any experiment.
        For example: exp_gen(t, [1, 2, 1, 3, 1], 500) = exp_narrow(t)
        Note: we have also exp_gen(t, [1, 2, 1, 3], 500) = exp_narrow(t)

        'stim_index_t' is the index of stim_list we are at time step t.
        For example if bin_size is 100, we are at t = 150, and our
        'stim_list' is [6,7,8], 'stim_index_t' will be 1 and stim_t will be 7.
        """

        stim_index_t = int(t / bin_size)
        stim_t = stim_list[stim_index_t]
        try:
            # if exp_or is given to exp_gen transmit it to update.
            IOR_n = exp_or.IOR_n
            self.update(dt=dt, stim=stim_t, Iother=IOR_n[:, t], **param_noise)
        except:
            self.update(dt=dt, stim=stim_t, **param_noise)

