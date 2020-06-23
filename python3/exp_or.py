import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import helper_function as hf

### Definition of the class ExpOR containing all the objects (ORs, compounds)
### necessary for the simulation of OR stimulation.

### For efficiency reason, this script is called only when
### the flag OR_flag = True in set_parameter.py

n_or = 100

class OR:

    """
    Definition of the OR class containing all the characteristics about
    ORs of our experiments.
    """

    ### Contains the affinity of all the ORs with all the odorants.
    ### Current is max for an affinity of 0? 0 is not working yet
    ### 0.5 is a good value (for a single compound) to have IOR values between
    ### 0 and Imax_compound if the cOd is between 0 and 1.
    affin_dic = {
        "A" : {"0": 0.5, "1": 0.1, "2": 1e+10, "3": 1e+10, "4": 1e+10, "5": 1e+10},
        "B" : {"0": 1e+10, "1": 1e-10, "2": 1e-5, "3": 1e+10, "4": 1e+10, "5": 1e+10},
        "C" : {"0": 0.1, "1": 0.2, "2": 0.15, "3": 1e+10, "4": 1e+10, "5": 1e+10},
        "D" : {"0": 0.2, "1": 0.1, "2": 0.15, "3": 1e+10, "4": 1e+10, "5": 1e+10},
        "OR2W1": {"DNT": 0.1, "TNG": 0.1, "TNT": 10.},
        "OR1A1": {"DNT": 0.1, "TNG": 10., "TNT": 0.2},
        "OR5K1": {"DNT": 0.25, "TNG": 0.1, "TNT": 10.},
        "189-1": {"DNT": 0.1, "TNG": 0.4, "TNT":0.2}

    }

#"OR2W1": {"DNT": 0.1, "TNG": 0.1, "TNT": 10.},
#"OR1A1": {"DNT": 0.1, "TNG": 10., "TNT": 0.2},
#"OR5K1": {"DNT": 0.25, "TNG": 0.1, "TNT": 10.}

    ### Contains the Imax for each ORs.

    Imax_dic = {"A": 100., "B": 100., "C": 100., "D": 50.,
                "OR2W1": 100., "OR1A1": 50., "OR5K1": 100., "189-1": 50.}


    def __init__(self, NOr = 3, Nodor = 8, name_ORs = [], name_odors = []):


        # Contains the name of each OR.
        if len(name_ORs) > 0:
            self.name_ORs = name_ORs
        else:
            self.name_ORs = [chr(i) for i in range(ord('A'), NOr + ord('A'))]
        # Vector containing Imax for each ORs of 'name_OR's'.
        self.Imax = np.array([self.Imax_dic[OR] for OR in self.name_ORs])
        # Contain the name of each odor.
        if len(name_odors) > 0:
            self.name_odors = name_odors
        else:
            self.name_odors = [str(i) for i in range(Nodor)]

        # if (NOr != len(self.name_ORs)) or (Nodor != len(self.name_odors)):
        #     raise Exception("Number of ORs/odors does not fit the name")
        self.NOr = len(self.name_ORs)
        self.Nodor = len(self.name_odors)

        # Contains the affinity of all ORs for all odorants in a matrix with
        # NOr rows (for the ORs) and Nodor columns (for the odors).

        self.affin = fill_matrix_from_lists_if_in_dict(self.name_ORs,
                                                       self.name_odors,
                                                       self.affin_dic)
        self.affin_inv = 1. / self.affin

###########################
# HELPER FUNCTIONS FOR OR #
###########################

def fill_matrix_from_lists_if_in_dict(row_list, column_list, dic_of_dic):
    """
    For each element r_i of row_list, if dic_of_dic[r_i] exists, it look for
    each element r_j of column_list if dic_of_dic[r_i][c_j] exists. If it does
    we stock it in out_tmp[i][j].
    """
    Nrow = len(row_list)
    Ncolumn = len(column_list)
    out_tmp = []
    for i in range(Nrow):
        r_i = row_list[i]
        # print r_i
        try:
            dic_i = dic_of_dic[r_i]
            # If dic_i exists, let's continue, otherwise let's try dic_i+1
            # print dic_i
            r_i_list = []
            for j in range(Ncolumn):
                r_j = column_list[j]
                try:
                    r_i_list.append(dic_i[r_j])
                    # if dic_i[r_j] exists, let's continue, otherwise let's
                    # try dic_i[r_j+1]
                    # print dic_i[r_j]
                except:
                    r_i_list.append(np.inf)
            out_tmp.append(r_i_list)
        except:
            raise Exception("[defOR]: Characteristics of '" + str(r_i) +
                            "' are not defined.")
    return np.array(out_tmp)



class ExpOR:
    """
    1) IOR = Inmax * ([X]/Knx + [Y]/Kny + ...)/(1 + [X]/Knx + [Y]/Kny + ...)
    """

    def __init__(self, dt, expParameters, Nneuron=n_or, rCoupl=0.2):

        self.param = expParameters
        self.name = expParameters["name"]
        self.stimSeq = expParameters["stimSeq"]
        self.bin_size = expParameters["bin_size"]

        self.ORs = expParameters["ORs"]
        self.Nodor = self.ORs.Nodor
        self.NOr = self.ORs.NOr


        self.cOd = stimSeq2mat(self.stimSeq, self.Nodor, self.bin_size)
        # Number of 'bins' of the experiment.
        self.T = self.cOd.shape[1]

        ### Contains 'sum_i [Xi]/Ki'. (dim = NOr x T)
        self.tmp = np.dot(self.ORs.affin_inv, self.cOd)

        ### Contains intensities produced by of all receptors (dim = NOr x T)
        Imax_tmp = self.ORs.Imax
        Imax_tmp = Imax_tmp.reshape((len(Imax_tmp), 1)) * np.ones((1, self.T))
        self.IOR = Imax_tmp / (1 + 1 / self.tmp)

        ### Linking the experiment to a Neuron population.
        self.link_with_popu(Nneuron, rCoupl)


    def link_with_popu(self, Nneuron, rCoupl):

        ### Which neuron express which OR in which proportion (0 to 1).
        ### Dim = (Nneuron x NOr)
        ### For now, the expression is binary (a neuron express or not an OR).
        ### XXX Change rCoupl into a vector (dim = NOr) so we can have
        ### different concentrations of each OR express.
        self.expression_OR = np.random.choice([0, 1],
                                              size=(Nneuron, self.NOr),
                                              p=[1 - rCoupl, rCoupl])

        ### Contains the ORs' intensities received by all the neurons.
        ### dim = (Nneuron x T)
        self.IOR_n = np.dot(self.expression_OR, self.IOR)


    def save(self, dirPath):
        """
        Save different parameter of the experiment in dirPath.
        """
        filePath = os.path.join(dirPath, "ORparameters.csv")
        tmp = self.param
        del tmp["stimSeq"]
        del tmp["ORs"]
        tmpdf = pd.DataFrame(tmp, index = [0])
        Nodor = self.ORs.Nodor
        tmpdf["Nodor"] = Nodor
        tmpdf["NOr"] = self.ORs.NOr
        tmpdf.to_csv(filePath, index = False)
        ### Write the Concentration of odorants.
        filePath = os.path.join(dirPath, "cOd.csv")
        self.cOd.tofile(filePath, sep = ",")
        ### Plot the concentration of odorants.
        plt.close()
        cOdPlot = self.cOd.transpose()
        t_sample = 12000
        for comp in range(Nodor):
            hf.plot_values(cOdPlot[:t_sample, comp], label = comp,
                           show_flag = False)
        plt.legend()
        filePath = os.path.join(dirPath, "cOd.png")
        #plt.savefig(dirPath + "cOd.png")
        plt.savefig(filePath)





def stim2vec(stim, Nodor, bin_size):
    """
    For now, this version handle only compounds name include in [|0,9|].
    """
    cOd = np.zeros((Nodor, bin_size))

    try:
        # If the stim start by a string it means that we want a continuous
        # stimulation of all the integers mentionned.
        int(stim[0])
        for c in stim:
        # if c is an integer change the int(C) value of the vector by 1.
            idx = int(c)
            cOd[idx,:] = 1.
    except:
        # If it start by a letter, we want either a baseline, either to apply
        # a function on all the integers mentionned afterward.
        if stim == 'bl':
            ### At baseline, concentration stays at 0.
            cOd = cOd
        else:
            if stim[0] == '+':
                # Linear increase on the bin for every compounds mentionned.
                cOd_tmp = np.linspace(0, 1, bin_size)
            elif stim[0] == '-':
                # Linear decrease on the bin for every compounds mentionned.
                cOd_tmp = np.linspace(1, 0, bin_size)
            # elif stim[0] == 's':
                # Create a sinusoid on the bin for every compounds mentionned.
            for c in stim[1:]:
                idx = int(c)
                cOd[idx,:] = cOd_tmp


    return cOd


def stimSeq2mat(stimSeq, Nodor, bin_size):
    tmp_list = [stim2vec(stim, Nodor, bin_size) for stim in stimSeq]
    return np.concatenate(tmp_list, 1)



#########
# DRAFT #
#########

# class ExpOR:
#     """
#     1) IOR = Inmax * ([X]/Knx + [Y]/Kny + ...)/(1 + [X]/Knx + [Y]/Kny + ...)
#     """

#     def __init__(self, stimSeq, Nneuron, rCoupl = 0.2):

#         if stimSeq == oneOd_oneOR:
#             """
#             First experience where we stimulate neurons with OR 'A' with
#             odorant '1'. Around 20% of the neurons express 'A'
#             """
#             self.Nodor = 1
#             self.NOr = 1
#             ### Initialization of the ORs.
#             # self.ORs = OR(name_ORs = ["A"], name_odors = ["1"])
#             self.ORs = OR(Nodor = self.Nodor, NOr = self.NOr)


#         if stimSeq == oneOd_oneOR_linVar:
#             """
#             First experience where we stimulate neurons with OR 'A' with
#             odorant '1' that concentration varies as a sinusoid.
#             Around 20% of the neurons express 'A'
#             """


#         if stimSeq == twoOd_twoOR:
#             """
#             First experience where we stimulate neurons with OR 'A' and 'B'
#             with odorant '1' and '2'. Around 20% of the neurons express 'A' or
#             'B'. Note: Few neurons will exprime both 'A' and 'B' (4%).
#             """

#             self.Nodor = 2
#             self.NOr = 2
#             ### Initialization of the ORs.
#             # self.ORs = OR(name_ORs = ["A"], name_odors = ["1"])
#             self.ORs = OR(Nodor = self.Nodor, NOr = self.NOr)

#         ### Initialization of the odors (Nodor x T).
#         # self.cOd_list = []
#         # XXX Replace stim2vec by a function stimSeq2Vec.
#         # self.cOd_list = [stim2vec(stim, self.Nodor) for stim in stimSeq]
#         # self.cOd = (np.array(self.cOd_list)).transpose()

#         self.cOd = stimSeq2mat(stimSeq, self.Nodor)
#         # Number of 'bins' of the experiment.
#         self.T = self.cOd.shape[1]

#         ### Contains 'sum_i [Xi]/Ki'. (dim = NOr x T)
#         self.tmp = np.dot(self.ORs.affin_inv, self.cOd)

#         ### Contains intensities produced by of all receptors (dim = NOr x T)
#         Imax_tmp = self.ORs.Imax
#         Imax_tmp = Imax_tmp.reshape((len(Imax_tmp), 1)) * np.ones((1, T))
#         self.IOR = Imax_tmp / (1 + 1 / self.tmp)

#         ### Linking the experiment to a Neuron population.
#         self.link_with_popu(Nneuron, rCoupl)


#     def link_with_popu(self, Nneuron, rCoupl):

#         ### Which neuron express which OR in which proportion (0 to 1).
#         ### Dim = (Nneuron x NOr)
#         ### For now, the expression is binary (a neuron express or not an OR).
#         ### XXX Change rCoupl into a vector (dim = NOr) so we can have
#         ### different concentrations of each OR express.
#         self.expression_OR = np.random.choice([0, 1],
#                                               size=(Nneuron, self.NOr),
#                                               p=[1 - rCoupl, rCoupl])

#         ### Contains the ORs' intensities received by all the neurons.
#         ### dim = (Nneuron x T)
#         self.IOR_n = np.dot(self.expression_OR, self.IOR)



# def stimSeqToExp(stimSeq, N, rCoupl = 0.2):



#     def stim2vec(stim, Nodor):
#         if stim == 'a':
#             return np.zeros(Nodor)
#         if stim == 'b':
#             return np.ones(Nodor)


#     if stimSeq == oneOd_oneOR:
#         """
#         First experience where we stimulate neurons with OR 'A' with
#         odorant '1'. Around 20% of the neurons express 'A'
#         """
#     if stimSeq == twoOd_twoOR:
#         """
#         First experience where we stimulate neurons with OR 'A' and 'B'
#         with odorant '1' and '2'. Around 20% of the neurons express 'A' or
#         'B'. Note: Few neurons will exprime both 'A' and 'B' (4%).
#         """

#             self.Nodor = 2
#             self.NOr = 2
#             ### Initialization of the ORs.
#             # self.ORs = OR(name_ORs = ["A"], name_odors = ["1"])
#             self.ORs = OR(Nodor = self.Nodor, NOr = self.NOr)

#         ### Initialization of the odors (Nodor x T).
#         # self.cOd_list = []
#         self.cOd_list = [stim2vec(stim, self.Nodor) for stim in stimSeq]
