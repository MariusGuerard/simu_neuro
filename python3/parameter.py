#from def_exp_or import *
import exp_or

# Script containing all the global parameters of the model.

### Simulation parameters
# 'multithread_flag' must be True to operate parallel computing.
# Note: for now multithreading should be ran from the terminal with:
# > python init.py
# There will be as many processor as they are of Nrepet.
multithread_flag = False

# 'or_flag' must be True to start a simulation where Neurons are stimulated by
# or-chemicals bindings.
or_flag = False

# Time-step length (Note: Izhi chose 1ms in some article...).
dt = 0.2
#dt = 1.


###############
# EXPERIMENTS #
###############

### Length of each stimulation.
### Note: Best if the total length of an experiment is equal to length_simu.
# Length (in ms) of each stimulation.
#bin_size_t = 200
bin_size_t = 50

# Length (in timesteps) of each bin.
bin_size = int(bin_size_t / dt)


exp_bl = [1, 1]
exp_exclusive = [1, 2, 1, 3, 1]
exp_narrow = [1, 4, 1, 5, 1]
exp_broad = [1, 6, 1, 7, 1]


### For OR.

oneOd_oneOR = {"name": "oneOd_oneOR",
               "stimSeq": ['bl', '0'],
               "bin_size_t": bin_size_t,
               "ORs": exp_or.OR(Nodor = 1, NOr = 1)}

# Many odor, one OR (compartimented device).

twoOd_oneOR_broad2x = {"name": "twoOd_oneOR_broad2x",
                       "bin_size_t": bin_size_t,
                       "stimSeq": ['0', 'bl', 'bl', '1'] * 20,
                       "ORs": exp_or.OR(name_ORs = ["C"],
                                 name_odors = ["0", "1"])}
twoOd_oneOR_broad2x = {"name": "twoOd_oneOR_broad2x",
                       "bin_size_t": bin_size_t,
                       "stimSeq": ['0', 'bl', '1', 'bl'] * 20,
                       "ORs": exp_or.OR(name_ORs = ["C"],
                                 name_odors = ["0", "1"])}


# AB-test: threeOd+BL_manyORs

abTesting0 = {"name": "abTesting0",
          "bin_size_t": bin_size_t,
          "stimSeq": ['0', 'bl', 'bl', '1', 'bl', '2', 'bl'] * 20,
          "ORs": exp_or.OR(name_ORs = ["C", "D"],
                    name_odors = ["0", "1", "2"])}

abTesting1 = {"name": "abTesting1",
              "bin_size_t": bin_size_t,
              "stimSeq": ['bl', '0', 'bl', '1', 'bl', '2'] * 20,
              "ORs": exp_or.OR(name_ORs = ["OR2W1", "OR1A1", "OR5K1"],
                    name_odors = ["DNT", "TNG", "TNT"])}



# Others.
twoOd_twoOR = {"name": "twoOd_twoOR",
               "bin_size_t": bin_size_t,
               "stimSeq": ['0', 'bl', 'bl', '1'] * 20,
               "ORs": exp_or.OR(Nodor = 2, NOr = 2)}

twoOd_twoOR_broad = {"name": "twoOd_twoOR_broad",
                     "bin_size_t": bin_size_t,
                     "stimSeq": ['0', 'bl', '1', 'bl'],
                     "ORs": exp_or.OR(name_ORs = ["C", "D"],
                               name_odors = ["0", "1"])}


threeOd_twoOR = {"name": "threeOd_twoOR",
                 "stimSeq": ['0', 'bl', 'bl', '1', 'bl', '2'],
                 "bin_size_t": bin_size_t,
                 "ORs": exp_or.OR(Nodor = 3, NOr = 2)}

oneOdInc_oneOR = {"name": "oneOdInc_oneOR",
                  "stimSeq": ['+0', 'bl', 'bl'] * 20,
                  "bin_size_t": bin_size_t,
                  "ORs": exp_or.OR(Nodor = 1, NOr = 1)}

twoOdInc_oneOR = {"name": "twoOdInc_oneOR",
                  "stimSeq": ['+0', 'bl', '+1', 'bl'] * 20,
                  "bin_size_t": bin_size_t * 2,
                  "ORs": exp_or.OR(Nodor = 2, NOr = 1)}


exp_or_param = oneOd_oneOR
#exp_or_param = twoOd_twoOR_broad


exp_or_param["bin_size"] = int(exp_or_param["bin_size_t"] / dt)

T = len(exp_or_param["stimSeq"])

### Experience to conduct for serial simulations.
# Number of repetition of a stimulation scheme.
Nrepet = 1
if or_flag:
    exp_sr = exp_or_param["stimSeq"]
    bin_size = exp_or_param["bin_size"]
else:
    Nrepet = 1
    exp_name = 'exp_broad'
    exp_sr = exp_bl * Nrepet


######################
# GENERAL PARAMETERS #
######################

# length_simu = len(exp_sr) * bin_size_t

# # Length of the simulation in time steps.
# tsimu = int(length_simu / dt)

# Length of the simulation in time steps.
tsimu = int(len(exp_sr) * bin_size)

# Length of the simulation in ms.
length_simu = tsimu * dt












#########
# DRAFT #
#########

### Network parameters.
# Number of neurons in the population.
#N = 100

# Ratio of inhibitory neurons in the population.
# ratio_inhib = 0.2
#ratio_inhib = 0.

# # Each neurons is in average connected to 'ratio_connec' * N other neurons.
# ratio_connec = 0.2

# Ratio of neurons in which we input a random current.
# rin = 0.01
# rin = 0.1

# Selective stimulation. Nsel represent the number of neuron we are
# selectively stimulating.
#rin_sel = 0.1
#Nsel = int(rin_sel * N)


############
# SYNAPSES #
############

# synapse_flag = 'basic'
#synapse_flag = 'Renaud'

# ### For Renaud synapses.
# # Weight for excitatory synapses.
# W_RR_e = 40.
# # W_RR_e = 60.
# # Weight for inhibitory synapses.
# W_RR_i = 40.
# #W_RR = 300.

### STDP (XXX To be checked).
# Aeep, Aeem = 0.05, -0.05
# Aeip, Aeim = 0.05, -0.05
# Aiep, Aiem = 0., 0.
# Aiip, Aiim = 0., 0.

# Weep = 2000
# Weip = 2000
# Wiep = 2000
# Wiip = 2000

###########
# NEURONS #
###########


### Intensity Noise.
#param_noise = {'mean_noise': 70, 'sigma_noise': 30}
