"""
Definition of a Synapse.

    A synapse is represented by:
    - its post-synaptic neuron: the 'target'.
    - its pre-synaptic neuron: the 'origin'.
    - its weight.
    - its delays.
    - its post-synaptic current PSC.

    - its neurotransmitter pool.
    - its STDP parameters.

Author: Marius Guerard
"""

import numpy as np

### For basic synapses.
PSC_BASIC = 300.

### For Renaud synapses.
# Weight for excitatory synapses.
W_E = 40.
# Weight for inhibitory synapses.
W_I = 40.



class Basic:
    """
    """

    def __init__(self, target_number, origin_number, typ, psc=PSC_BASIC):
        """
        Careful, for now all the delay are 0 because all the neurons have
        the same position (0, 0, 0).
        """

        # Permanent parameters.
        self.target_n = target_number
        self.origin_n = origin_number
        self.name = str(origin_number) + "-" + str(target_number)
        self.delay = 0
        self.typ = typ

        ### Dynamic parameters (likely to change during the simulation).

        # Post Synaptic Current.
        self.PSC = 0.


    def update(self):
        """
        In this first simple version of Synapses, the update is only a reset to 0 function.
        """
        self.PSC = 0.


    def fire(self):
        """
        In this first simple version of Synapses, the fire is only a constant.
        """
        #print "FIRE " + self.name
        self.PSC = psc * self.typ


class Renaud:
    """
    """

    def __init__(self, target_number, origin_number, typ, w_e=W_E, w_i=W_I):
        """
        Careful, for now the name of the synapse is not unique.
        All synapses between the same neurons hold the same name.

        Careful, for now all the delay are 0 because all the neurons have
        the same position (0, 0, 0).
        """

        # Permanent parameters.
        self.target_n = target_number
        self.origin_n = origin_number
        self.name = str(origin_number) + "-" + str(target_number)
        self.delay = 0
        self.typ = typ

        ### Dynamic parameters (likely to change during the simulation).

        # Post Synaptic Current.
        self.PSC = 0.
        # Quantity of neurotransmitter (nt) in the pool.
        self.P = 0.
        # Fraction of nt delivered by the pool after a spike.
        # Rename r_nt
        self.U = 0.5
        # self.U = 0.9
        # PSC of the last spike: rename PSC_last
        self.truc = 0
        # Time Since Last Spike.
        self.tsls = 0.
        # Characteristic time of the pool filling. Rename tauF.
        self.F = 30.

        # Weight
        if self.typ == 1:
            self.W = w_e
        else:
            self.W = w_i
        # Characteristic time of PSC decay in a synapse. Rename tauPSC.
        self.tau = 2
        # self.tau = 10.


    def update(self, dt):
        """
        Called everytime a neuron with this output synapse is updated.
        """
        self.PSC = self.typ * self.truc * np.exp(-self.tsls/self.tau)
        self.tsls += dt

    def fire(self):
        """
        Called everytime a neuron with this output synapse fires.
        """
        # print "FIRE " + self.name
        # Update the nt quantity in the pool (increase).
        self.P += (1 - self.P) * (1 - np.exp(-self.tsls/self.F))
        # Quantity of neurotransmitter delivered by the pool. Rename DP
        u = self.P * self.U
        # Record the PSC of the spike.
        self.truc = self.PSC + self.W * u
        # Update the PSC of the synapse with the PSC of the spike.
        self.PSC = self.truc
        # Update the nt quantity in the pool (decrease).
        self.P -= u
        # Update tsls.
        self.tsls = 0
