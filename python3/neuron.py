"""
Definition of the class coresponding to the different neurons' model.

Note: For now, a neuron is represented by the class Izhikevich.
Modifications in the other classes might be necessary to change the model.

Author: Marius Guerard
"""

import synapse

#synapse_flag = 'basic'
synapse_flag = 'Renaud'


class Izhikevich:
    """
    Model the spike dynamic of the neuron following Izhikevich (2003).

    A neuron can be represented by:
    - its type: it can be inhibitor or excitator.
    - its position (x, y, z) in um (not used for now).
    - its input synapses.
    - its output synapses.
    - its spikes parameters.
    """

    def __init__(self, number=0, typ=1, spike_typ='rs',
                 position=(0, 0, 0)):

        self.number = number

        ### Type of Neuron
        # Excitatory or Inhibitory (inhibitory: -1 ; excitatory: 1).
        self.typ = typ
        # Type of spiking (Fast spiking, ...).
        self.spike_typ = spike_typ

        ### Geometry.
        # Position is defined by a triplet (x, y, z).
        self.position = position

        ### Synapses.
        # List of input synapses of this neuron.
        self.inSyn = []
        self.outSyn = []

        ### Spike Parameters.
        # 'on' == True if the Neuron is currently firing.
        self.on = 0
        # 'Iin' is the input current in the neuron.
        self.Iin = 0.
        # Time since last spike (initialized to a big number).
        self.tsls = 10000
        # Number of spikes.
        self.nSpikes = 0


        ### Dynamic variables.
        self.Isyn = 0

        self.Vm = -65.
        self.Vm_buffer = 0.
        self.Vthresh = 30.

        self.u = 0.
        self.u_buffer = 0.



        # Parameters for an excitatory regular spiking neuron.
        if (self.typ, self.spike_typ) == (1, 'rs'):
            self.a = 0.02
            self.b = 0.2
            self.c = -65.
            self.Vreset = -65.
            self.d = 8.0

        # Parameters for an inhibitory fast spiking neuron.
        elif (self.typ, self.spike_typ) == (-1, 'fs'):
            self.a = 0.1
            self.b = 0.2
            self.c = -65.
            self.Vreset = -65.
            self.d = 2.0

        else:
            raise NameError("typ or spike_typ not defined.")


    ### Network construction.

    def connectTo(self, target_neuron):
        """
        If distance is time consuming, maybe compute the delay of all synapses
        in once, using for example scipy.spatial.distance.
        """
        if synapse_flag == 'basic':
            new_synapse = synapse.Basic(target_neuron.number, self.number, self.typ)
        else:
            new_synapse = synapse.Renaud(target_neuron.number, self.number, self.typ)
        self.outSyn.append(new_synapse)
        target_neuron.inSyn.append(new_synapse)


    ### Dynamics.

    def fire(self):
        """
        When the neuron fire, we update the variables 'Vm', 'u',
        'nSpikes', 'on', and 'tsls'.
        """
        self.nSpikes += 1
        self.on = 1
        self.Vm = self.Vreset
        self.u = self.u + self.d
        self.tsls = 0



    def update(self, dt, Iext):
        """
        Update the buffer of the neuron in function of the intensity delivered
        by the synapses (+ an external current Iext), and the Izhikevich model.
        """

        # If the membrane potential is above the threshold: Fire!
        if self.Vm > self.Vthresh:
            self.fire()

        # Otherwise, "load the soma".
        else:
            # If the neuron fired in the last time step, inactivate it now.
            if self.on:
                self.on = 0
            # Let's add a timestep to the time since last spike.
            self.tsls += dt

            ### Compute input intensity.
            # Compute the input intensity coming from the synapses.
            for syn in self.inSyn:
                self.Isyn += syn.PSC
            # Add a random gaussian noise.
            self.Iin = Iext + self.Isyn

            # Reset the intensity received by the synapses.
            self.Isyn = 0

            ### Update the neuron's parameters Vm and u according to Izhikevich
            ### model with updated value of Vm, u, a, b, Iin.
            self.update_v_u_izhi(dt)


    ### Helper functions.

    def update_v_u_izhi(self, dt):
        """ Return v and u at the next time step according to the
        diff. eq. of Izhikevich (2003).
        v' = 0.04v^2 + 5v + 140 -u + I
        u' = a(bv - u)
        """
        self.Vm = self.Vm + dt * ((0.04 * self.Vm + 5) * self.Vm + 140 - self.u + self.Iin)
        self.u = self.u + dt * (self.a * (self.b * self.Vm - self.u))



    def update_output_synapses(self, dt):
        """
        Update the all the output Synapses of the neuron.
        Update can be different according to the synapses' class chosen.
        """
        for syn in self.outSyn:
            syn.update(dt)
