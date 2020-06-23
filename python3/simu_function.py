import numpy as np
import matplotlib.pyplot as plt
image_dir = '.'

#################
# VISUALIZATION #
#################

def plot_Vm(Vm_vec, save_path = "./Vm_last_simu.png"):
    plt.close()
    Ntimesteps = len(Vm_vec)
    Nneuron = Vm_vec.shape[1]
    for i in range(Nneuron):
        plot_values(Vm_vec[:,i], show_flag = False, label = str(i))
    plt.legend()
    plt.savefig(save_path)

# Create the Graph of the Network.
def draw_graph_connection(connection_vec, Ninhib, on_vec,
                          position_flag = False, position_vec = np.array([]),
                          save_fig = False, name_fig = "0"):
    plt.close()
    N = len(connection_vec)
    G = nx.Graph()
    ### Neurons.
    # Add neurons.
    G.add_nodes_from(list(range(N)))
    # Green represents excitatory neuron, red inhibitory (only when they fire).
    # When they don't fire, neuron are blue)
    # color_list = ['g' for i in range(N - Ninhib)] \
    #              + ['r' for i in range(N - Ninhib, N)]
    color_list = ['b' for i in range(N)]
    for i in np.where(on_vec)[0]:
        if i < N - Ninhib: # ! Careful to the limit case (<= better?)
            color_list[i] = 'g'
        else:
            color_list[i] = 'r'
    # Add labels (used for understanding small networks)

    # Add synapses.
    w = np.where(connection_vec == 1)
    G.add_edges_from([(w[0][i], w[1][i])
                      for i in range(int(connection_vec.sum()))])

    # pos to be changed when we affect one to neurons.
    if position_flag:
        pos = dict()
        for i in range(N):
            #     G.node[i]['pos'] = self.n[i].position[:2]
            # pos[i] = self.n[i].position[:2]
            pos[i] = position_vec[i][:2]
    else:
        pos = nx.random_layout(G)
    # draw graph
    nx.draw(G, pos, node_size = 100, node_color = color_list)

    # show graph
    plt.show()
    if save_fig:
        plt.savefig(name_fig + ".png")
    plt.close()
    return G


# Functions to generate visualization from saved simulations.

def anim_from_on(on_vec, connection_vec, Ninhib,
                 position_flag = True, position_vec = np.array([]),
                 save_flag = True, save_path = "./tmp/"):
    Ntimesteps = on_vec.shape[0]
    Nneurons = on_vec.shape[1]

    for t in range(Ntimesteps):
        on_t = on_vec[t,:]
        draw_graph_connection(connection_vec, Ninhib, on_t,
                              position_flag = position_flag, position_vec = position_vec,
                              save_fig = save_flag, name_fig = save_path + str(t))


#def Vm_from_V(V_vec,


def divide_and_scatter(on_vec_long, save_path, new_length, dt, inhib_list = []):
    # In case of long simulation, we plot it every 'new_length' time steps.
    Ntimesteps = len(on_vec_long)
    if Ntimesteps <= new_length:
        scatter_from_on(on_vec_long, dt=dt, inhib_list = inhib_list,
                        save_path = save_path)
    else:
        Nscatter = int(Ntimesteps / new_length)
        print(Ntimesteps)
        print(Nscatter)
        for i in range(Nscatter):
            scatter_from_on(on_vec_long[i * new_length: (i + 1) * new_length],
                            dt=dt,
                            origin = new_length * i * dt,
                            inhib_list = inhib_list,
                            save_path = save_path + "_" + str(i))
            # To generate the last scatter plot if tsimu % 1000 != 0
            if ((Ntimesteps % new_length) != 0):
                scatter_from_on(on_vec_long[Nscatter * new_length:],
                                dt=dt,
                                origin = new_length * Nscatter * dt,
                                inhib_list = inhib_list,
                                save_path = save_path + "_" + str(i + 1))


# def scatter_from_on(on_vec, origin = 0., inhib_list = [],
#                     show_flag = False, save_flag = True,
#                     save_path = image_dir + "scatter_last_simu",
#                     size_marker = 10., size_fig = 10.):

def scatter_from_on(on_vec, save_path, dt,
                    origin=0., inhib_list=[],
                    show_flag = False, save_flag = True,
                    size_marker = 10., size_fig = 10.):
    """
    Main function for plotting scatter plots after simulations.
    Add: Plot x as time in ms and not timestep.
    Add: Color red for inhibitor neurons.
    """

    Ntimesteps = len(on_vec)
    Nneurons = on_vec.shape[1]
    size_simu = float(Ntimesteps)
    if save_flag:
        plt.close('all')
        s_x = size_fig * (1 + 0.001 * size_simu)
        s_y = size_fig + (1 + 0.001 * Nneurons)
        plt.figure(figsize=(s_x, s_y))
    for t in range(Ntimesteps):
        who_fire = np.where(on_vec[t])[0]
        t_vec = np.ones(len(who_fire)) * t * dt + origin
        # plt.plot(t_vec, num_fire, marker='o', linestyle='none',
        #          color='b', markersize = 1)
        plt.plot(t_vec, who_fire, marker='o', linestyle='none',
                 color='b', markersize = 100. * size_marker / size_simu)

        # Note: Might be better to plot only excitator neurons on the above.
        # Add red color for inhibitor neurons.
        who_fire_inhib = np.intersect1d(inhib_list, who_fire, True)
        t_vec = np.ones(len(who_fire_inhib)) * t * dt + origin
        plt.plot(t_vec, who_fire_inhib, marker='o', linestyle='none',
                 color='r', markersize = 100. * size_marker / size_simu)

        plt.xlabel("t(ms)")
        plt.ylabel("Neuron label")
        # plt.xlim((0, Ntimesteps))
        # plt.ylim((0, Nneuron))
    if show_flag:
        plt.show()
    if save_flag:
        plt.savefig(save_path + ".png", dpi = 200)
        # plt.savefig(save_path + ".svgz", format='svgz')
        # plt.savefig(save_path + ".eps", format='eps')
        # plt.savefig(save_path + ".svg", format='svg')



##########
# PATENT #
##########

def I_Rn(conc_vec, kn_vec, In_max):
    """
    Implement Renaud formula:

    """

    tmp = conc_vec/kn_vec
    return In_max * tmp / (1 + tmp)


####### DRAFT
