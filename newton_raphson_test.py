import numpy as np

# TODO : DELETE START
from examples_example_input import parameter_set1
nodes,lines,transformers = parameter_set1
from network import Node, Line, Transformer, Network
from network_utils import check_network, form_inductance_array, print_inductance_array

nodes = [Node.from_str(node.strip())
         for node in nodes.strip().split("\n")]
transformers = [Transformer.from_str(transformer.strip())
                for transformer in transformers.strip().split("\n")]
transformers = [transformer for transformer in transformers
                if transformer is not None]
lines = [Line.from_str(line.strip())
         for line in lines.strip().split("\n")]
branches = transformers + lines

node_num = check_network(branches, nodes)
Y = form_inductance_array(branches)
# print_inductance_array(Y)
network = Network(nodes, Y)

eps = 1e-5
max_num_iter = 15
# TODO : DELETE END


# TODO revert START
# def NR(network,eps,max_num_iter):
'''
won't consider Y as a sparse matrix
:param max_num_iter: if iteration times > max_num_iter, return None
:param eps: iteration termination criteria, epsilon

:return: (1) node voltage(amplitude and phase angle) array if iteration converges,
         else None
         # TODO: do not return network
         (2) iteration times taken
'''
num_PQ_nodes = network.get_num_PQ_nodes()
num_PV_nodes = network.get_num_PV_nodes()
num_Vtheta_nodes = network.get_num_Vtheta_nodes()

for num_iter in range(max_num_iter):
    print("***************Starting num_iter: {}".format(num_iter))
    fx_stack = []
    Pix_stack = []
    Qix_stack = []
    for node_i in network.PQ_nodes + network.PV_nodes:
        Pi = 0.0
        adjacent_nodes= network.get_adjacent_nodes(node_i)
        for node_j in adjacent_nodes:
            Gij, Bij = network.get_Gij_Bij(node_i, node_j)
            thetaij = network.get_thetaij(node_i, node_j)
            aaa= np.cos(thetaij)
            aab= np.sin(thetaij)
            d_Pi = node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
            Pi += d_Pi
        delta_P = node_i.PG - node_i.PL - Pi
        # TODO: MAKE SURE REALLY UNDERSTAND INDEX
        fx_stack.append(delta_P)
        Pix_stack.append(Pi)

    for node_i in network.PQ_nodes:
        Qi = 0.0
        for node_j in network.get_adjacent_nodes(node_i):
            Gij, Bij = network.get_Gij_Bij(node_i, node_j)
            thetaij = network.get_thetaij(node_i, node_j)
            Qi += node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
        delta_Q = node_i.QG - node_i.QL - Qi
        fx_stack.append(delta_Q)

    for node_i in network.PQ_nodes + network.PV_nodes:
        Qi = 0.0
        for node_j in network.get_adjacent_nodes(node_i):
            Gij, Bij = network.get_Gij_Bij(node_i, node_j)
            thetaij = network.get_thetaij(node_i, node_j)
            Qi += node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
        Qix_stack.append(Qi)
    abs_fx_max = np.max(np.abs(fx_stack))
    if abs_fx_max <= eps:
        # TODO: REVERT
        # return network, num_iter
        exit(0)
    else:
        print("------------fx and max_of_abs_fx")
        print("{}".format(fx_stack))
        print("max of fx: {}".format(abs_fx_max))
        # calculate Jacobian matrix J
        J = np.zeros(dtype=float, shape=(num_PQ_nodes + num_PV_nodes + num_PQ_nodes,
                                         num_PQ_nodes + num_PV_nodes + num_PQ_nodes))
        '''
        caluculate J
        J: H N
           M L
        '''
        for i in range(num_PQ_nodes + num_PV_nodes):
            # H
            node_i = network.get_node(i)
            for j in range(num_PQ_nodes + num_PV_nodes):
                node_j = network.get_node(j)
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                if i == j:
                    J[i, j] = node_i.V ** 2 * Bij + Qix_stack[i]
                else:
                    J[i, j] = - node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
            # N
            for j in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
                node_j = network.get_node(j - (num_PQ_nodes + num_PV_nodes))
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                if i + (num_PQ_nodes + num_PV_nodes) == j:
                    J[i, j] = - node_i.V ** 2 * Gij - Pix_stack[i]
                else:
                    J[i, j] = - node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))

        for i in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
            # M
            node_i = network.get_node(i - (num_PQ_nodes + num_PV_nodes))
            for j in range(num_PQ_nodes + num_PV_nodes):
                node_j = network.get_node(j)
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                if i == j + num_PQ_nodes + num_PV_nodes:
                    J[i, j] = node_i.V ** 2 * Gij - Pix_stack[j]
                else:
                    J[i, j] = node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
            # L
            for j in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
                node_j = network.get_node(j - (num_PQ_nodes + num_PV_nodes))
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                if i == j:
                    J[i, j] = node_i.V ** 2 * Bij - Qix_stack[i-(num_PQ_nodes + num_PV_nodes)]
                else:
                    J[i, j] = - node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
        print("------------J")
        print(J)
        # update network parameters
        inv_J = np.linalg.inv(J)
        fx =np.reshape(fx_stack,(len(fx_stack),1))
        dx = np.matmul(inv_J, fx)
        #
        dx = -dx
        print("dx:")
        print(dx)
        d_theta = dx[:num_PQ_nodes + num_PV_nodes]
        assert len(d_theta) == num_PQ_nodes + num_PV_nodes
        # TODO: assert out of possible misunderstanding of numpy grammar. need to be removed
        dV_over_V = dx[num_PQ_nodes + num_PV_nodes:]
        assert len(dV_over_V) == num_PQ_nodes

        for node in network.PQ_nodes:
            # TODO: try += on object's attributes
            node.theta = node.theta + d_theta[node.index]
            node.V = node.V + node.V * dV_over_V[node.index]
            network.update_node(node)
        for node in network.PV_nodes:
            node.theta = node.theta + d_theta[node.index]
            network.update_node(node)

        print("after num_iter: {}".format(num_iter))
        print("V")
        print([node.V for node in network.all_nodes])
        print("theta")
        print([node.theta for node in network.all_nodes])
        print("\n")
# if iteration times > max_num_iter, return None(means power flow calculation divergence)
# return None, max_num_iter
# TODO: REVERT END


print(num_iter)
