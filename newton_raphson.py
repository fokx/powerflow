import numpy as np
from network import Node, Line, Transformer, Network

def NR(network, eps=1e-8, max_num_iter=15, show_steps=False):
    '''
    won't consider Y as a sparse matrix
    :param max_num_iter: if iteration times > max_num_iter, return None
    :param eps: iteration termination criteria, epsilon

    :return: (1) node voltage(amplitude and phase angle) array if iteration converges,
             else None
             (2) iteration times taken
    '''
    global num_iter
    num_PQ_nodes = network.get_num_PQ_nodes()
    num_PV_nodes = network.get_num_PV_nodes()

    for num_iter in range(max_num_iter):

        fx_stack = []
        Pix_stack = []
        Qix_stack = []
        for node_i in network.PQ_nodes + network.PV_nodes:
            Pi = 0.0
            adjacent_nodes = network.get_adjacent_nodes(node_i)
            for node_j in adjacent_nodes:
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                Pi += node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
            delta_P = node_i.PG - node_i.PL - Pi
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
            if show_steps:
                print("Calculation converges in required number of iteration({} < {}).".
                      format(num_iter, max_num_iter))

            # calculate S of Vtheta nodes and Q of PV nodes
            Vtheta_node: Node
            for Vtheta_node in network.Vtheta_nodes:
                P = 0.
                Q = 0.
                for node in network.all_nodes:
                    Gij, Bij = network.get_Gij_Bij(Vtheta_node, node)
                    thetaij = network.get_thetaij(Vtheta_node, node)
                    Q += Vtheta_node.V * node.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
                    P += Vtheta_node.V * node.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
                Vtheta_node.S = P + 1j * Q
                network.update_node(Vtheta_node)

            PV_node: Node
            for PV_node in network.PV_nodes:
                P = 0.
                Q = 0.
                for node in network.all_nodes:
                    Gij, Bij = network.get_Gij_Bij(PV_node, node)
                    thetaij = network.get_thetaij(PV_node, node)
                    Q += PV_node.V * node.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
                PV_node.Q = Q
                network.update_node(PV_node)

            # return the updated network and number of iterations
            return network, num_iter + 1  # beacause num_iter starts from 0, add 1 to `num_iter`

        else:
            if show_steps:
                print("***************Starting num_iter: {}**********************".format(num_iter))
                print("delta P and delta Q array:")
                print("{}".format(fx_stack))
                print("max time after abs: {}".format(abs_fx_max))
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
                        J[i, j] = node_i.V ** 2 * Bij - Qix_stack[i - (num_PQ_nodes + num_PV_nodes)]
                    else:
                        J[i, j] = - node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
            if show_steps:
                print("-" * 10)
                print("J")
                print(J)

            # update network parameters
            inv_J = np.linalg.inv(J)
            fx = np.reshape(fx_stack, (len(fx_stack), 1))
            dx = -np.matmul(inv_J, fx)
            d_theta = dx[:num_PQ_nodes + num_PV_nodes]
            dV_over_V = dx[num_PQ_nodes + num_PV_nodes:]
            if show_steps:
                print("-" * 10)
                print("dx")
                print(dx)
            for node in network.PQ_nodes:
                node.theta += d_theta[node.index]
                node.V += node.V * dV_over_V[node.index]
                network.update_node(node)
            for node in network.PV_nodes:
                node.theta += d_theta[node.index]
                network.update_node(node)
            if show_steps:
                print("-" * 10)
                print("after num_iter: {}".format(num_iter))
                print("V")
                print([node.V for node in network.all_nodes])
                print("theta(degrees)")
                print([node.theta * 180 / np.pi for node in network.all_nodes])
                print("\n")

    # return None if the power flow calculation diverges
    return None, num_iter + 1
