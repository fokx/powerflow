import numpy as np


def NR(network,eps,max_num_iter):
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

    delta_P_array = np.zeros(dtype=np.float, shape=(num_PQ_nodes + num_PV_nodes, 1))
    delta_Q_array = np.zeros(dtype=np.float, shape=(num_PQ_nodes, 1))

    # to make them live after iteration in all nodes, declare the below two arrays here.
    fx = np.vstack((delta_P_array, delta_Q_array))
    dx = np.vstack((delta_P_array, delta_Q_array))

    for num_iter in range(max_num_iter):
        for node_i in network.PQ_nodes:
            Pi = 0.0
            Qi = 0.0
            for node_j in network.get_adjacent_nodes(node_i):

                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                assert not (Gij == 0.0 and Bij == 0.0)
                thetaij = network.get_thetaij(node_i, node_j)
                Pi += node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
                Qi += node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
            delta_P = node_i.PG - node_i.PL - Pi
            delta_Q = node_i.QG - node_i.QL - Qi

            # TODO: MAKE SURE REALLY UNDERSTAND INDEX
            # print(node_i.index)
            delta_P_array[node_i.index] = delta_P
            delta_Q_array[node_i.index] = delta_Q

        for node_i in network.PV_nodes:
            Pi = 0.0
            Qi = 0.0
            for node_j in network.get_adjacent_nodes(node_i):
                Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                thetaij = network.get_thetaij(node_i, node_j)
                Pi += node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
            delta_P = node_i.PG - node_i.PL - Pi

            # TODO: MAKE SURE REALLY UNDERSTAND INDEX
            # print(node_i.index)
            delta_P_array[node_i.index] = delta_P

        fx = np.vstack((delta_P_array, delta_Q_array))

        if np.max(fx) <= eps:
            # TODO: REVERT
            # return network, num_iter
            exit(0)
        else:
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
                        J[i, j] = node_i.V ** 2 * Bij + node_i.QG - node_i.QL
                    else:
                        J[i, j] = - node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
                # N
                for j in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
                    node_j = network.get_node(j - (num_PQ_nodes + num_PV_nodes))
                    Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                    thetaij = network.get_thetaij(node_i, node_j)
                    if i - (num_PQ_nodes + num_PV_nodes) == j:
                        J[i, j] = - node_i.V ** 2 * Gij - node_i.PG + node_i.PL
                    else:
                        J[i, j] = - node_i.V * node_j.V * (Gij * np.cos(thetaij) - Bij * np.sin(thetaij))

            for i in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
                # M
                node_i = network.get_node(i - (num_PQ_nodes + num_PV_nodes))
                for j in range(num_PQ_nodes + num_PV_nodes):
                    node_j = network.get_node(j)
                    Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                    thetaij = network.get_thetaij(node_i, node_j)
                    if i == j - num_PQ_nodes + num_PV_nodes:
                        J[i, j] = node_i.V ** 2 * Gij - node_i.PG + node_i.PL
                    else:
                        J[i, j] = node_i.V * node_j.V * (Gij * np.cos(thetaij) + Bij * np.sin(thetaij))
                # L
                for j in range(num_PQ_nodes + num_PV_nodes, num_PQ_nodes + num_PV_nodes + num_PQ_nodes):
                    node_j = network.get_node(j - (num_PQ_nodes + num_PV_nodes))
                    Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                    thetaij = network.get_thetaij(node_i, node_j)
                    if i == j:
                        J[i, j] = node_i.V ** 2 * Bij - node_i.QG + node_i.QL
                    else:
                        J[i, j] = - node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))

            # update network parameters
            inv_J = np.linalg.inv(J)
            dx = np.matmul(inv_J, fx)

            d_theta = dx[:num_PQ_nodes + num_PV_nodes]
            dV_over_V = dx[num_PQ_nodes + num_PV_nodes:]

            # TODO: assert out of possible misunderstanding of numpy grammar. need to be removed
            assert len(d_theta) == num_PQ_nodes + num_PV_nodes
            assert len(dV_over_V) == num_PQ_nodes

            for node in network.PQ_nodes:
                # TODO: try += on object's attributes
                node.theta = node.theta + d_theta[node.index]
                node.V = node.V + dV_over_V[node.index]
                network.update_node(node)
            for node in network.PV_nodes:
                node.theta = node.theta + d_theta[node.index]
                network.update_node(node)

    # if iteration times > max_num_iter, return None(means power flow calculation divergence)
    # return None, max_num_iter
