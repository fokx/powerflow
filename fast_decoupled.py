import numpy as np
import scipy.linalg
from network import Node, Line, Transformer, Network
from network_utils import check_branches_and_nodes, print_inductance_array
from examples_example_input import parameter_set3 as parameter_set

def FDPF(network, eps=1e-8, max_num_iter=15, show_steps=False):
    if show_steps:
        print("\nStarting FDPF calculation\n")

    num_PQ_nodes = network.num_PQ_nodes
    num_PV_nodes = network.num_PV_nodes
    num_Vtheta_ndoes = network.num_Vtheta_nodes

    # B1 denotes B'
    B1 = network.form_B_prime_array()
    p, l, u = scipy.linalg.lu(B1)
    L1 = l
    D1 = np.diag(np.diag(u))

    # B2 denotes B''
    B2 = -np.imag(network.Y)[:num_PQ_nodes, :num_PQ_nodes]
    p, l, u = scipy.linalg.lu(B2)
    L2 = l
    D2 = np.diag(np.diag(u))

    for num_iter in range(max_num_iter):
        for num_iter in range(max_num_iter):
            if show_steps:
                print("***************Starting iteration {}**********************".format(num_iter + 1))

            delta_P_stack = []
            fx_stack = []
            delta_Q_stack = []

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
                delta_P_stack.append(delta_P)
                fx_stack.append(delta_P)

            assert len(delta_P_stack) == num_PQ_nodes + num_PV_nodes
            delta_P_over_V = np.zeros(shape=(num_PQ_nodes + num_PV_nodes, 1))
            for node in network.PQ_nodes + network.PV_nodes:
                delta_P_over_V[node.index, 0] = delta_P_stack[node.index] / node.V

            delta_theta = np.matmul(np.linalg.inv(B1), delta_P_over_V)
            # # Alternatively, no approximation
            # delta_theta_times_V = np.matmul(np.linalg.inv(B1), delta_P_over_V)
            # delta_theta = np.zeros(shape=(num_PQ_nodes + num_PV_nodes, 1))
            # for node in network.PQ_nodes + network.PV_nodes:
            #     delta_theta[node.index,1] = delta_theta_times_V / node.V

            for node in network.PQ_nodes + network.PV_nodes:
                node.theta += float(delta_theta[node.index])
                network.update_node(node)

            for node_i in network.PQ_nodes:
                Qi = 0.0
                for node_j in network.get_adjacent_nodes(node_i):
                    Gij, Bij = network.get_Gij_Bij(node_i, node_j)
                    thetaij = network.get_thetaij(node_i, node_j)
                    Qi += node_i.V * node_j.V * (Gij * np.sin(thetaij) - Bij * np.cos(thetaij))
                delta_Q = node_i.QG - node_i.QL - Qi
                delta_Q_stack.append(delta_Q)
                fx_stack.append(delta_Q)

            assert len(delta_Q_stack) == num_PQ_nodes
            delta_Q_over_V = np.zeros(shape=(num_PQ_nodes, 1))
            for node in network.PQ_nodes:
                delta_Q_over_V[node.index, 0] = delta_Q_stack[node.index] / node.V

            delta_V = np.matmul(np.linalg.inv(B2), delta_Q_over_V)

            for node in network.PQ_nodes:
                node.V += float(delta_V[node.index])
                network.update_node(node)

            abs_delta_P_max = np.max(np.abs(delta_P_stack))
            abs_delta_Q_max = np.max(np.abs(delta_Q_stack))
            abs_fx_max = np.max(np.abs(fx_stack))
            assert abs_fx_max == max([abs_delta_P_max, abs_delta_Q_max])

            if show_steps:
                print("delta P array in interation {}: ".format(num_iter + 1))
                print("{}".format(delta_P_stack))
                print("delta P array in interation {}: ".format(num_iter + 1))
                print("{}".format(delta_P_stack))
                print("max delta P after abs: {}".format(abs_delta_P_max))
                print("max delta Q after abs: {}".format(abs_delta_Q_max))

            if abs_fx_max <= eps:
                if show_steps:
                    print("\nCalculation converges in required number of iterations({} < {}) and epsilon {}.\n".
                          format(num_iter + 1, max_num_iter, eps))
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
    # return None if the power flow calculation diverges
    return None, num_iter + 1