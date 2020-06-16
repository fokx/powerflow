import numpy as np

'''
including simple functions not limited to Network Class
'''

def print_inductance_array(Y):
    '''
    print non-zero elements

    :param Y: inductance_array
    :return: None
    '''
    print('导纳矩阵非零元素如下：')
    nonzero_i, nonzero_j = np.nonzero(Y)
    print('      行号        列号          实部              虚部      ')
    for i, j in zip(nonzero_i, nonzero_j):
        print('      %2d          %2d        %12.6f     %12.6f '
              % (i, j, np.real(Y[i, j]), np.imag(Y[i, j])))
    print("\n")


def check_branches_and_nodes(branches, nodes):
    '''

    :param branches: a list branches. All bracnches,
        including nodes and transformers
    :param nodes: a list of nodes.
    '''
    # Check 1,
    # assert all branches's i != its j, i.e. not loop for every node
    for branch in branches:
        assert branch.i != branch.j, "Node {} has a loop to itself".format(branch.i)

    # Check 2,
    # assert all nodes do have valid data
    node_index = [node.index for node in nodes]
    node_num = max(node_index)
    for i in range(1, node_num + 1):
        assert i in node_index, "Not all node has valid data"

    # Check 3,
    # assert all branch's index i and j are legal node index
    branch_index_i = [branch.i for branch in branches]
    branch_index_j = [branch.j for branch in branches]
    all_node_index = set(node_index)

    for i in branch_index_i + branch_index_j:
        assert i in all_node_index, "Branch terminal" \
                                    " {} in not a valid node.".format(i)
