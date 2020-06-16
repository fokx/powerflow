from enum import Enum
import copy
import numpy as np
import unittest


class ValueRange:
    def __init__(self, min, max):
        assert isinstance(min, float) or isinstance(min, int)
        assert isinstance(max, float) or isinstance(max, int)
        self.min = min
        self.max = max

    def __contains__(self, item):
        assert isinstance(item, float)
        if self.min < item < self.max:
            return True
        else:
            return False

    def __lt__(self, other):
        if other < self.min:
            return True
        else:
            return False

    def __gt__(self, other):
        if other > self.max:
            return True
        else:
            return False


class NodeType(Enum):
    PV = -1
    PQ = 1
    Vtheta = 0


class Node:
    def __init__(self, index, node_type, PG, QG, PL, QL, V0,
                 PGmin, PGmax, QGmin, QGmax, Vmin, Vmax, theta0=0.0):
        self.raw_index = int(index)
        self.index = self.raw_index
        self.node_type = NodeType(int(node_type))
        self.PG = float(PG)
        self.QG = float(QG)
        self.PL = float(PL)
        self.QL = float(QL)
        self.V0 = float(V0)  # init value

        self.S = complex(0)  # only used when this is a Vtheta node
        self.Q = complex(0)  # only used when this is a PV node

        self.theta0 = theta0  # init value
        self.V = self.V0
        self.theta = self.theta0

        self.PG_range = ValueRange(float(PGmin), float(PGmax))
        self.QG_range = ValueRange(float(QGmin), float(QGmax))
        self.V_range = ValueRange(float(Vmin), float(Vmax))
        self.S = self.PG + 1j * self.QG - self.PL - 1j * self.QL

        self.adjacent_nodes = None

    @classmethod
    def from_array(self, arr):
        '''
        :param arr:
            1-d array,
            object dtype can be str or number
            节点号，节点类型，发电机有功，发电机无功，负荷有功，负荷无功，节点电压初值
            发电机有功下限   发电机有功上限    发电机无功下限   发电机无功上限    节点最低电压    节点最高电压
            其中，节点类型为-1表示PV节点，1表示PQ节点，0表示平衡机节点
            其中，对于平衡机和PV节点，节点电压初值即为给定的节点电压，在计算中不变化

            An example of `arr`: 2   1    0.0000    0.0000    0.0000    0.0000    1.0000
                                0.00    0.00    0.00    0.00    0.95    1.07
        :return:
            Node instance
        '''
        assert len(arr) == 13, \
            "Array data unacceptable because of length {} != 13".format(len(arr))

        return Node(arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6],
                    arr[7], arr[8], arr[9], arr[10], arr[11], arr[12])

    @classmethod
    def from_str(cls, line):
        arr = line.strip().split()
        arr = [i.strip() for i in arr]
        assert len(arr) == 13
        return Node.from_array(arr)

    def value_range_check(self, type, value):
        '''

        :param type: "PG", "V", "QG"
        :param value:
        :return:
        '''
        type_to_name = {"PG": "有功功率",
                        "V": "电压",
                        "QG": "无功功率"}
        type_name = type_to_name[type]
        too_high = 'Warning：发电机： 节点%3d 超出%s上限。实际值：%7.3f, 上限值：%7.3f. \r' % (
            self.i, type_name, value, self.PG_range.max)
        too_low = 'Warning：发电机： 节点%3d 低于%s下限。实际值：%7.3f, 下限值：%7.3f. \r' % (
            self.i, type_name, value, self.PG_range.min)
        type_to_ranges = {"PG": self.PG_range,
                          "V": self.V_range,
                          "QG": self.QG_range}
        specific_range = type_to_ranges[type]
        warning = ''
        if value not in specific_range:
            # warning = too_high if value > specific_range.max else too_low
            warning = too_high if value > specific_range else too_low
        print(warning)


class Network():
    # TODO: implement Network graph (power flow graph) in visulize() method
    '''
     :param list_of_nodes: list of Node(index, node_type, PG, QG, PL, QL, V0,
                               PGmin, PGmax, QGmin, QGmax, Vmin, Vmax)
    :param Y: inductance matrix if the first row and column of `Y` is stripped;
              calculated by 'form_inductance_array'
    '''

    def __init__(self, list_of_nodes, branches, extra_branches=None):
        '''

        :param list_of_nodes: [Node]
        :param branches: [Branch]
        :param extra_branches: [(node_index, inductance)]
        '''
        self.all_nodes = []

        # rearrange list of nodes in correct order: PQ...PV...Vtheta

        list_of_PQ_nodes = [node for node in list_of_nodes
                            if node.node_type == NodeType.PQ]
        list_of_PQ_nodes = sorted(list_of_PQ_nodes, key=lambda node: node.index)
        for new_index, node in enumerate(list_of_PQ_nodes):
            # Actually this writes to the instance of Node, thus `list_of_nodes` itself is changed!
            # so, to generate old_to_current_index_map, it is the same to use either `list_of_nodes` or `self.all_nodes`
            # see the assert statement below
            node.index = new_index
            self.all_nodes.append(node)

        list_of_PV_nodes = [node for node in list_of_nodes
                            if node.node_type == NodeType.PV]
        list_of_PV_nodes = sorted(list_of_PV_nodes, key=lambda node: node.index)
        for new_index, node in enumerate(list_of_PV_nodes):
            node.index = new_index + len(list_of_PQ_nodes)
            self.all_nodes.append(node)

        # write to attributes
        list_of_Vtheta_nodes = [node for node in list_of_nodes
                                if node.node_type == NodeType.Vtheta]

        for new_index, node in enumerate(list_of_Vtheta_nodes):
            node.index = new_index + len(list_of_PQ_nodes) + len(list_of_PV_nodes)
            self.all_nodes.append(node)

        self.PQ_nodes = []
        self.PV_nodes = []
        self.Vtheta_nodes = []
        for node in self.all_nodes:
            if node.node_type == NodeType.PQ:
                self.PQ_nodes.append(node)
            elif node.node_type == NodeType.PV:
                self.PV_nodes.append(node)
            elif node.node_type == NodeType.Vtheta:
                self.Vtheta_nodes.append(node)

        self.num_all_nodes = len(self.all_nodes)
        self.num_PQ_nodes = len(self.PQ_nodes)
        self.num_PV_nodes = len(self.PV_nodes)
        self.num_Vtheta_nodes = len(self.Vtheta_nodes)

        # old <-> current index map
        # if Vtheta node is put at last initially, both order and elements are equal
        # otherwise, only elements equal
        for node in self.all_nodes:
            assert node in list_of_nodes
        assert len(self.all_nodes) == len(list_of_nodes)

        list_of_raw_index = [node.raw_index for node in self.all_nodes]
        list_of_current_index = [node.index for node in self.all_nodes]

        self.raw_to_current_index_map = dict(zip(list_of_raw_index, list_of_current_index))
        self.current_to_raw_index_map = dict(zip(list_of_current_index, list_of_raw_index))

        # re-index branches and extra_branches
        self.branches = []
        branch: Branch
        for branch in branches:
            # deep copy is not needed
            # branch_copy = copy.deepcopy(branch)
            # branch_copy.i = self.raw_to_current_index_map[branch.i]
            # branch_copy.j = self.raw_to_current_index_map[branch.j]

            branch.i = self.raw_to_current_index_map[branch.i]
            branch.j = self.raw_to_current_index_map[branch.j]

            # assert branch.i == branch_copy.i
            # assert branch.j == branch_copy.j

            self.branches.append(branch)

        if extra_branches is not None:
            if len(extra_branches) != 0:
                self.extra_branches = []
                for extra_branch in extra_branches:
                    extra_branch = (self.raw_to_current_index_map[extra_branch[0]], extra_branch[1])
                    self.extra_branches.append(extra_branch)
        else:
            self.extra_branches = None
        self.Y = self._form_inductance_array(self.branches, self.extra_branches)

    def get_node(self, index):
        # BE ATTENTION ABOUT `update_node`'s index
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        return self.all_nodes[index]

    def get_node_using_raw_index(self, index):
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        return self.all_nodes[index]

    def get_adjacent_nodes(self, node_i: Node):
        '''

        :param node_i: a Node
        :return: list of adjacent nodes
        '''
        adjacent_nodes = []
        list_row_i = list(self.Y[node_i.index, :])
        for j, Yij in enumerate(list_row_i):
            if Yij != 0.0:
                # assert j<len(self.all_nodes), "{} > {} of {}".format(j, len(self.all_nodes), Yij)
                adjacent_nodes.append(self.get_node(j))

        return adjacent_nodes

    def get_Yij(self, node_i: Node, node_j: Node):
        return self.Y[node_i.index, node_j.index]

    def get_Gij_Bij(self, node_i: Node, node_j: Node):
        Yij = self.get_Yij(node_i, node_j)
        Gij = np.real(Yij)
        Bij = np.imag(Yij)
        return Gij, Bij

    def get_thetaij(self, node_i: Node, node_j: Node):
        return node_i.theta - node_j.theta

    def update_node(self, node: Node):
        assert node.index >= 0
        assert node.index < self.num_all_nodes, "{}, {}".format(node.index, self.num_all_nodes)
        self.all_nodes[node.index] = node

    def update_node_i_PG(self, index, PG):
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        node_i: Node = self.get_node(index)
        node_i.PG = PG
        self.update_node(node_i)

    def update_node_i_V(self, index, V):
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        node_i: Node = self.get_node(index)
        node_i.V = V
        self.update_node(node_i)

    def visualize(self):
        pass

    def statistics(self):

        network_loss = 0.
        node: Node
        for node in self.all_nodes:
            if node.node_type == NodeType.Vtheta:
                network_loss += np.real(node.S)
            else:
                network_loss += node.PG - node.PL

        all_V = [node.V for node in self.all_nodes]
        max_V_index = np.argmax(all_V)
        max_V = all_V[max_V_index]
        min_V_index = np.argmin(all_V)
        min_V = all_V[min_V_index]
        return network_loss, (max_V_index, max_V), (min_V_index, min_V)

    def _form_inductance_array(self, branches, extra_branches):
        '''
        :param extra_branches: list of tuple: (node_index, inductance_to_the_ground)
        :param branches: a list branches. All bracnches,
            including lines and transformers
        :return: inductance matrix
        '''
        branch_index_i = [branch.i for branch in branches]
        branch_index_j = [branch.j for branch in branches]
        all_node_index = set(branch_index_i + branch_index_j)
        node_num = len(all_node_index)

        Y = np.zeros(dtype=np.complex, shape=(node_num, node_num))

        for branch in branches:
            # Y1 is connected to node i
            # Y2 is connected to node j
            # self.Z = R + 1j * X
            # self.Y = 1 / self.Z
            i = branch.i
            j = branch.j

            Y[i, i] += branch.Y1 + branch.Y
            Y[j, j] += branch.Y2 + branch.Y

            Y[i, j] += -branch.Y
            Y[j, i] += -branch.Y

        if extra_branches is not None:
            if len(extra_branches) != 0:
                for extra_branch in extra_branches:
                    Y[extra_branch[0], extra_branch[0]] += extra_branch[1]

        return Y

    def form_B_prime_array(self, ignore_transformer_k=True):
        '''
        used to generate B' used in Fast Decoupled Power Flow
        :param ignore_transformer_k: use X_T instead of k_T * X_T
        :param branches: a list branches. All bracnches,
            including lines and transformers
        :return: B' matrix
        '''
        B_prime = np.zeros(dtype=np.float, shape=(self.num_all_nodes, self.num_all_nodes))

        if ignore_transformer_k:
            branch: Branch
            for branch in self.branches:
                i = branch.i
                j = branch.j
                if isinstance(branch, Line):
                    B_prime[i, j] -= 1 / branch.X
                    B_prime[j, i] -= 1 / branch.X
                    B_prime[i, i] += 1 / branch.X
                    B_prime[j, j] += 1 / branch.X
                elif isinstance(branch, Transformer):
                    B_prime[i, j] -= 1 / branch.X_T
                    B_prime[j, i] -= 1 / branch.X_T
                    B_prime[i, i] += 1 / branch.X_T
                    B_prime[j, j] += 1 / branch.X_T
        else:
            branch: Branch
            for branch in self.branches:
                i = branch.i
                j = branch.j
                B_prime[i, j] -= 1 / branch.X
                B_prime[j, i] -= 1 / branch.X
                B_prime[i, i] += 1 / branch.X
                B_prime[j, j] += 1 / branch.X

        return B_prime[:self.num_PQ_nodes + self.num_PV_nodes,
               :self.num_PQ_nodes + self.num_PV_nodes]


class Branch:
    # sub-class: Line, Transformer
    def __init__(self, i, j, R, X, Y1, Y2):
        # Y1 is connected to node i
        # Y2 is connected to node j
        self.i = i
        self.j = j
        self.R = R
        self.X = X
        self.Y1 = Y1
        self.Y2 = Y2
        self.G1 = np.real(Y1)
        self.B1 = np.imag(Y1)
        self.G2 = np.real(Y2)
        self.B2 = np.imag(Y2)
        self.Z = R + 1j * X
        if self.Z != 0:
            self.Y = 1 / self.Z
        else:
            self.Y = np.inf


class Line(Branch):
    def __init__(self, i, j, R, X, half_Y):
        # % Line_para = [始端节点号 末端节点号  电阻   电抗  对地导纳的一半]
        self.i = int(i)
        self.j = int(j)
        self.R = float(R)
        self.X = float(X)
        self.half_Y = 1j * float(half_Y)
        super().__init__(self.i, self.j, self.R, self.X, self.half_Y, self.half_Y)

    @classmethod
    def from_str(cls, line):
        # % Line_para = [始端节点号 末端节点号  电阻   电抗  对地导纳的一半]
        line = line.strip().split()
        line = [i.strip() for i in line]
        assert len(line) == 5
        return Line(line[0], line[1], line[2], line[3], line[4])


class Transformer(Branch):
    def __init__(self, i, j, R_T, X_T, k_T):
        # % Trans_para = [始端节点号 末端节点号  电阻   电抗  非标准变比]
        self.i = int(i)
        self.j = int(j)
        self.R_T = float(R_T)
        self.X_T = float(X_T)
        self.k_T = float(k_T)
        self.Z_T = self.R_T + 1j * self.X_T
        super().__init__(self.i, self.j, self.k_T * self.R_T, self.k_T * self.X_T,
                         (self.k_T - 1) / (self.k_T * self.Z_T),
                         (1 - self.k_T) / (self.k_T ** 2 * self.Z_T))

    @classmethod
    def from_str(cls, line):
        if len(line) != 0:
            # % Trans_para = [始端节点号 末端节点号  电阻   电抗  非标准变比]
            line = line.strip().split()
            line = [i.strip() for i in line]
            assert len(line) == 5, "{} with length {}".format(len(line), line)
            return Transformer(line[0], line[1], line[2], line[3], line[4])
        else:
            # TODO: HANDLE None better
            return None
