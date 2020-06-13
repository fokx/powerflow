from enum import Enum

import numpy as np


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
    # TODO: implement Network graph (power flow graph)
    '''
     :param list_of_nodes: list of Node(index, node_type, PG, QG, PL, QL, V0,
                               PGmin, PGmax, QGmin, QGmax, Vmin, Vmax)
    :param Y: inductance matrix if the first row and column of `Y` is stripped;
              calculated by 'form_inductance_array'
    '''

    def __init__(self, list_of_nodes, Y):
        self.all_nodes = []
        self.PQ_nodes = []
        self.PV_nodes = []
        self.Vtheta_nodes = []
        self.Y = Y[1:, 1:]

        # arrangelsit of nodes in correct order: PQ...PV...Vtheta

        list_of_nodes_without_Vtheta = [node for node in list_of_nodes
                                        if node.node_type != NodeType.Vtheta]
        list_of_nodes_without_Vtheta = sorted(list_of_nodes_without_Vtheta, key=lambda node: node.index)
        for new_index, node in enumerate(list_of_nodes_without_Vtheta):
            node.index = new_index
            self.all_nodes.append(node)

        list_of_Vtheta_nodes = [node for node in list_of_nodes
                                if node.node_type == NodeType.Vtheta]

        for new_index, node in enumerate(list_of_Vtheta_nodes):
            node.index = new_index + len(list_of_nodes_without_Vtheta)
            self.all_nodes.append(node)

        for node in self.all_nodes:
            if node.node_type == NodeType.PQ:
                self.PQ_nodes.append(node)
            elif node.node_type == NodeType.PV:
                self.PV_nodes.append(node)
            elif node.node_type == NodeType.Vtheta:
                self.Vtheta_nodes.append(node)

        self.num_all_nodes = len(self.all_nodes)

    def get_node(self, index):
        # BE ATTENTION ABOUT `update_node`'s index
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        return self.all_nodes[index]

    def get_node_using_raw_index(self, index):
        # TODO maybe useless
        assert index >= 0
        assert index < self.num_all_nodes, "{}, {}".format(index, self.num_all_nodes)
        return self.all_nodes[index]

    def get_num_PQ_nodes(self):
        return len(self.PQ_nodes)

    def get_num_PV_nodes(self):
        return len(self.PV_nodes)

    def get_num_Vtheta_nodes(self):
        return len(self.Vtheta_nodes)

    def get_adjacent_nodes(self, node_i: Node):
        '''

        :param node_i: a Node
        :return: list of adjacent nodes
        '''
        adjacent_nodes = []
        list_row_i = list(self.Y[node_i.index, :])
        for j, Yij in enumerate(list_row_i):
            if Yij != 0.0:
                # TODO: MAKE SURE REALLY UNDERSTAND INDEX
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
        self.Y = 1 / self.Z


class Line(Branch):
    def __init__(self, i, j, R, X, half_Y):
        # % Line_para = [始端节点号 末端节点号  电阻   电抗  对地导纳的一半]
        i = int(i)
        j = int(j)
        R = float(R)
        X = float(X)
        half_Y = 1j * float(half_Y)
        super().__init__(i, j, R, X, half_Y, half_Y)

    @classmethod
    def from_str(cls, line):
        # % Line_para = [始端节点号 末端节点号  电阻   电抗  对地导纳的一半]
        line = line.strip().split()
        line = [i.strip() for i in line]
        assert len(line) == 5
        return Line(line[0], line[1], line[2], line[3], line[4])


class Transformer(Branch):
    def __init__(self, i, j, R, X, k):
        # % Trans_para = [始端节点号 末端节点号  电阻   电抗  非标准变比]
        i = int(i)
        j = int(j)
        R = float(R)
        X = float(X)
        k = float(k)
        Z_T = R + 1j * X
        super().__init__(i, j, k * R, k * X, (k - 1) / (k * Z_T), (1 - k) / (k ** 2 * Z_T))

    @classmethod
    def from_str(cls, line):
        # % Trans_para = [始端节点号 末端节点号  电阻   电抗  非标准变比]
        line = line.strip().split()
        line = [i.strip() for i in line]
        assert len(line) == 5
        return Transformer(line[0], line[1], line[2], line[3], line[4])
