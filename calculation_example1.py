from examples_example_input import nodes, lines, transformers
from network import Node, Line, Transformer, Network
from network_utils import check_network, form_inductance_array, print_inductance_array
from newton_raphson import NR

nodes = [Node.from_str(node.strip())
         for node in nodes.strip().split("\n")]
transformers = [Transformer.from_str(transformer.strip())
                for transformer in transformers.strip().split("\n")]
lines = [Line.from_str(line.strip())
         for line in lines.strip().split("\n")]
branches = transformers + lines

node_num = check_network(branches, nodes)
Y = form_inductance_array(branches)
print(Y)
print_inductance_array(Y)
network = Network(nodes, Y)

network, num_iter = NR(network, eps=10 ^ -8, max_num_iter=15)

print(num_iter)
