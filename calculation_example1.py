from network import Node, Line, Transformer, Network
from network_utils import check_branches_and_nodes, print_inductance_array
from newton_raphson import NR

# Step1, define network parameters(lines and transformers)
# Here we use our example parameter set 2
from examples_example_input import parameter_set3 as parameter_set

nodes, lines, transformers, extra_branches = parameter_set

nodes = [Node.from_str(node.strip())
         for node in nodes.strip().split("\n")]

transformers = [Transformer.from_str(transformer.strip())
                for transformer in transformers.strip().split("\n")]
# drop them in case there are some blank lines accidentally
transformers = [transformer for transformer in transformers
                if transformer is not None]

lines = [Line.from_str(line.strip())
         for line in lines.strip().split("\n")]
# drop them in case there are some blank lines accidentally
lines = [line for line in lines
         if line is not None]

branches = transformers + lines  # Both are lists, simply added together

# Step 2, check the parameters and form a network
check_branches_and_nodes(branches, nodes)
network = Network(nodes, branches, extra_branches)
print("Inductance Array after node re-index:\n{}\n".format(network.Y))

# Step 3, calculate using N-R method
# All initial values are passed by `network`,
# you can modify it using network.update_node() or update_node_i_PG() or update_node_i_V() method
eps = 1e-5  # epsilon
max_num_iter = 15
# To see iteration process, set `show_steps` to True, default: False
network, num_iter = NR(network, eps, max_num_iter, show_steps=True)
if network is not None:  # if calculation diverges
    network_loss, (max_V_index, max_V), (min_V_index, min_V) = network.statistics()
    print("Network loss: {}".format(network_loss))
    print("Node {} has maximum voltage of {}.".format(max_V_index, max_V))
    print("Node {} has minimum voltage of {}.".format(min_V_index, min_V))

    network.visualize()
