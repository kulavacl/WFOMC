import typing
from typing import FrozenSet
from collections.abc import Set
from typing_extensions import Self
from copy import copy
from networkx import Graph
from networkx.algorithms.approximation import treewidth_min_degree
from logzero import logger

class TD_binary_evidence:
    #this will be used for so many things
    def __init__(self, graph: Graph):
        treewidth, graph = treewidth_min_degree(graph)
        logger.info(f"creating nice tree decomposition with treewidth {treewidth}")
        nodes = list(graph.nodes)
        edges = list(graph.edges)
        self.root = TDNode(frozenset(), None)
        self.root.children.append(TDNode(nodes[0], self.root))
        nodes_dict = {self.root.bag : self.root,
                    self.root.children[0].bag : self.root.children[0]}
        current = 0
        while len(edges) > 0:
            edge = edges[current]
            if edge[0] in nodes_dict:
                new_node = TDNode(edge[1], nodes_dict[edge[0]])
                nodes_dict[edge[0]].children.append(new_node)
                nodes_dict[edge[1]] = new_node
                edges.pop(current)
            elif edge[1] in nodes_dict:
                new_node = TDNode(edge[0], nodes_dict[edge[1]])
                nodes_dict[edge[1]].children.append(new_node)
                nodes_dict[edge[0]] = new_node
                edges.pop(current)
            else:
                current += 1
            if current >= len(graph.edges):
                current = 0
        self.be_nice()

    
    #method to make existing good format tree decomposition nice
    #good format tree decomposition not available yet
    def be_nice(self):
        q = []
        q.extend(self.root.children)
        while len(q) > 0:
            node = q[0]
            q.pop(0)
            node.separate_self()
            q.extend(node.children)
        
        q.append(self.root)
        while len(q) > 0:
            node = q[0]
            q.pop(0)
            node.create_child()
            node.combine_children()
            q.extend(node.children)
        
    def print(self, print_configs = False):
        self.root.print(print_configs=print_configs)

    def count_nodes(self):
        return self.root.count_nodes()

    def verify_niceness(self):
        q = []
        q.append(self.root)
        while len(q) > 0:
            node = q[0]
            q.pop(0)
            if len(node.children) > 2:
                self.root.print()
                raise(RuntimeError(f"Node {node.bag} has too many children!"))
            if len(node.children) == 2:
                if len(node.bag - (node.children[0].bag | node.children[1].bag)) > 0 or \
                    len((node.children[0].bag | node.children[1].bag) - node.bag) > 0:
                    self.root.print()
                    raise(RuntimeError(f"Join node {node.bag} is not union of its children!"))
            if len(node.children) == 1:
                if len(node.bag - (node.children[0].bag)) > 0 or \
                    len((node.children[0].bag) - node.bag) != 1:
                    self.root.print()
                    raise(RuntimeError(f"Separator node {node.bag} does not have exactly one element less than its child!"))
            q.extend(node.children)
        #print("Niceness verified.")


class TDNode:

    def __init__(self, bag: FrozenSet, parent: Self | None):
        self.bag = bag
        self.parent = parent
        self.children = []
        self.configs = None

    #create child from elements in self that are not in children
    def create_child(self):
        if len(self.children) == 0:
            return
        child_bag = copy(self.bag)
        for child in self.children:
            child_bag = child_bag - child.bag
        if len(child_bag) > 0:
            self.children.append(TDNode(child_bag, self))

    #separate elements from self
    def separate_self(self):
        to_separate = self.bag - self.parent.bag
        if len(self.parent.children) == 1 and \
            len(to_separate) == 1 and \
            self.bag - to_separate == self.parent.bag:
            return
        if len(to_separate) == 0:
            return
        self.bag = self.bag - to_separate
        if self.parent.bag == frozenset() and \
            len(self.parent.children) == 1:
            elem = next(iter(to_separate))
            self.bag = frozenset([elem])
            to_separate = to_separate - {elem}
        old_children = self.children
        current_node = self
        for elem in to_separate:
            child_bag = current_node.bag | {elem}
            current_node.children = [TDNode(child_bag, current_node)]
            current_node = current_node.children[0]
        current_node.children = old_children
        for child in current_node.children:
            child.parent = current_node

    #combine children if there is too many
    def combine_children(self):
        if len(self.children) <= 2:
            return
        child_bag = frozenset()
        for child in self.children[1:]:
            child_bag = child_bag | child.bag
        new_child = TDNode(child_bag, self)
        new_child.children = self.children[1:]
        self.children = [self.children[0], new_child]
        for child in new_child.children:
            child.parent = new_child

    def print(self, offset=0, print_configs = False):
        print("\t" * offset, end="")
        if not print_configs:
            print(self.bag)
        else:
            print(self.bag, " : ", self.configs)
        for child in self.children:
            child.print(offset + 1, print_configs)

    def count_nodes(self):
        if len(self.children) == 0:
            return 1
        return sum(child.count_nodes() for child in self.children) + 1

if __name__ == "__main__":
    node1 = TDNode(set((3, 4, 5)), None)
    node2 = TDNode(set((2, 4, 5)), node1)
    node3 = TDNode(set((1, 5)), node1)
    node4 = TDNode(set((0, 1, 5)), node3)
    node5 = TDNode(set((7, 9)), node2)
    node1.children = [node2, node3]
    node2.children = [node5]
    node3.children = [node4]
    bin_ev = TD_binary_evidence()
    bin_ev.root = node1
    bin_ev.print()
    bin_ev.be_nice()
    bin_ev.print()
    bin_ev.verify_niceness()
