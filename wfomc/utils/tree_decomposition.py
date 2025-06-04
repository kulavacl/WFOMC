import typing
from collections.abc import Set
from typing_extensions import Self
from copy import copy

class TD_binary_evidence:
    #this will be used for so man things
    def __init__(self):
        self.root = None
    
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
        
    def print(self):
        self.root.print()

    def verify_niceness(self):
        q = []
        q.append(self.root)
        while len(q) > 0:
            node = q[0]
            q.pop(0)
            if len(node.children) > 2:
                raise(RuntimeError(f"Node {node.bag} has too many children!"))
            if len(node.children) == 2:
                if len(node.bag - (node.children[0].bag | node.children[1].bag)) > 0 or \
                    len((node.children[0].bag | node.children[1].bag) - node.bag) > 0:
                    raise(RuntimeError(f"Join node {node.bag} is not union of its children!"))
            if len(node.children) == 1:
                if len(node.bag - (node.children[0].bag)) > 0 or \
                    len((node.children[0].bag) - node.bag) != 1:
                    raise(RuntimeError(f"Separator node {node.bag} does not have exactly one element less than its child!"))
            q.extend(node.children)
        print("Niceness verified.")


class TDNode:

    def __init__(self, bag: Set, parent: Self | None):
        self.bag = bag
        self.parent = parent
        self.children = []

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
            len(to_separate) == 1:
            return
        if len(to_separate) == 0:
            return
        self.bag = self.bag - to_separate
        old_children = self.children
        current_node = self
        for elem in to_separate:
            child_bag = copy(current_node.bag)
            child_bag.add(elem)
            current_node.children = [TDNode(child_bag, current_node)]
            current_node = current_node.children[0]
        current_node.children = old_children
        for child in current_node.children:
            child.parent = current_node

    #combine children if there is too many
    def combine_children(self):
        if len(self.children) <= 2:
            return
        child_bag = set()
        for child in self.children[1:]:
            child_bag = child_bag | child.bag
        new_child = TDNode(child_bag, self)
        new_child.children = self.children[1:]
        self.children = [self.children[0], new_child]
        for child in new_child.children:
            child.parent = new_child

    def print(self, offset=0):
        print("\t" * offset, end="")
        print(self.bag)
        for child in self.children:
            child.print(offset + 1)

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
