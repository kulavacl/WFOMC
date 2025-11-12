from typing import Callable
import typing
from functools import reduce

from wfomc.cell_graph import build_cell_graphs, Cell, CellGraph
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula, a, b, X, AtomicFormula
import math
from copy import copy
import numpy as np
from networkx import Graph
from networkx.algorithms.approximation import treewidth_min_degree
from wfomc.utils.tree_decomposition import TD_binary_evidence, TDNode
import itertools
from logzero import logger
from collections import defaultdict
from functools import lru_cache, reduce
import numpy as np
from operator import mul
from symengine import expand

from time import time

np.random.seed(42)
#defaultdicts can be replaces by dict.get to increase performance

def td_wfomc(context: WFOMCContext) -> RingElement:
    formula = context.formula
    domain = context.domain
    get_weight = context.get_weight
    leq_pred = context.leq_pred
    predecessor_pred = context.predecessor_pred
    successor_pred = context.successor_pred
    first_pred = context.first_pred
    last_pred = context.last_pred
    evidence = context.evidence
    closed_world = context.closed_world
    res = Rational(0, 1)
    domain_size = len(domain)

    #Handle base case evidence for closed world evidence
    if closed_world is None:
        binary_evidence = defaultdict(lambda: None)
    else:
        if any(ev.pred == pred and ev.positive == False for ev in evidence for pred in closed_world):
            raise ValueError("Evidence on closed world predicates can only be positive!")
        BASE_BINARY_EVIDENCE = frozenset(AtomicFormula(pred=pred, args=(arg1, arg2), positive=False) \
            for pred in closed_world for arg1, arg2 in ((a, b), (b, a)))
        binary_evidence = defaultdict(
            lambda : BASE_BINARY_EVIDENCE)

    if closed_world is None:
        unary_evidence = {}
    else:
        unary_evidence = {}
    evidence_copy = copy(evidence)

    #handle unary evidence
    while len(evidence_copy) > 0:
        ev = evidence_copy.pop()
        if ev.pred.arity == 2:
            if (ev.args[0] == ev.args[1]):
                if ev.args[0] not in unary_evidence:
                    unary_evidence[ev.args[0]] = []
                unary_evidence[ev.args[0]].append(ev.substitute({ev.args[0]: X}))
                continue
            if (ev.args[0], ev.args[1]) not in binary_evidence:
                binary_evidence[(ev.args[0], ev.args[1])] = []
                binary_evidence[(ev.args[1], ev.args[0])] = []
            binary_evidence[(ev.args[1], ev.args[0])].append(ev.substitute({ev.args[0]: b, ev.args[1]: a}))
            binary_evidence[(ev.args[0], ev.args[1])].append(ev.substitute({ev.args[0]: a, ev.args[1]: b}))
        elif ev.pred.arity == 1:
            if ev.args[0] not in unary_evidence:
                unary_evidence[ev.args[0]] = []
            unary_evidence[ev.args[0]].append(ev.substitute({ev.args[0]: X}))
    
    #handle binary evidence
    for key, value in binary_evidence.items():
        if closed_world is not None:
            for pred in closed_world:
                if all(ev.pred != pred for ev in value):
                    binary_evidence[key].extend( \
                        [AtomicFormula(pred=pred, args=(a, b), positive=False), \
                         AtomicFormula(pred=pred, args=(b, a), positive=False)])
                elif AtomicFormula(pred=pred, args=(a, b), positive=True) not in value:
                    binary_evidence[key].append(AtomicFormula(pred=pred, args=(a, b), positive=False))
                elif AtomicFormula(pred=pred, args=(b, a), positive=True) not in value:
                    binary_evidence[key].append(AtomicFormula(pred=pred, args=(b, a), positive=False))
        binary_evidence[key] = frozenset(value)
    
    #add reflexive atoms if closed world evidence
    if closed_world is not None:
        for elem, ev in unary_evidence.items():
            unary_evidence[elem].extend(AtomicFormula(pred=pred, args=(X, X), positive=False) for pred in closed_world if \
                                        AtomicFormula(pred=pred, args=(X, X), positive=True) not in unary_evidence[elem])

    logger.info(f"binary evidence found: {evidence}")

    #create gaifman graph
    gaifman_graph = Graph()
    for element in domain:
        gaifman_graph.add_node(element)

    for elem1, elem2 in binary_evidence.keys():
        gaifman_graph.add_edge(elem1, elem2)

    #start WFOMC procedure
    for cell_graph, weight in build_cell_graphs(
        formula, get_weight,
        leq_pred=leq_pred,
        predecessor_pred=predecessor_pred,
        successor_pred=successor_pred,
        first_pred=first_pred,
        last_pred=last_pred,
        domain_size=domain_size
    ):
        #create tree decomposition
        td_bin_ev = TD_binary_evidence(gaifman_graph)
        #td_bin_ev.print()
        td_bin_ev.verify_niceness()
        #cell_graph.show()
        
        cells = cell_graph.get_cells()

        n_cells = len(cells)
        domain_size = len(domain)
        stack = [td_bin_ev.root]

        #prepare unary evidence
        if unary_evidence != {}:
            if closed_world is None:
                BASE_POSSIBLE_CELLS = tuple(i for i in range(n_cells))
            else:
                BASE_POSSIBLE_CELLS = tuple(i for i, cell in enumerate(cells) \
                                             if all(not cell.is_positive(pred) \
                                             for pred in closed_world))
            elem2cells = defaultdict(lambda : BASE_POSSIBLE_CELLS)
            for elem, ev in unary_evidence.items():
                elem2cells[elem] = tuple(i for i, cell in enumerate(cells) if \
                                         all(cell.is_positive(_ev.pred) == _ev.positive for _ev in ev))
        elif closed_world is not None:
            BASE_POSSIBLE_CELLS = tuple(i for i, cell in enumerate(cells) \
                                         if all(not cell.is_positive(pred) \
                                         for pred in closed_world))
            elem2cells = defaultdict(lambda : BASE_POSSIBLE_CELLS)
        else:
            elem2cells = None

        #print(elem2cells)

        #start recursion
        while len(stack) > 0:
            to_process = stack[0]
            if any(child.configs is None for child in to_process.children):
                    stack = [child for child in to_process.children if child.configs is None] + stack
            else:
                if len(to_process.children) == 0:
                    #process_leaf(to_process, n_cells, elem2cells)
                    process_leaf_with_satisfiability_check(to_process, cells, n_cells, cell_graph, elem2cells, binary_evidence)
                elif len(to_process.children) == 1:
                    process_separator(to_process, cell_graph, cells, binary_evidence)
                else:
                    process_join_with_satisfiability_check(to_process, cells, cell_graph, binary_evidence)
                    #process_join(to_process, cells, cell_graph, binary_evidence)
                #print(to_process.bag, to_process.configs)
                stack.pop(0)
        #wfomc = sum in root (not yet, now it is not reuired that root.bag = {})
        res += sum(count for count in td_bin_ev.root.configs[()].values()) * weight
    return res

def process_leaf(node: TDNode, n_cells, elem2cells = None):
    node.configs = dict()
    if elem2cells is None:
        for nu in itertools.product(list(range(n_cells)), repeat=len(node.bag)):
            node.configs[nu] = {tuple(0 for i in range(n_cells)) : Rational(1, 1)}
    else:
        ordered_bag = sorted(node.bag)
        for nu in itertools.product(*[elem2cells[elem] for elem in ordered_bag]):
            node.configs[nu] = {tuple(0 for i in range(n_cells)) : Rational(1, 1)}
    if node.configs == {}:
        raise RuntimeError(f"Formula is unsatisfiable due to unary evidence on bag {node.bag}")
    
def process_leaf_with_satisfiability_check(node: TDNode, cells, n_cells, cell_graph, elem2cells = None, binary_evidence=None):
    #handle leaf nodes
    node.configs = dict()
    ordered_bag = sorted(node.bag)
    #print("leaf", sum(len(child.configs[key]) for child in node.children for key in child.configs.keys()))
    if elem2cells is None:
        for nu in itertools.product(list(range(n_cells)), repeat=len(node.bag)):
            if any(cell_graph.get_two_table_weight((cells[nu1], cells[nu2]), binary_evidence[elem1, elem2]) == 0 \
                   if elem1 != elem2 else False \
                   for elem1, nu1 in zip(ordered_bag, nu) for elem2, nu2 in zip(ordered_bag, nu)):
                continue
            node.configs[nu] = {tuple(0 for i in range(n_cells)) : Rational(1, 1)}
    else:
        #ordered_bag = sorted(node.bag)
        for nu in itertools.product(*[elem2cells[elem] for elem in ordered_bag]):
            if any(cell_graph.get_two_table_weight((cells[nu1], cells[nu2]), binary_evidence[elem1, elem2]) == 0 \
                   if elem1 != elem2 else False \
                   for elem1, nu1 in zip(ordered_bag, nu) for elem2, nu2 in zip(ordered_bag, nu)):
                continue
            node.configs[nu] = {tuple(0 for i in range(n_cells)) : Rational(1, 1)}
    if node.configs == {}:
        raise RuntimeError(f"Formula is unsatisfiable due to unary evidence on bag {node.bag}")
    
def get_combined_nu(left_nu, left_indices, right_nu, right_indices, nu_length):
    new_nu = [-1] * nu_length
    for index, val in zip(left_indices, left_nu):
        new_nu[index] = val
    for index, val in zip(right_indices, right_nu):
        if new_nu[index] == -1:
            new_nu[index] = val
        elif new_nu[index] != val:
            return None
    return tuple(new_nu)


def process_join(node: TDNode, cells, cell_graph, elem2evidence):
    left_bag = node.children[0].bag
    right_bag = node.children[1].bag
    left_configs = node.children[0].configs
    right_configs = node.children[1].configs
    ordered_bag = sorted(node.bag)
    left_indices = tuple(ordered_bag.index(x) for x in sorted(left_bag))
    right_indices = tuple(ordered_bag.index(x) for x in sorted(right_bag))
    node.configs = dict()
    for left_nu, left_WFOMCConfig in left_configs.items():
        for right_nu, right_WFOMCConfig in right_configs.items():
            new_nu = get_combined_nu(left_nu, left_indices, right_nu, right_indices, len(node.bag))
            if new_nu is None:
                continue
            for left_config, left_weight in left_WFOMCConfig.items():
                for right_config, right_weight in right_WFOMCConfig.items():
                    new_config = tuple(i + j for i, j in zip(left_config, right_config))
                    w_12 = 1
                    w_21 = 1
                    for elem, cell_idx in zip(ordered_bag, new_nu):
                        if elem not in left_bag:
                            # combine with left elements
                            for count, cell_k in zip(left_config, cells):
                                w_12 *= cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]) ** count
                        elif elem not in right_bag:
                            # combine with right elements
                            for count, cell_k in zip(right_config, cells):
                                w_21 *= cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]) ** count
                    w_ij = 1
                    for i, cell1 in enumerate(cells):
                        for j, cell2 in enumerate(cells):
                                w_ij *= cell_graph.get_two_table_weight((cell1, cell2), elem2evidence[None])**(left_config[i] * right_config[j])
                    #new_weight = expand(left_weight * right_weight * w_12 * w_21 * w_ij)
                    new_weight = left_weight * right_weight * w_12 * w_21 * w_ij
                    if new_weight != 0:
                        if new_nu not in node.configs:
                            node.configs[new_nu] = {new_config : new_weight}
                        else:
                            node.configs[new_nu][new_config] = node.configs[new_nu].get(new_config, 0) + new_weight

def process_join_with_satisfiability_check(node: TDNode, cells, cell_graph, elem2evidence):
    #handle join nodes
    left_bag = node.children[0].bag
    right_bag = node.children[1].bag
    left_configs = node.children[0].configs
    right_configs = node.children[1].configs
    ordered_bag = sorted(node.bag)
    left_indices = tuple(ordered_bag.index(x) for x in sorted(left_bag))
    right_indices = tuple(ordered_bag.index(x) for x in sorted(right_bag))
    node.configs = dict()
    #print("join", sum(len(child.configs[key]) for child in node.children for key in child.configs.keys()))
    
    @lru_cache(maxsize = None)
    def cached_power(a, b):
        return a ** b
    
    for left_nu, left_WFOMCConfig in left_configs.items():
        for right_nu, right_WFOMCConfig in right_configs.items():
            new_nu = get_combined_nu(left_nu, left_indices, right_nu, right_indices, len(node.bag))
            if new_nu is None:
                continue
            if any(cell_graph.get_two_table_weight((cells[nu1], cells[nu2]), elem2evidence[elem1, elem2]) == 0 \
                   if elem1 != elem2 else False \
                   for elem1, nu1 in zip(ordered_bag, new_nu) for elem2, nu2 in zip(ordered_bag, new_nu)):
                continue

            for left_config, left_weight in left_WFOMCConfig.items():
                for right_config, right_weight in right_WFOMCConfig.items():
                    if any(cell_graph.get_two_table_weight((cells[i], cells[j]), elem2evidence[None]) == 0 \
                          if left_config[i] != 0 and \
                             right_config[j] != 0 else False \
                          for i in range(len(cells)) for j in range(len(cells))):
                        continue
                    new_config = tuple(i + j for i, j in zip(left_config, right_config))
                    w_12 = 1
                    w_21 = 1
                    for elem, cell_idx in zip(ordered_bag, new_nu):
                        if elem not in left_bag:
                            # combine with left elements
                            for count, cell_k in zip(left_config, cells):
                                #w_12 *= cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]) ** count
                                w_12 *= cached_power(cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]), count)
                        elif elem not in right_bag:
                            # combine with right elements
                            for count, cell_k in zip(right_config, cells):
                                #w_21 *= cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]) ** count
                                w_21 *= cached_power(cell_graph.get_two_table_weight((cells[cell_idx], cell_k), elem2evidence[None]), count)
                    if w_12 == 0 or w_21 == 0:
                        continue
                    

                    w_ij = 1
                    # combine cells without evidence
                    for i, cell1 in enumerate(cells):
                        for j, cell2 in enumerate(cells):
                                two_table_weight = cell_graph.get_two_table_weight((cell1, cell2), elem2evidence[None])
                                #if two_table_weight != 1:
                                w_ij *= cached_power(two_table_weight, left_config[i] * right_config[j])
                                #w_ij *= cell_graph.get_two_table_weight((cell1, cell2), elem2evidence[None])**(left_config[i] * right_config[j])
                    
                    
                    new_weight = left_weight * right_weight * w_12 * w_21 * w_ij
                    #if new_weight != 0:
                    if new_nu not in node.configs:
                        node.configs[new_nu] = {new_config : new_weight}
                    else:
                        node.configs[new_nu][new_config] = node.configs[new_nu].get(new_config, 0) + new_weight

def process_separator(node: TDNode, cell_graph: CellGraph, cells: list[Cell], elem2evidence: dict):
    #handle separator nodes
    node.configs = dict()
    separated_elem = next(iter(node.children[0].bag - node.bag))
    ordered_bag = sorted(node.children[0].bag)
    separated_idx = ordered_bag.index(separated_elem)
    for nu, WFOMCConfigs in node.children[0].configs.items():
        separated_cell = cells[nu[separated_idx]]
        new_nu = (*nu[:separated_idx], *nu[separated_idx + 1:])
        for config, weight in WFOMCConfigs.items():
            new_weight = cell_graph.get_cell_weight(separated_cell) * weight #cell graph not as argument yet
            for elem, i in zip(ordered_bag, nu):
                if elem == separated_elem:
                    continue
                new_weight *= cell_graph.get_two_table_weight((separated_cell, cells[i]), elem2evidence[separated_elem, elem])
            #new_weight = expand(new_weight)
            if new_weight != 0:
                new_config = tuple(num + 1 if i == nu[separated_idx] else num for i, num in enumerate(config))
                if new_nu not in node.configs:
                    node.configs[new_nu] = {new_config : new_weight}
                else:
                    node.configs[new_nu][new_config] = node.configs[new_nu].get(new_config, 0) + new_weight
    if node.configs == {}:
        raise RuntimeError(f"Formula is unsatisfiable due to binary evidence on bag {node.children[0].bag}")

