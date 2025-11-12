from typing import Callable
from functools import reduce

from wfomc.cell_graph import build_cell_graphs
from wfomc.context.wfomc_context import WFOMCContext
from wfomc.utils import RingElement, Rational
from wfomc.fol.syntax import Const, Pred, QFFormula
import math


def incremental_wfomc(context: WFOMCContext, t) -> RingElement:
    formula = context.formula
    domain = context.domain
    get_weight = context.get_weight
    leq_pred = context.leq_pred
    predecessor_pred = context.predecessor_pred
    successor_pred = context.successor_pred
    first_pred = context.first_pred
    last_pred = context.last_pred
    if first_pred is not None or last_pred is not None:
        if successor_pred is None:
            raise RuntimeError("FIRST and LAST are only supported with SUC at the moment!")
    if predecessor_pred is not None and successor_pred is not None:
        raise RuntimeError("PRED is only supported without SUC at the moment!")
    res = Rational(0, 1)
    domain_size = len(domain)
    for cell_graph, weight in build_cell_graphs(
        formula, get_weight,
        leq_pred=leq_pred,
        predecessor_pred=predecessor_pred,
        successor_pred=successor_pred,
        first_pred=first_pred,
        last_pred=last_pred,
        domain_size=domain_size
    ):
        #cell_graph.show()
        cells = cell_graph.get_cells()
        n_cells = len(cells)
        domain_size = len(domain)

        if first_pred is not None and last_pred is not None and domain_size > 1:
            first_cells = []
            body_cells = []
            last_cells = []
            for i, cell in enumerate(cells):
                if cell.is_positive(first_pred):
                    first_cells.append((i, cell))
                elif cell.is_positive(last_pred):
                    last_cells.append((i, cell))
                else:
                    body_cells.append((i, cell))
        elif first_pred is not None and domain_size > 1:
            first_cells = []
            body_cells = []
            for i, cell in enumerate(cells):
                if cell.is_positive(first_pred):
                    first_cells.append((i, cell))
                else:
                    body_cells.append((i, cell))
        elif last_pred is not None and domain_size > 1:
            last_cells = []
            body_cells = []
            for i, cell in enumerate(cells):
                if cell.is_positive(last_pred):
                    last_cells.append((i, cell))
                else:
                    body_cells.append((i, cell))
        elif first_pred is not None:
            first_cells = [(i, cell) for i, cell in enumerate(cells)]

        if successor_pred is not None and first_pred is not None:
            table = dict(
                (
                    (
                        tuple(int(k == i) for k in range(n_cells)),
                        tuple(int(i * n_cells + i == k) for k in range(n_cells**2))
                    ),
                    cell_graph.get_cell_weight(cell),
                )
                for i, cell in first_cells
            )
        elif successor_pred is not None and last_pred is not None:
            table = dict(
                (
                    (
                        tuple(int(k == i) for k in range(n_cells)),
                        tuple(int(i * n_cells + i == k) for k in range(n_cells**2))
                    ),
                    cell_graph.get_cell_weight(cell),
                )
                for i, cell in enumerate(body_cells)
            )
        elif successor_pred is not None:
            table = dict(
                (
                    (
                        tuple(int(k == i) for k in range(n_cells)),
                        tuple(int(i * n_cells + i == k) for k in range(n_cells**2))
                    ),
                    cell_graph.get_cell_weight(cell),
                )
                for i, cell in enumerate(cells)
            )
        elif predecessor_pred is not None:
            table = dict(
                (
                    (
                        tuple(int(k == i) for k in range(n_cells)),
                        cell
                    ),
                    cell_graph.get_cell_weight(cell),
                )
                for i, cell in enumerate(cells)
            )
        else:
            table = dict(
                (
                    tuple(int(k == i) for k in range(n_cells)),
                    cell_graph.get_cell_weight(cell),
                )
                for i, cell in enumerate(cells)
            )

        if successor_pred is None:
            for _ in range(domain_size - 1):
                old_table = table
                table = dict()
                for j, cell in enumerate(cells):
                    w = cell_graph.get_cell_weight(cell)
                    for key, w_old in old_table.items():
                        if predecessor_pred is None:
                            ivec = key
                            w_new = w_old * w * reduce(
                                lambda x, y: x * y,
                                (
                                    cell_graph.get_two_table_weight((cell, cells[k]))
                                    ** int(ivec[k]) for k in range(n_cells)
                                ),
                                Rational(1, 1)
                            )
                        else:
                            ivec, last_cell = key
                            w_new = w_old * w
                            for k, other_cell in enumerate(cells):
                                if other_cell == last_cell:
                                    w_new = (
                                        w_new * cell_graph.get_two_table_with_pred_weight((cell, other_cell))
                                        * cell_graph.get_two_table_weight((cell, other_cell)) ** max(ivec[k] - 1, 0)
                                    )
                                else:
                                    w_new = w_new * cell_graph.get_two_table_weight((cell, other_cell)) ** ivec[k]
                        ivec = list(ivec)
                        ivec[j] += 1
                        ivec = tuple(ivec)
                        if predecessor_pred is None:
                            w_new = w_new + table.get(ivec, Rational(0, 1))
                            table[ivec] = w_new
                        else:
                            w_new = w_new + table.get((ivec, cell), Rational(0, 1))
                            table[(tuple(ivec), cell)] = w_new
            res = res + weight * sum(table.values())
            continue
        def merge(h_old, sigma, rho, cell, m):
            out = []
            for a, cell_a in enumerate(cells):
                for b, cell_b in enumerate(cells):
                    if rho[a * n_cells + b] <= 0:
                        continue
                    for c, cell_c in enumerate(cells):
                        #merge to as if connecting a -'-,-> b to b-'-,-> c
                        #merge2
                        #----------------------------------------------------------------
                        if rho[b * n_cells + c] > (0 if a != b or b != c else 1):
                            eta = Rational(rho[a * n_cells + b] * rho[b * n_cells + c], 1) \
                                if a != b or b != c else \
                                Rational(rho[a * n_cells + b] * (rho[a * n_cells + b] - 1), 1)
                            W_cell = cell_graph.get_cell_weight(cell)
                            lambda_ = cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type = 2) * \
                                cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type = 3) * \
                                (cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type = 1) ** (sigma[b] - 2)) * \
                                math.prod(1 if (e == b) else \
                                (cell_graph.get_two_table_with_suc_weight((cell_s, cell), suc_type = 1) ** sigma[e])\
                                for e, cell_s in enumerate(cells))
                            new_rho = list(rho)
                            new_rho[a * n_cells + b] -= 1
                            new_rho[b * n_cells + c] -= 1
                            new_rho[a * n_cells + c] += 1
                            new_sigma = list(sigma)
                            new_sigma[m] += 1
                            new_key = (tuple(new_sigma), tuple(new_rho))
                            if h_old * eta * W_cell * lambda_ != 0:
                                out.append((new_key, h_old * eta * W_cell * lambda_))
                        #----------------------------------------------------------------
                        if b == c:
                            continue
                        for d, cell_d in enumerate(cells):
                            if rho[c * n_cells + d] <= (0 if a != c or b != d else 1):
                                continue
                            #merge1
                            eta = Rational(rho[a * n_cells + b] * rho[c * n_cells + d], 1) \
                                if a != c or b != d else \
                                Rational(rho[a * n_cells + b] * (rho[a * n_cells + b] - 1), 1)
                            W_cell = cell_graph.get_cell_weight(cell)
                            lambda_ = cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type = 2) * \
                                cell_graph.get_two_table_with_suc_weight((cell_c, cell), suc_type = 3) * \
                                (cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type = 1) ** (sigma[b] - 1)) * \
                                (cell_graph.get_two_table_with_suc_weight((cell_c, cell), suc_type = 1) ** (sigma[c] - 1)) * \
                                math.prod(1 if (e == b or e == c) else \
                                (cell_graph.get_two_table_with_suc_weight((cell_s, cell), suc_type = 1) ** sigma[e])\
                                for e, cell_s in enumerate(cells))
                            new_rho = list(rho)
                            new_rho[a * n_cells + b] -= 1
                            new_rho[c * n_cells + d] -= 1
                            new_rho[a * n_cells + d] += 1
                            new_sigma = list(sigma)
                            new_sigma[m] += 1
                            new_key = (tuple(new_sigma), tuple(new_rho))
                            if h_old * eta * W_cell * lambda_ != 0:
                                out.append((new_key, h_old * eta * W_cell * lambda_))
            return out
        
        def head(h_old, sigma, rho, cell, m):
            out = []
            for a, cell_a in enumerate(cells):
                for b, cell_b in enumerate(cells):
                    if rho[a * n_cells + b] <= 0:
                        continue
                    W_cell = cell_graph.get_cell_weight(cell)
                    lambda_ = cell_graph.get_two_table_with_suc_weight((cell_a, cell), suc_type=3) * \
                        (cell_graph.get_two_table_with_suc_weight((cell_a, cell), suc_type=1) ** (sigma[a] - 1)) * \
                        math.prod(1 if e == a else cell_graph.get_two_table_with_suc_weight((cell_s, cell), suc_type = 1) ** sigma[e] \
                        for e, cell_s in enumerate(cells))
                    new_rho = list(rho)
                    new_rho[a * n_cells + b] -= 1
                    new_rho[m * n_cells + b] += 1
                    new_sigma = list(sigma)
                    new_sigma[m] += 1
                    new_key = (tuple(new_sigma), tuple(new_rho))
                    if h_old * rho[a * n_cells + b] * W_cell * lambda_ != 0:
                        out.append((new_key, h_old * rho[a * n_cells + b] * W_cell * lambda_))
            return out
        
        def tail(h_old, sigma, rho, cell, m):
            out = []
            for a, cell_a in enumerate(cells):
                for b, cell_b in enumerate(cells):
                    if rho[a * n_cells + b] <= 0:
                        continue
                    W_cell = cell_graph.get_cell_weight(cell)
                    lambda_ = cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type=2) * \
                        (cell_graph.get_two_table_with_suc_weight((cell_b, cell), suc_type=1) ** (sigma[b] - 1)) * \
                        math.prod(1 if e == b else cell_graph.get_two_table_with_suc_weight((cell_s, cell), suc_type = 1) ** sigma[e] \
                        for e, cell_s in enumerate(cells))
                    new_rho = list(rho)
                    new_rho[a * n_cells + b] -= 1
                    new_rho[a * n_cells + m] += 1
                    new_sigma = list(sigma)
                    new_sigma[m] += 1
                    new_key = (tuple(new_sigma), tuple(new_rho))
                    if h_old * rho[a * n_cells + b] * W_cell * lambda_ != 0:
                        out.append((new_key, h_old * rho[a * n_cells + b] * W_cell * lambda_))
            return out
        
        def only(h_old, sigma, rho, cell, m):
            out = []
            W_cell = cell_graph.get_cell_weight(cell)
            lambda_ = math.prod(cell_graph.get_two_table_with_suc_weight((cell_s, cell), suc_type = 1) ** sigma[e] \
                        for e, cell_s in enumerate(cells))
            new_rho = list(rho)
            new_rho[m * n_cells + m] += 1
            new_sigma = list(sigma)
            new_sigma[m] += 1
            new_key = (tuple(new_sigma), tuple(new_rho))
            if (h_old * W_cell * lambda_ != 0):
                out.append((new_key, h_old * W_cell * lambda_))
            return out
        
        for i in range(1, domain_size - (0 if last_pred is None else 1)):
            old_table = table
            table = dict()
            if i >= (domain_size + 1) // 2:
                max_rho_size = domain_size - i
            else:
                max_rho_size = domain_size
            for j, cell in enumerate(cells) if first_pred is None and last_pred is None else body_cells:
                w = cell_graph.get_cell_weight(cell)
                for key, h_old in old_table.items():
                    sigma, rho = key
                    rho_size = sum(list(rho))
                    if rho_size > 1:
                    #if merge is possible, add merges
                        h_news = merge(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
                    
                    if rho_size <= max_rho_size:
                    #if |rho_old| = |rho_new| is possible, do head and tail
                        h_news = tail(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
                        h_news = head(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new

                    if rho_size < max_rho_size:
                    #if |rho_old| + 1 = |rho_new| is posible, do only
                    #e. g. at no time we need h(m, sigma, rho)
                    #where |rho| > domain_size - m + 1
                        h_news = only(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
        if last_pred is not None and domain_size > 1:
            max_rho_size = 1
            for j, cell in last_cells:
                old_table = table
                table = dict()
                w = cell_graph.get_cell_weight(cell)
                for key, h_old in old_table.items():
                    sigma, rho = key
                    rho_size = sum(list(rho))
                    if rho_size == 2:
                    #if merge is possible, add merges
                        h_news = merge(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
                    
                    if rho_size <= max_rho_size:
                    #if |rho_old| = |rho_new| is possible, do head and tail
                        h_news = tail(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
                        h_news = head(h_old, sigma, rho, cell, j)
                        for new_key, h_new in h_news:
                            table[new_key] = table.get(new_key, Rational(0, 1)) + h_new
        res = res + weight * sum(table.values())
                
    if leq_pred is not None:
            res *= Rational(math.factorial(domain_size), 1)
    return res

#