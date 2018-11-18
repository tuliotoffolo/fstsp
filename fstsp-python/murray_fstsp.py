#!/usr/bin/env python
"""
Implementation of the compact formulation for the Flying Sidekick Traveling Salesman Problem
(FSTSP) proposed by Murray and Chu (2014).
"""
__author__ = "Tulio Toffolo and Julia Caria"
__copyright__ = "Copyright 2018, UFOP"

import time
from gurobipy import *
from arguments import Arguments
from instance import Instance
from solution import Solution

EPS = 1e-4


class MurrayFormulation:
    def __init__(self, inst, big_m=5000):
        model = Model()

        x = {(i, j): model.addVar(obj=0, vtype=GRB.BINARY,
                                  name='x({i},{j})'.format(**locals()))
             for (i, j) in inst.A}

        y = {(i, k, j): model.addVar(obj=0, vtype=GRB.BINARY,
                                     name='y({i},{k},{j})'.format(**locals()))
             for (i, k, j) in inst.D}

        u = {i: model.addVar(obj=0, vtype=GRB.INTEGER, lb=1, ub=len(inst.V_) + 2,
                             name="u({i})".format(**locals()))
             for i in inst.V}

        p = {(i, j): model.addVar(obj=0, vtype=GRB.BINARY, lb=1 if i == 0 else 0,
                                  name="p({i},{j})".format(**locals()))
             for (i, j) in inst.A}

        t = {i: model.addVar(obj=0 if i != inst.V[-1] else 1, lb=0,
                             name='t({i})'.format(**locals()))
             for i in inst.V}
        t[0].ub = 0

        t_ = {i: model.addVar(obj=0, lb=0,
                              name='t_prime({i})'.format(**locals()))
              for i in inst.V}
        t_[0].ub = 0

        model.addConstrs((quicksum(x[i, j]
                                   for i in inst.V[:-1] if (i, j) in inst.A) +
                          quicksum(y[i, j, k]
                                   for i in inst.V[:-1]
                                   for k in inst.V[1:] if (i, j, k) in inst.D)
                          == 1 for j in inst.V_), "c2")

        model.addConstr(quicksum(x[0, j]
                                 for j in inst.V[1:] if (0, j) in inst.A)
                        == 1, "c3")

        model.addConstr(quicksum(x[i, inst.V[-1]]
                                 for i in inst.V[:-1] if (i, inst.V[-1]) in inst.A)
                        == 1, "c4")

        model.addConstrs((u[i] - u[j] + 1 <= (len(inst.V_) + 2) * (1 - x[i, j])
                          for (i, j) in inst.A if i != 0), "c5")

        model.addConstrs((quicksum(x[i, j]
                                   for i in inst.V[:-1] if (i, j) in inst.A)
                          == quicksum(x[j, k]
                                      for k in inst.V[1:] if (j, k) in inst.A)
                          for j in inst.V_), "c6")

        model.addConstrs((quicksum(y[i, j, k]
                                   for j in inst.V_
                                   for k in inst.V[1:]
                                   if i != j != k and (i, j, k) in inst.D) <= 1
                          for i in inst.V[:-1]), "c7")

        model.addConstrs((quicksum(y[i, j, k]
                                   for i in inst.V[:-1]
                                   for j in inst.V_
                                   if i != j != k and (i, j, k) in inst.D) <= 1
                          for k in inst.V[1:]), "c8")

        model.addConstrs((2 * y[i, j, k]
                          <= quicksum(x[h, i]
                                      for h in inst.V[:-1] if (h, i) in inst.A)
                          + quicksum(x[l, k]
                                     for l in inst.V_ if (l, k) in inst.A)
                          for i in inst.V_
                          for j in inst.V_
                          for k in inst.V[1:] if (i, j, k) in inst.D), "c9")

        model.addConstrs((y[0, j, k] <= quicksum(x[h, k]
                                                 for h in inst.V[:-1] if (h, k) in inst.A)
                          for j in inst.V_
                          for k in inst.V[1:] if (0, j, k) in inst.D), "c10")

        model.addConstrs((u[k] - u[i]
                          >= 1 - (len(inst.V_) + 2) * (1 - quicksum(y[i, j, k]
                                                                    for j in inst.V_ if (i, j, k) in inst.D))
                          for i in inst.V_
                          for k in inst.V[1:] if (i, k) in inst.A), "c11")

        model.addConstrs((t_[i]
                          >= t[i] - big_m * (1 - quicksum(y[i, j, k]
                                                          for j in inst.V_
                                                          for k in inst.V[1:]
                                                          if i != j != k and (i, j, k) in inst.D))
                          for i in inst.V_), "c12")

        model.addConstrs((t_[i]
                          <= t[i] + big_m * (1 - quicksum(y[i, j, k]
                                                          for j in inst.V_
                                                          for k in inst.V[1:]
                                                          if i != j != k and (i, j, k) in inst.D))
                          for i in inst.V_), "c13")

        model.addConstrs((t_[k]
                          >= t[k] - big_m * (1 - quicksum(y[i, j, k]
                                                          for i in inst.V[:-1]
                                                          for j in inst.V_
                                                          if i != j != k and (i, j, k) in inst.D))
                          for k in inst.V[1:]), "c14")

        model.addConstrs((t_[k]
                          <= t[k] + big_m * (1 - quicksum(y[i, j, k]
                                                          for i in inst.V[:-1]
                                                          for j in inst.V_
                                                          if i != j != k and (i, j, k) in inst.D))
                          for k in inst.V[1:]),
                         "c15")

        model.addConstrs((t[k] >= t[h] + inst.tau_truck[h][k]
                          + inst.sl * (quicksum(y[k, l, m]
                                                for l in inst.V_
                                                for m in inst.V[1:]
                                                if k != l != m and (k, l, m) in inst.D))
                          + inst.sr * (quicksum(y[i, j, k]
                                                for i in inst.V[:-1]
                                                for j in inst.V_
                                                if i != j != k and (i, j, k) in inst.D))
                          - (big_m * (1 - x[h, k]))
                          for h in inst.V[:-1]
                          for k in inst.V[1:] if k != h), "c16")

        model.addConstrs((t_[j] >= t_[i] + inst.tau_drone[i][j]
                          - big_m * (1 - quicksum(y[i, j, k]
                                                  for k in inst.V[1:] if (i, j, k) in inst.D))
                          for j in inst.V_drone
                          for i in inst.V[:-1] if i != j), "c17")

        model.addConstrs((t_[k] >= t_[j] + inst.tau_drone[j][k] + inst.sr
                          - big_m * (1 - quicksum(y[i, j, k]
                                                  for i in inst.V[:-1] if (i, j, k) in inst.D))
                          for j in inst.V_drone
                          for k in inst.V[1:] if j != k), "c18")

        model.addConstrs((t_[k] - (t_[j] - inst.tau_drone[i][j])
                          <= inst.E + big_m * (1 - y[i, j, k])
                          for (i, j, k) in inst.D),
                         "c19")

        model.addConstrs((u[i] - u[j] >= 1 - (len(inst.V_) + 2) * p[i, j]
                          for i in inst.V_
                          for j in inst.V_ if j != i), "c20")

        model.addConstrs((u[i] - u[j] <= -1 + (len(inst.V_) + 2) * (1 - p[i, j])
                          for i in inst.V_
                          for j in inst.V_ if j != i), "c21")

        model.addConstrs((p[i, j] + p[j, i] == 1
                          for i in inst.V_
                          for j in inst.V_ if j != i), "c22")

        model.addConstrs((t_[l] >= t_[k]
                          - big_m * (3 - quicksum(y[i, j, k]
                                                  for j in inst.V_ if j != l and (i, j, k) in inst.D)
                                     - quicksum(y[l, m, n]
                                                for m in inst.V_
                                                for n in inst.V[1:]
                                                if m != i and m != k and m != l
                                                and n != i and n != k and (l, m, n) in inst.D)
                                     - p[i, l])
                          for i in inst.V[:-1]
                          for k in inst.V[1:]
                          for l in inst.V_
                          if k != i and l != i and l != k),
                         "c23")

        # creating class variables
        self.inst = inst
        self.big_m = big_m
        self.murray_rules = True
        self.model = model
        self.x = x
        self.y = y
        self.p = p
        self.u = u
        self.t = t
        self.t_ = t_
        self.solution = None

    def optimize(self, timelimit, sol=None):
        """
        Solves the compact formulation considering the time limit given as argument.
        """
        if sol:
            # reading initial solution file
            initial_solution = Solution(self.inst, murray_rules=self.murray_rules)
            initial_solution.read(sol)

            # setting initial truck path
            i = initial_solution.truck_path[0]
            for ell, j in enumerate(initial_solution.truck_path[1:]):
                j = j if j > 0 else self.inst.V[-1]
                self.x[i, j].start = 1
                i = j

            # setting initial drone paths
            for (i, j, k) in initial_solution.drone_paths:
                k = k if k > 0 else self.inst.V[-1]
                self.y[i, j, k].start = 1

        self.model.setParam('timelimit', timelimit)
        self.model.update()
        self.model.write('logs/murray_fstsp.lp')
        self.model.optimize()

        # creating final solution
        self.solution = Solution(self.inst, cost=self.model.objval, murray_rules=self.murray_rules)
        for key in [key for key in self.x.keys() if self.x[key].x > EPS]:
            self.solution.add_truck_arc((key[0], key[1] if key[1] != self.inst.V[-1] else 0))
        for key in [key for key in self.y.keys() if self.y[key].x > EPS]:
            self.solution.add_drone_visit((key[0], key[1], key[2] if key[2] != self.inst.V[-1] else 0))


if __name__ == "__main__":
    arg = Arguments(sys.argv, prefix="murray_")

    if not arg.E:
        inst = Instance(arg.instance)
    else:
        inst = Instance(arg.instance, E=arg.E)

    # quick instance adaptation for Murray's format (with an extra index for the depot)
    inst.V = inst.V + [len(inst.V)]
    inst.A = set([(i, j)
                  for i in inst.V[:-1]
                  for j in inst.V[1:] if i != j])
    inst.D = set([(i, j, k)
                  for i in inst.V[:-1]
                  for j in inst.V_drone
                  for k in inst.V[1:]
                  if inst.tau_drone[i][j] + inst.tau_drone[j][k] <= inst.E - inst.sr
                  and i != j != k])

    if arg.sol and arg.validate_only:
        solution = Solution(inst, murray_rules=True)
        solution.read(arg.sol)
        solution.print_solution()
        solution.write(arg.out)
        exit(0)

    start_time = time.time()

    # creating folders 'logs' and 'solutions' if they don't exist already
    if not os.path.exists('logs'):
        os.makedirs('logs')
    if not os.path.exists('solutions'):
        os.makedirs('solutions')

    # creating and solving the formulation using Gurobi
    formulation = MurrayFormulation(inst)
    formulation.optimize(arg.timelimit, sol=arg.sol)
    formulation.solution.print_solution()
    formulation.solution.write(arg.out)

    print("Total runtime: %.2f seconds" % (time.time() - start_time))
