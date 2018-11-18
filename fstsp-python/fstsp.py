# !/usr/bin/env python
"""
Compact formulation for the Flying Sidekick Traveling Salesman Problem (FSTSP)
employing O(|V|^5) variables. Note that this formulation employs the rules proposed
by Ponza (2016), except when the flag "murray_rules" is set, which applies the rules
from Murray et al (2015).
"""

__author__ = "Tulio Toffolo"
__copyright__ = "Copyright 2018, UFOP"

import time
from gurobipy import *
from arguments import Arguments
from instance import Instance
from solution import Solution

EPS = 1e-5


class Ell:
    """
    This class represents the 'moments' used within the formulation.
    """

    def __init__(self, inst):
        self.list = list(range(len(inst.V)))
        self.max_l = len(inst.V) - 1

    def list_ell(self, ell_):
        if ell_ < len(self.list):
            return self.list[max(0, ell_ - self.max_l):ell_]
        return []

    def list_ell_prime(self, ell):
        return self.list[ell + 1:min(len(self.list), ell + 1 + self.max_l)]


class CompactFormulation:
    """
    This class implements the proposed formulation using Gurobi.
    """

    def __init__(self, inst, big_m=None, murray_rules=False):
        if not big_m:
            big_m = inst.big_m if inst.big_m else 1000

        model = Model()
        L = Ell(inst)

        x = {(i, j, ell): model.addVar(obj=0, vtype=GRB.BINARY,
                                       name="x({i},{j}__{ell})".format(**locals()))
             for (i, j) in inst.A
             for ell in L.list}

        y = {(i, j, k, ell, ell_): model.addVar(obj=0,
                                                vtype=GRB.BINARY,
                                                name="y({i},{j},{k}__{ell},{ell_})".format(**locals()))
             for (i, j, k) in inst.D
             for ell in L.list
             for ell_ in L.list_ell_prime(ell)}

        t = [model.addVar(obj=0,
                          name="t({ell})".format(**locals()))
             for ell in L.list]
        t[0].ub = 0
        t.append(model.addVar(obj=1, name="t_total".format(**locals())))

        # depot constraints
        model.addConstr(quicksum(x[0, j, 0]
                                 for j in inst.V if (0, j) in inst.A)
                        == quicksum(x[j, 0, ell]
                                    for j in inst.V if (j, 0) in inst.A
                                    for ell in L.list[1:]),
                        "depot")
        model.addConstr(quicksum(x[0, j, 0]
                                 for j in inst.V if (0, j) in inst.A) == 1,
                        "depot_eq_1")

        # single customer visit
        for i in inst.V:
            model.addConstr(quicksum(x[i, j, ell]
                                     for j in inst.V if (i, j) in inst.A
                                     for ell in L.list) <= 1,
                            "single_moment_in({i})".format(**locals()))
            model.addConstr(quicksum(x[j, i, ell]
                                     for j in inst.V if (j, i) in inst.A
                                     for ell in L.list) <= 1,
                            "single_moment_out({i})".format(**locals()))

        # flow preservation constraints
        for k in inst.V_:
            for ell in L.list[1:]:
                model.addConstr(quicksum(x[j, k, ell - 1]
                                         for j in inst.V if (j, k) in inst.A)
                                == quicksum(x[k, j, ell]
                                            for j in inst.V if (k, j) in inst.A),
                                "flow({k},{ell})".format(**locals()))

        # one arc per moment constraints
        for ell in L.list:
            model.addConstr(quicksum(x[i, j, ell]
                                     for (i, j) in inst.A)
                            <= 1, "one_arc({ell})".format(**locals()))

        # one drone per moment constraints
        for ell in L.list:
            model.addConstr(quicksum(y[i, k, j, l, l_]
                                     for (i, k, j) in inst.D
                                     for l in L.list[:ell + 1]
                                     for l_ in L.list[ell + 1:l + 1 + L.max_l])
                            <= 1, "one_drone({ell})".format(**locals()))

        # all customers must be visited constraints
        for k in inst.V_:
            model.addConstr(quicksum(x[k, j, l]
                                     for j in inst.V if (k, j) in inst.A
                                     for l in L.list)
                            + quicksum(y[i, k, j, l, l_]
                                       for i in inst.V
                                       for j in inst.V if (i, k, j) in inst.D
                                       for l in L.list
                                       for l_ in L.list_ell_prime(l))
                            == 1, "all_customers({k})".format(**locals()))

        # drone launch constraints
        for i in inst.V:
            for ell in L.list:
                model.addConstr(quicksum(y[i, k, j, ell, l_]
                                         for k in inst.V_
                                         for j in inst.V if (i, k, j) in inst.D
                                         for l_ in L.list_ell_prime(ell))
                                <= quicksum(x[i, j, ell]
                                            for j in inst.V if (i, j) in inst.A),
                                "drone_launch({i},{ell})".format(**locals()))

        # drone return constraints
        for j in inst.V:
            for ell_ in L.list[1:]:
                model.addConstr(quicksum(y[i, k, j, l, ell_]
                                         for i in inst.V
                                         for k in inst.V_ if (i, k, j) in inst.D
                                         for l in L.list_ell(ell_))
                                <= quicksum(x[i, j, ell_ - 1]
                                            for i in inst.V if (i, j) in inst.A),
                                "drone_return({j},{ell_})".format(**locals()))

        # maximum drone endurance constraints
        for idx, ell in enumerate(L.list[1:]):
            if murray_rules:  # and idx == 0:
                continue
            for ell_ in L.list_ell_prime(ell):
                model.addConstr(t[ell_] - t[ell]
                                <= inst.E
                                + big_m * (1 - quicksum(y[i, k, j, ell, ell_]
                                                        for (i, k, j) in inst.D)),
                                "endurance({ell},{ell_})".format(**locals()))

        # time constraints considering only the truck and drone setup time
        for ell in L.list[1:] + [len(L.list)]:
            model.addConstr(t[ell] >= t[ell - 1]
                            + quicksum(inst.tau_truck[i][j] * x[i, j, ell - 1]
                                       for (i, j) in inst.A)
                            + inst.sl * (quicksum(y[i, k, j, ell - 1, ell_]
                                                  for (i, k, j) in inst.D for ell_ in L.list_ell_prime(ell - 1) if i != 0))
                            + inst.sr * (quicksum(y[i, k, j, ell_, ell]
                                                  for (i, k, j) in inst.D for ell_ in L.list_ell(ell))),
                            "time_truck_only({ell})".format(**locals()))

        # time constraints accounting drone's time
        for ell_ in L.list[1:]:
            for ell in L.list_ell(ell_):
                model.addConstr(t[ell_] >= t[ell]
                                + quicksum((inst.tau_drone[i][k] + inst.tau_drone[k][j] + (0 if murray_rules else inst.sl) + inst.sr)
                                           * y[i, k, j, ell, ell_]
                                           for (i, k, j) in inst.D),
                                "time_drones({ell},{ell_})".format(**locals()))

        # creating class variables
        self.inst = inst
        self.big_m = big_m
        self.murray_rules = murray_rules
        self.model = model
        self.L = L
        self.x = x
        self.y = y
        self.t = t
        self.solution = None

    def optimize(self, timelimit, sol=None):
        """
        Solves the compact formulation considering the time limit and the initial solution
        passed as arguments
        """
        if sol:
            # reading initial solution file
            initial_solution = Solution(self.inst, self.murray_rules)
            initial_solution.read(sol)

            # setting initial truck path
            i = initial_solution.truck_path[0]
            for ell, j in enumerate(initial_solution.truck_path[1:]):
                self.x[i, j, ell].start = 1
                i = j

            # setting initial drone paths
            for (i, j, k) in initial_solution.drone_paths:
                ell = initial_solution.truck_path.index(i)
                ell_ = initial_solution.truck_path.index(k) if k != 0 else len(initial_solution.truck_path) - 1
                self.y[i, j, k, ell, ell_].start = 1

        self.model.setParam("timelimit", timelimit)
        self.model.update()
        self.model.write("logs/fstsp.lp")
        self.model.optimize()
        self.model.write("logs/fstsp.sol")

        # creating final solution
        self.solution = Solution(self.inst, cost=self.model.objval, murray_rules=self.murray_rules)
        for key in [key for key in self.x.keys() if self.x[key].x > EPS]:
            self.solution.add_truck_arc((key[0], key[1]))
        for key in [key for key in self.y.keys() if self.y[key].x > EPS]:
            self.solution.add_drone_visit((key[0], key[1], key[2]))


if __name__ == "__main__":
    arg = Arguments(sys.argv)

    if not arg.E:
        inst = Instance(arg.instance, murray_rules=arg.murray_rules)
    else:
        inst = Instance(arg.instance, E=arg.E, murray_rules=arg.murray_rules)

    if arg.sol and arg.validate_only:
        solution = Solution(inst)
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
    formulation = CompactFormulation(inst, murray_rules=arg.murray_rules)
    formulation.optimize(arg.timelimit, sol=arg.sol)
    formulation.solution.print_solution()
    formulation.solution.write(arg.out)

    print("Total runtime: %.2f seconds" % (time.time() - start_time))
