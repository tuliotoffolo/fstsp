#!/usr/bin/env python
"""
Representation of a Flying Sidekick Traveling Salesman Problem (FSTSP) solution.
"""
from __future__ import print_function

__author__ = "Tulio Toffolo"
__copyright__ = "Copyright 2018, UFOP"

import os


class Solution:
    """
    This class represents a solution for the Flying Sidekick Traveling Salesman Problem (FSTSP)
    """

    def __init__(self, inst, cost=None, murray_rules=False):
        self.inst = inst
        self.murray_rules = murray_rules

        self.cost = cost
        self.truck_path = []
        self.drone_paths = []
        self.truck_arcs = set()
        self.drone_visits = set()

        self.computed = False

    def _compute_arcs(self):
        """
        Computes which arcs are traversed by truck and drone, assuming that their
        paths were already given.
        """
        if self.computed:
            return

        # computing truck arcs
        for i in range(1, len(self.truck_path)):
            self.truck_arcs.add((self.truck_path[i - 1], self.truck_path[i]))

        # computing drone visits
        for visit in self.drone_paths:
            self.drone_visits.add((visit[0], visit[1], visit[2]))

        self._compute_cost()
        self.computed = True

    def _compute_cost(self):
        """
        Computes the cost of this solution.
        """
        costs = [0] * (len(self.truck_path))

        # filling cost_truck list
        i = self.truck_path[0]
        t = 1
        for j in self.truck_path[1:]:
            # getting first the cost for the truck
            costs[t] = costs[t - 1] + self.inst.tau_truck[i][j]

            # getting SL additional cost
            if [key for key in self.drone_visits if key[0] == i and (not self.murray_rules or i != 0) and i != 0]:
                costs[t] += self.inst.sl

            # getting SR additional cost
            if [key for key in self.drone_visits if key[2] == j]:
                costs[t] += self.inst.sr

            # now comparing against cost given by drones
            for (k, l, m) in [key for key in self.drone_visits if key[2] == j]:
                ell = self.truck_path.index(k)
                costs[t] = max(costs[t], costs[ell] + self.inst.tau_drone[k][l] + self.inst.tau_drone[l][m]
                               + (0 if self.murray_rules else self.inst.sl) + self.inst.sr)

            # updating i and t vars
            i = j
            t = t + 1

        self.costs = costs
        self.cost = costs[-1]
        self.computed = True

    def _compute_paths(self):
        """
        Computes the paths of both truck and drone, assuming that truck arcs and
        drone visits were already given.
        :return:
        """
        if self.computed:
            return

        # computing truck path
        self.truck_path = [0] + [j for j in self.inst.V if (0, j) in self.truck_arcs]
        while self.truck_path[-1] != 0:
            i = self.truck_path[-1]
            self.truck_path += [j for j in self.inst.V if (i, j) in self.truck_arcs]

        # computing drone paths
        self.drone_paths = []
        for i in self.truck_path[:-1]:
            drone_path = [[i, k, j]
                          for k in self.inst.V_
                          for j in self.inst.V if (i, k, j) in self.drone_visits]
            if drone_path:
                self.drone_paths += drone_path

        self._compute_cost()
        self.computed = True

    def add_truck_arc(self, arc):
        """
        Includes an arc (i,j) traversed by the truck in the solution
        """
        i, j = arc
        self.truck_arcs.add((i, j))

    def add_drone_visit(self, arcs):
        """
        Includes a drone customer visit (i,k,j) in the solution.
        """
        i, k, j = arcs
        self.drone_visits.add((i, k, j))

    def print_solution(self):
        """
        Prints the solution value followed by the truck and drone's paths.
        """
        self._compute_paths()
        print("")
        print("    Solution cost = %.6f" % self.cost)

        # printing truck's path
        print("    Truck's path  = ", end="")
        for ell, i in enumerate(self.truck_path):
            if ell == len(self.truck_path) - 1:
                print("{i:2} ".format(**locals()), end="")
            else:
                print("{i:2}  - ".format(**locals()), end="")
        print("")

        # printing drone's paths
        for n, y in enumerate(self.drone_paths):
            i, k, j = y
            ell = self.truck_path.index(i)
            ell_ = self.truck_path[1:].index(j) + 1
            if n == 0:
                print("    Drone's paths = ", end="")
            else:
                print("                    ", end="")
            print((" " * 6 * ell + "{i:2} {k:2} " + " " * 6 * (ell_ - ell - 1) + "{j:2}").format(**locals()))
        print("")

    def read(self, path):
        """
        Reads a solution from the path given as argument.
        :param path: text file to read the solution from.
        """
        with open(path, "r") as f:
            import ast

            f.readline()  # ignoring instance name
            self.cost = float(f.readline().split()[-1])
            self.truck_path = ast.literal_eval(f.readline()[6:].strip())
            self.drone_paths = ast.literal_eval(f.readline()[6:].strip())

            self._compute_arcs()

    def write(self, path):
        """
        Writes the solution to the file given by argument
        :param path: text file to write the solution to
        """
        inst = self.inst
        with open(path, "w") as f:
            f.write("Instance: %s\n" % os.path.basename(self.inst.path))
            f.write("Cost: %.6f\n" % self.cost)

            # computing truck and drone paths
            self._compute_paths()

            # writing solution data
            f.write("Truck: %s\n" % str(self.truck_path))
            f.write("Drone: %s\n" % str(self.drone_paths))
            f.write("\n")

            # writing solution visual representation
            f.write("Visual representation:\n")

            # writing truck path
            f.write("    Truck: ")
            for ell, i in enumerate(self.truck_path):
                if ell == len(self.truck_path) - 1:
                    f.write("{i:2} ".format(**locals()))
                else:
                    f.write("{i:2}  - ".format(**locals()))
            f.write("\n")

            # writing drone path(s)
            for n, y in enumerate(self.drone_paths):
                i, k, j = y
                ell = self.truck_path.index(i)
                ell_ = self.truck_path[1:].index(j) + 1
                if n == 0:
                    f.write("    Drone: ")
                else:
                    f.write(" " * 11)
                f.write((" " * 6 * ell + "{i:2} {k:2} " + " " * 6 * (ell_ - ell - 1) + "{j:2}\n").format(**locals()))
            f.write("\n")

            # writing instance data
            f.write("Instance data:\n")
            f.write("    %-7s\t%10.4f\n" % ("Endur.", self.inst.E))
            f.write("    %-7s\t%10.4f\n" % ("SL", self.inst.sl))
            f.write("    %-7s\t%10.4f\n" % ("SR", self.inst.sr))
            f.write("\n")

            # writing details about the cost
            if not self.murray_rules:
                cost_sl_sr = sum([self.inst.sl for y in self.drone_visits] + [self.inst.sl for y in self.drone_visits])
                f.write("Cost components:\n")
                f.write("    %-7s\t%10.4f\n" % ("Arcs", self.cost - cost_sl_sr))
                f.write("    %-7s\t%10.4f\n" % ("SL+SR", cost_sl_sr))
                f.write("\n")

            # writing details about the whole trip
            f.write("Trip times:\n")
            for i in range(len(self.truck_path[1:])):
                f.write("    %-7s\t%10.4f\n" % ("T({i:2})".format(**locals()), self.costs[i]))
            f.write("    %-7s\t%10.4f\n" % ("Total ", self.costs[-1]))
            f.write("\n")

            # writing details about the truck trip
            f.write("Truck times (ignoring eventual waiting times):\n")
            i = self.truck_path[0]
            dist = 0
            for n, j in enumerate(self.truck_path[1:]):
                sl = inst.sl if [1 for key in self.drone_visits if key[0] == i and (not self.murray_rules or i != 0) and i != 0] else 0
                sr = inst.sr if [1 for key in self.drone_visits if key[2] == j] else 0
                dist += self.inst.tau_truck[i][j] + sl + sr

                if sl:
                    f.write("    %-7s\t%10.4f\n" % ("SL".format(**locals()), sl))
                f.write("    %-7s\t%10.4f\n" % ("[{i:2},{j:2}]".format(**locals()), self.inst.tau_truck[i][j]))
                if sr:
                    f.write("    %-7s\t%10.4f\n" % ("SR".format(**locals()), sr))
                i = j
            f.write("    %-7s\t%10.4f\n" % ("Total ", dist))
            f.write("\n")

            # writing details about the drone trip
            for n, y in enumerate(self.drone_paths):
                i, k, j = y
                f.write("Drone trip #{n}:\n".format(**locals()))
                if not self.murray_rules:
                    f.write("    %-7s\t%10.4f\n" % ("SL".format(**locals()), self.inst.sl))
                f.write("    %-7s\t%10.4f\n" % ("[{i:2},{k:2}]".format(**locals()), self.inst.tau_drone[i][k]))
                f.write("    %-7s\t%10.4f\n" % ("[{k:2},{j:2}]".format(**locals()), self.inst.tau_drone[k][j]))
                f.write("    %-7s\t%10.4f\n" % ("SR".format(**locals()), self.inst.sr))
                f.write("    %-7s\t%10.4f\n" % ("Total ", (0 if self.murray_rules else self.inst.sl)
                                                + self.inst.tau_drone[i][k] + self.inst.tau_drone[k][j]
                                                + self.inst.sr))
                f.write("\n")
