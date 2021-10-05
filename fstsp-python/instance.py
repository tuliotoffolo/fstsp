#!/usr/bin/env python
"""
Representation of a Flying Sidekick Traveling Salesman Problem (FSTSP) instance.
"""
__author__ = "Tulio Toffolo and Julia Caria"
__copyright__ = "Copyright 2018, UFOP"

import importlib
import itertools
import os


class Instance:
    """
    This class represents an instance of the Flying Sidekick Traveling Salesman Problem (FSTSP).
    """

    def __init__(self, path, E=40, murray_rules=False):
        self.path = path
        self.V = []
        self.V_ = []
        self.V_drone = []
        self.tau_truck = []
        self.tau_drone = []
        self.A = []
        self.D = []
        self.big_m = None

        # drone launch, return and autonomy
        self.sl = 1
        self.sr = 1
        self.E = E

        # selecting the instance format: python or murray?
        if ".dat" in path:
            self.read_ponza_instance()
        else:
            self.read_murray_instance()

        # computing valid edges
        if not self.A:
            self.A = set([(i, j)
                          for i in self.V
                          for j in self.V if i != j or i == 0])

        # computing valid drone visits
        if not self.D:
            self.D = set([(i, k, j)
                          for i in self.V
                          for k in self.V_drone
                          for j in self.V
                          if self.tau_drone[i][k] + self.tau_drone[k][j] <=
                          self.E - (0 if murray_rules else self.sl) - self.sr
                          and i != k and k != j and (i != j or i == 0)])
        else:
            self.D = set([v for v in self.D
                          if self.tau_drone[v[0]][v[1]] + self.tau_drone[v[1]][v[2]] <=
                          self.E - (0 if murray_rules else self.sl) - self.sr])

        # computing sets for sub-route elimination constraints
        # self.sub = []
        # for subset in all_subsets(list(self.V_)):
        #     self.sub.append(subset)

    def read_murray_instance(self):
        """
        Reads instance data from a folder with the files specified by Murray and Chu (2015).
        """
        # reading nodes
        nodes = []
        nodes_csv = open(os.path.join(self.path, "nodes.csv"), "r")
        for line in nodes_csv.readlines():
            data = line.split(",")
            nodes.append({
                "id": int(data[0]),
                "x": float(data[1]),
                "y": float(data[2]),
                "drone_compatible": False if round(float(data[3])) else True,
            })
        nodes_csv.close()

        self.V = [node["id"] for node in nodes[:-1]]
        self.V_ = [node["id"] for node in nodes[1:-1]]
        self.V_drone = [node["id"] for node in nodes[1:-1] if node["drone_compatible"]]

        # reading truck travel times
        tau_csv = open(os.path.join(self.path, "tau.csv"), "r")
        for i, line in enumerate(tau_csv.readlines()):
            data = line.split(",")
            self.tau_truck.append([float(t) for t in data])
        tau_csv.close()

        # reading drone travel times
        tauprime_csv = open(os.path.join(self.path, "tauprime.csv"), "r")
        for i, line in enumerate(tauprime_csv.readlines()):
            data = line.split(",")
            self.tau_drone.append([float(t) for t in data])
        tauprime_csv.close()

    def read_ponza_instance(self):
        """
        Reads instance data from a '.dat' file specified by Ponza (2016).
        """
        dat_file = open(self.path, "r")
        dat_lines = dat_file.readlines()
        l = -1
        while l < len(dat_lines) - 1:
            l = l + 1
            line = dat_lines[l].replace(";", "").replace("\t", "  ")
            data = [v.strip() for v in line.split()]
            if not line or not data: continue

            # reading initial problem data
            if data[0] == "c":
                n = int(data[2])
                self.V = list(range(n + 1))
                self.V_ = list(range(1, n + 1))
            elif data[0] == "SL":
                self.sl = int(data[2])
            elif data[0] == "SR":
                self.sr = int(data[2])
            elif data[0] == "E":
                self.E = float(data[2])
            elif data[0] == "M":
                self.big_m = int(data[2])

            # reading drone customers
            elif data[0] == "Cdrones":
                data = line.replace("{", "").replace("}", "").split()
                self.V_drone = [int(v.strip()) for v in data[2:]]

            # reading distances
            elif data[0] in ("distance", "tauDrone", "tauTruck"):
                matrix_ = []
                for i in range(n + 2):
                    l = l + 1
                    line_ = dat_lines[l]
                    data_ = line_.replace("[", "").replace("]", "").split()
                    matrix_.append([float(v.strip()) for v in data_])

                if data[0] == "distance":
                    self.distances = matrix_
                elif data[0] == "tauDrone":
                    self.tau_drone = matrix_
                elif data[0] == "tauTruck":
                    self.tau_truck = matrix_

            # reading set F (equivalent to inst.D)
            elif data[0] == "F":
                matrix_ = []
                while "}" not in line:
                    l = l + 1
                    line = dat_lines[l]
                    data = line.replace("{", "").replace("}", "").replace(";", "").replace("\r", "") \
                        .replace("<", "").replace(">", "").replace(",", " ").split()
                    matrix_ += [int(v.strip()) for v in data]

                for i in range(0, len(matrix_), 3):
                    self.D.append((matrix_[i], matrix_[i + 1], matrix_[i + 2] if matrix_[i + 2] <= n else 0))


def all_subsets(ss):
    """
    Calculates all subsets of a given set
    """
    subsets = itertools.chain(*map(lambda x: itertools.combinations(ss, x),
                                   range(0, len(ss) + 1)))
    return [S for S in subsets if len(S) >= 2]
