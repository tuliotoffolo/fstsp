import os


class Arguments:
    """
    This simple class stores data passed by the user via command-line arguments
    """

    def __init__(self, argv, prefix=""):
        if len(argv) < 2:
            usage(argv)

        # initializing default values
        self.instance = argv[1]
        self.big_m = None
        self.E = None
        self.murray_rules = False
        self.sol = None
        self.out = None
        self.timelimit = 900
        self.validate_only = False

        # reading command-line arguments
        i = 2
        while i < len(argv):
            if argv[i].lower() == '-e':
                i += 1
                self.E = float(argv[i])
            elif argv[i].lower() == '-m':
                i += 1
                self.big_m = float(argv[i])
            elif argv[i].lower() == '-murray':
                self.murray_rules = True
            elif argv[i].lower() == '-sol':
                i += 1
                self.sol = argv[i]
            elif argv[i].lower() == '-out':
                i += 1
                self.out = argv[i]
            elif argv[i].lower() == '-time':
                i += 1
                self.timelimit = float(argv[i])
            elif argv[i].lower() == '-validate':
                self.validate_only = True
            else:
                usage(argv)

            i += 1

        if not self.out:
            self.out = os.path.join("solutions", prefix + os.path.basename(argv[1])) + \
                       (("_e%d" % self.E) if self.E else "") + ".txt"
            print(self.out)


def usage(argv):
    """
    Prints the usage of the software, including all possible arguments.
    :param argv: arguments passed by the user.
    """
    exe = os.path.basename(argv[0])
    print("")
    print("Usage: gurobi.sh " + exe + " <instance_file> [options]")
    print("    <instance_file>   : path of instance input file(s).")
    print("")
    print("Options:")
    print("    -e <drone_endurance> : specify drone battery endurance (default: 40).")
    print("    -m <big_m_value>     : specify maximum objective value (default: 1000 or instance's value).")
    print("    -murray              : use murray rules within formulation.")
    print("    -sol <solution_file> : read initial solution.")
    print("    -out <solution_file> : specify output solution file name.")
    print("    -time <time_limit>   : runtime limit in seconds (default: 900).")
    print("    -validate            : use formulation only to validate initial solution.")
    print("")
    exit(1)
