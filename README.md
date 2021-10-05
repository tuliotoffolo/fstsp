## Flying Sidekick Traveling Salesman Problem (FSTSP)

### MIP formulations

This repository contains multiple MIP formulations for the problem. For more information, see the article [Exact and Heuristic Approaches to Drone Delivery Problems](https://arxiv.org/abs/2108.01996).

To cite our work, please use the following:

```bib
@misc{freitas2021exact,
    title={Exact and Heuristic Approaches to Drone Delivery Problems}, 
    author={J\'ulia C. Freitas and Puca Huachi V. Penna and T\'ulio A. M. Toffolo},
    year={2021},
    eprint={2108.01996},
    archivePrefix={arXiv},
    primaryClass={cs.AI}
}
```

## How to run

To run the code, use the syntax presented below:

```
Usage: gurobi.sh <formulation.py> <instance_file> [options]
    <formulation.py>     : formulation python file.
    <instance_file>      : path of instance input file(s).

Options:
    -e <drone_endurance> : specify drone battery endurance (default: 40).
    -m <big_m_value>     : specify maximum objective value (default: 5000).
    -murray              : use murray rules within formulation.
    -sol <solution_file> : read initial solution.
    -out <solution_file> : specify output solution file name.
    -time <time_limit>   : runtime limit in seconds (default: 900).
    -validate            : use formulation only to validate the initial solution.
    
Examples:
    gurobi.sh fstsp.py data/5a
    gurobi.sh murray_fstsp.py data/5b
    gurobi.sh fstsp.py data/5a -e 50 -m 100
    gurobi.sh fstsp.py data/20140810T123437v1 -e 20 -sol initial_solution.txt -out final_solution.txt -time 3600 
```
