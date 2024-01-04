## Flying Sidekick Traveling Salesman Problem (FSTSP)

### MIP formulations

This repository contains multiple MIP formulations for the problem. For more information, see the article [Exact and Heuristic Approaches to Drone Delivery Problems](https://arxiv.org/abs/2108.01996).

To cite our work, please use the following:

```bib
@article{Freitas2023,
    title = {Exact and heuristic approaches to Truck–Drone Delivery Problems},
    journal = {EURO Journal on Transportation and Logistics},
    volume = {12},
    pages = {100094},
    year = {2023},
    issn = {2192-4376},
    doi = {https://doi.org/10.1016/j.ejtl.2022.100094},
    url = {https://www.sciencedirect.com/science/article/pii/S219243762200019X},
    author = {J\'{u}lia C. Freitas and Puca Huachi V. Penna and T\'{u}lio A.M. Toffolo},
}
```

## Instance Files

Instance files in folder `input` were generated by **Ponza (2016)**. We thank the author for his kindness. He provided the files and granted us permission for publishing them online.

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
    python3 fstsp.py data/5a
    python3 murray_fstsp.py data/5b
    python3 fstsp.py data/5a -e 50 -m 100
    python3 fstsp.py data/20140810T123437v1 -e 20 -sol initial_solution.txt -out final_solution.txt -time 3600 
```

## References

- **Ponza, A. (2016)**. Optimization of drone-assisted parcel delivery. Ph.D. thesis, Universit a Degli Studi Di Padova.
