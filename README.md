# CoSo
This repository contains the code associated to the paper Lifted Reasoning for Combinatorics Problems.
In order to reproduce the tests, presented in the paper, the script ``install_tools.sh`` downloads MiniZinc, creates a python 3.8 virtual environment with the ProbLog package and the necessary packages for CoSo ``run_tests.sh`` verifies that the requirements are installed and launches ``tester.py`` on the folders in ``tests/benchmarks``.

## CoLa
The input language for the solver is CoLa, the syntax is as follows:

### Domains. 
Domains can be expressed as intervals: ``set of entities = [1:n];`` or as enumeration: ``set of entities = {e1,e2,e3,e4,e5};``.
Indistinguishability is denoted by the keyword ``indist``:  ``set of indist entities = {e1,e2,e3,e4,e5};``.

### Configuration.
Configurations are expressed as ``conf in (! set)``, where ``( )`` is either ``[ ]`` (sequences/permutations) or ``{ }`` (subsets/multisubsets)
 and ``!`` is either ``|`` (permutations/subsets) or ``||`` (sequences/multisubsets).
 Partitions and compositions are expressed as ``conf in {{set}}`` and ``conf in [{set}]``

### Set formulas.
Set formulas are the following: ``¬ sf`` negates a set name or set formula ``sf``. ``sf1 & sf2`` expresses the intersection, ``sf1 + sf2`` expresses an union of set formulas (or sets). 

### Positional constraints.
Positional constraints are expressed with ``name[i] = sf`` where ``name`` is the name of the configuration, ``i`` is the position, ``sf`` is a set formula or an entity.

### Size constraints.
Size constraints are of the form ``#set op n`` where ``set`` is a name of a set or set formula, ``op`` is one of ``=,<,>,>=,=<,\=`` and ``n`` is a natural number.