# CoSo
This repository contains the code associated to the paper Lifted Reasoning for Combinatorics Problems.
In order to reproduce the tests, presented in the paper, the script ``install_tools.sh`` downloads MiniZinc, creates a python 3.8 virtual environment with the ProbLog package and the necessary packages for CoSo ``run_tests.sh`` verifies that the requirements are installed and launches ``tester.py`` on the folders in ``tests/paper``.

## CoLa
The input language for the solver is CoLa, the syntax is as follows:

### Domains. 
Domains can be expressed as intervals: ``student([1,n]).`` or as enumeration: ``student(5,s1,s2,s3,s4,s5).``

### Structure.
Structures are expressed as ``structure(name,type,variant,domain)``, where ``name`` is the label associated to the combinatorial structure, type is either ``sequence`` or ``subset``, type is ``false`` for the aforementioned kind of structures, ``true`` if the sequence is a permutation or the subset is a multi-subset.

``size(name,n)`` declares the size of the combinatorial structure.

### Domain formulas.
Domain formulas are the following: ``not(df)`` negates a domain name or domain formula ``df``. ``inter(df1,df2)`` expresses the intersection, ``union(df1,df2)`` expresses an union of domain formulas (or domains). 

### Choice constraints.
Choice constraints are expressed with ``pos(name,n,df)`` where ``name`` is the name of the structure, ``n`` is the position, ``df`` is a domain formula. ``in(e)`` states that element ``e`` from one of the declared domains belongs to the subset.

### Count constraints.
Count constraints are of the form ``count(name,expr)`` where ``name`` is the name of the structure and ``expr`` is an expression of the form ``df op n`` where ``df`` is a domain formula, ``op`` is one of ``=,<,>,>=,=<,\=`` and ``n`` is a natural number.