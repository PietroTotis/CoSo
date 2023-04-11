# About 

This repository contains the code associated to the paper [Lifted Reasoning for Combinatorical Counting](https://jair.org/index.php/jair/article/view/14062). 

## Installation 

``install.sh`` installs the python environment with all the requirements listed in ``requirements.txt``.

``install-tools.sh`` installs optional third-party tools used in the experiments of the paper, but they are not a necessary requirement to run CoSo.

## Run

Activate the python environment pyenv and run
``python coso.py [filename]``, where ``filename`` is a file containing a CoLa problem (see below).

Add the option ``-v`` to generate a visual representation of the reasoning in CoSo. You can specify a filename (.html) if you want to save the output in a custom location.

# CoSo
CoSo is a lifted solver for combinatorics problems expressed in CoLa. Lifted reasoning is a form of automated reasoning which is based on the manipulation of groups of objects and their sizes to answer counting-related queries. CoSo can thus answer efficiently to problems expressed in CoLa such as:

> A kit of toy shapes contains five triangles and two squares. One triangle and one square are red. Another triangle and the other square are blue, and the remaining triangles are green. In how many different rows of four objects can the shapes be arranged if the two squares are included and the second object is green?

> Given the same set of shapes, in how many ways can the objects be divided into three(non-empty) groups such that the green objects all belong to the same group?

# CoLa
The input language for the solver is CoLa: a CoLa program is defined by a multiset of objects, a configuration, and optional constraints. 

The multiset of objects, the universe, can be expressed by enumeration, for example:
``shapes ={square_red, square_blue, triangle_red, triangle_blue, triangle_green, triangle_green, triangle_green};``.
Indistinguishable elements are denoted by the repetition of the same label as for ``triangle_green``.
Objects can have properties, which make them distinguishable: 
``property red ={triangle_red, square_red};``

``property blue ={triangle_blue, square_blue};``

``property green ={triangle_green};``

``property triangle ={triangle_red, triangle_blue, triangle_green};``

``property square ={square_blue square_red};``

Alternatively, the universe can be declared by listing the different properties and the corresponding number of objects with that property:

``property red;  property blue;  property green;``

``#red=2;  #blue=2;  #green=3;``

``property triangle; property square;#triangle=7; #square=2;``

``#square&red=1; #square&blue=1; #triangle&red=1;``

``#triangle&blue=1;  #triangle&green=3;``


Set formulas allow us to refer to combinations of property according to the usual set-operations: ``¬prop`` is the complement of a property or a set formula ``prop`` with respect to the universe. ``prop1 & prop2`` expresses the intersection, ``prop1 + prop2`` expresses the union of set formulas/properties. 

The keyword ``labelled`` can be prepended to the property declaration as a shortcut for expressingthat all objects having the property can be distinguished from one another, for instance:
``labelled property people;``
``#people=4;``



### Configuration.
A configuration is declared by using the following notation:
- ``[universe]`` denotes the set of all *permutations* of the universe.
- ``[repeated universe]`` denotes the set of all *sequences* of the universe.
- ``{universe}`` denotes the set of all *subsets* of the universe.
- ``{repeated universe}`` denotes the set of all *multisubsets* of the universe.
- ``[{universe}]`` denotes the set of all *compositions* of the universe.
- ``{{universe}}`` denotes the set of all *partitions* of the universe.

The set ``universe`` can be replaced with any name of a declared property.
To refer to a valid configuration, a label is associated with the keyword ``in``, for instance:

``row in [shapes];``

``groups in {{shapes}};``

## Constraints

Constraints restrict the set of valid configurations according to some additional requirement. In CoLa there are currently three types of constraints: size constraints, positional constraints, and counting constraints.

### Size Constraints
Size constraints define either the size of a configuration or a part (a subset in a partition/composition):

``#row = 4``

defines that only the rows of length 4 are valid. Any arithmetic comparison can be used: ``=,!=,<,<=,>,>=``.

### Positional Constraints.
Positional constraints are expressed with ``name[i] = prop`` where ``name`` is the name of the configuration, ``i`` is the position, ``prop`` is a property for example:

``row[2] = green;``

### Counting Constraints.
Counting constraints count either entities or parts that satisfy some property in the configuration. For example:

``#row & squares = 2``

states that the number of squares in a row is 2. Here we can also use any of ``=,!=,<,<=,>,>=``.
Counting constraints can be nested in partitions and compositions:

`` #(#part & green = 3) = 1;``

the keyword ``part`` denotes a generic subset of a partition/configuration.