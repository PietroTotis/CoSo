#!/bin/bash


py="pyenv/bin/python3.8"
source pyenv/bin/activate

# asp
$py src/tester.py --asp --test-folder tests/benhmarks/sequence > tests/results/sequences_asp.txt
$py src/tester.py --asp --test-folder tests/benhmarks/permutation > tests/results/permutations_asp.txt
$py src/tester.py --asp --test-folder tests/benhmarks/subset > tests/results/subsets_asp.txt
$py src/tester.py --asp --test-folder tests/benhmarks/multisubset > tests/results/multisubsets_asp.txt
$py src/tester.py --asp --test-folder tests/benhmarks/composition > tests/results/compositions_asp.txt

# #sat
$py src/tester.py --sat --test-folder tests/benhmarks/sequence > tests/results/sequences_sat.txt
$py src/tester.py --sat --test-folder tests/benhmarks/permutation > tests/results/permutations_sat.txt
$py src/tester.py --sat --test-folder tests/benhmarks/subset > tests/results/subsets_sat.txt
$py src/tester.py --sat --test-folder tests/benhmarks/multisubset > tests/results/multisubsets_sat.txt
$py src/tester.py --sat --test-folder tests/benhmarks/composition > tests/results/compositions_sat.txt

# essence
$py src/tester.py --essence --test-folder tests/benhmarks/sequence > tests/results/sequences_essence.txt
$py src/tester.py --essence --test-folder tests/benhmarks/permutation > tests/results/permutations_essence.txt
$py src/tester.py --essence --test-folder tests/benhmarks/subset > tests/results/subsets_essence.txt
$py src/tester.py --essence --test-folder tests/benhmarks/multisubset > tests/results/multisubsets_essence.txt
$py src/tester.py --essence --test-folder tests/benhmarks/composition > tests/results/compositions_essence.txt