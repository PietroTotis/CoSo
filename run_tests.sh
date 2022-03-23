#!/bin/bash


py="pyenv/bin/python3.8"
source pyenv/bin/activate
path=realpath .

# asp
$py src/tester.py --asp --test-folder tests/benhmarks/sequence > "$path/tests/results/sequences_asp.txt"
$py src/tester.py --asp --test-folder tests/benhmarks/permutation > "$path/tests/results/permutations_asp.txt"
$py src/tester.py --asp --test-folder tests/benhmarks/subset > "$path/tests/results/subsets_asp.txt"
$py src/tester.py --asp --test-folder tests/benhmarks/multisubset > "$path/tests/results/multisubsets_asp.txt"
$py src/tester.py --asp --test-folder tests/benhmarks/composition > "$path/tests/results/compositions_asp.txt"

# #sat
$py src/tester.py --sat --test-folder tests/benhmarks/sequence > "$path/tests/results/sequences_sat.txt"
$py src/tester.py --sat --test-folder tests/benhmarks/permutation > "$path/tests/results/permutations_sat.txt"
$py src/tester.py --sat --test-folder tests/benhmarks/subset > "$path/tests/results/subsets_sat.txt"
$py src/tester.py --sat --test-folder tests/benhmarks/multisubset > "$path/tests/results/multisubsets_sat.txt"
$py src/tester.py --sat --test-folder tests/benhmarks/composition > "$path/tests/results/compositions_sat.txt"

# essence
$py src/tester.py --essence --test-folder tests/benhmarks/sequence > "$path/tests/results/sequences_essence.txt"
$py src/tester.py --essence --test-folder tests/benhmarks/permutation > "$path/tests/results/permutations_essence.txt"
$py src/tester.py --essence --test-folder tests/benhmarks/subset > "$path/tests/results/subsets_essence.txt"
$py src/tester.py --essence --test-folder tests/benhmarks/multisubset > "$path/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence --test-folder tests/benhmarks/composition > "$path/tests/results/compositions_essence.txt"