#!/bin/bash


py="pyenv/bin/python3.8"
source pyenv/bin/activate
path=$(realpath .)

echo "Testing ASP"
# asp
$py src/tester.py --asp --test-folder "${path}/tests/benchmarks/sequence" > "${path}/tests/results/sequences_asp.txt"
$py src/tester.py --asp --test-folder "${path}/tests/benchmarks/permutation" > "${path}/tests/results/permutations_asp.txt"
$py src/tester.py --asp --test-folder "${path}/tests/benchmarks/subset" > "${path}/tests/results/subsets_asp.txt"
$py src/tester.py --asp --test-folder "${path}/tests/benchmarks/multisubset" > "${path}/tests/results/multisubsets_asp.txt"
$py src/tester.py --asp --test-folder "${path}/tests/benchmarks/composition" > "${path}/tests/results/compositions_asp.txt"

echo "Testing #sat"
# #sat
$py src/tester.py --sat --test-folder "${path}/tests/benchmarks/sequence" > "${path}/tests/results/sequences_sat.txt"
$py src/tester.py --sat --test-folder "${path}/tests/benchmarks/permutation" > "${path}/tests/results/permutations_sat.txt"
$py src/tester.py --sat --test-folder "${path}/tests/benchmarks/subset" > "${path}/tests/results/subsets_sat.txt"
$py src/tester.py --sat --test-folder "${path}/tests/benchmarks/multisubset" > "${path}/tests/results/multisubsets_sat.txt"
$py src/tester.py --sat --test-folder "${path}/tests/benchmarks/composition" > "${path}/tests/results/compositions_sat.txt"

echo "Testing essence"
# essence
$py src/tester.py --essence --test-folder "${path}/tests/benchmarks/sequence" > "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence --test-folder "${path}/tests/benchmarks/permutation" > "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence --test-folder "${path}/tests/benchmarks/subset" > "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence --test-folder "${path}/tests/benchmarks/multisubset" > "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence --test-folder "${path}/tests/benchmarks/composition" > "${path}/tests/results/compositions_essence.txt"