#!/bin/bash


py="pyenv/bin/python3.8"
source pyenv/bin/activate
path=$(realpath .)

echo "Testing ASP"
# asp
# $py src/tester.py -t 300 --asp --test-folder "${path}/tests/benchmarks/sequence" &> "${path}/tests/results/sequences_asp.txt"
# $py src/tester.py -t 300 --asp --test-folder "${path}/tests/benchmarks/permutation" &> "${path}/tests/results/permutations_asp.txt"
# $py src/tester.py -t 300 --asp --test-folder "${path}/tests/benchmarks/subset" &> "${path}/tests/results/subsets_asp.txt"
# $py src/tester.py -t 300 --asp --test-folder "${path}/tests/benchmarks/multisubset" &> "${path}/tests/results/multisubsets_asp.txt"
# $py src/tester.py -t 300 --asp --test-folder "${path}/tests/benchmarks/composition" &> "${path}/tests/results/compositions_asp.txt"

echo "Testing #sat"
# #sat
# $py src/tester.py -t 300 --sat --test-folder "${path}/tests/benchmarks/sequence" &> "${path}/tests/results/sequences_sat.txt"
# $py src/tester.py -t 300 --sat --test-folder "${path}/tests/benchmarks/permutation" &> "${path}/tests/results/permutations_sat.txt"
# $py src/tester.py -t 300 --sat --test-folder "${path}/tests/benchmarks/subset" &> "${path}/tests/results/subsets_sat.txt"
# $py src/tester.py -t 300 --sat --test-folder "${path}/tests/benchmarks/multisubset" &> "${path}/tests/results/multisubsets_sat.txt"
# $py src/tester.py -t 300 --sat --test-folder "${path}/tests/benchmarks/composition" &> "${path}/tests/results/compositions_sat.txt"

echo "Testing essence"
# essence
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_10_5.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_15_5.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_20_5.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_15_10.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_20_10.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/sequence/sequence_20_15.test" >> "${path}/tests/results/sequences_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_10_5.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_15_5.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_20_5.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_15_10.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_20_10.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/permutation/permutation_20_15.test" >> "${path}/tests/results/permutations_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_10_5.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_15_5.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_20_5.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_15_10.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_20_10.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/subset/subset_20_15.test" >> "${path}/tests/results/subsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_10_5.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_15_5.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_20_5.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_15_10.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_20_10.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/multisubset/multisubset_20_15.test" >> "${path}/tests/results/multisubsets_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_10_5.test" >> "${path}/tests/results/compositions_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_15_5.test" >> "${path}/tests/results/compositions_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_20_5.test" >> "${path}/tests/results/compositions_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_15_10.test" >> "${path}/tests/results/compositions_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_20_15.test" >> "${path}/tests/results/compositions_essence.txt"
$py src/tester.py --essence -f "${path}/tests/benchmarks/composition/composition_20_15.test" >> "${path}/tests/results/compositions_essence.txt"