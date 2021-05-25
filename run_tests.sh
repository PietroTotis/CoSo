#!/bin/bash


py="pyenv/bin/python3.8"
source pyenv/bin/activate

$py src/tester.py --test-folder tests/examples

echo "h730: distinguishable positions"
time "$py" src/tester.py -f tests/experiments/dist_h730.pl
echo "h730: indistinguishable positions"
time "$py" src/tester.py -f tests/experiments/ind_h730.pl
echo "h830: distinguishable positions"
time "$py" src/tester.py -f tests/experiments/dist_h830.pl
echo "h830: indistinguishable positions"
time "$py" src/tester.py -f tests/experiments/ind_h830.pl
echo "l305: distinguishable positions"
time "$py" src/tester.py -f tests/experiments/dist_l305.pl
echo "l305: indistinguishable positions"
time "$py" src/tester.py -f tests/experiments/ind_l305.pl

echo ""
echo "Testing domain sizes..."
echo ""

echo "h730: domain size 1x"
time "$py" src/tester.py -f tests/experiments/domains_h730/h730_1.pl
echo "h730: domain size 10x"
time "$py" src/tester.py -f tests/experiments/domains_h730/h730_10.pl
echo "h730: domain size 100x"
time "$py" src/tester.py -f tests/experiments/domains_h730/h730_100.pl
echo "h730: domain size 1000x"
time "$py" src/tester.py -f tests/experiments/domains_h730/h730_1000.pl
echo "h730: domain size 10000x"
time "$py" src/tester.py -f tests/experiments/domains_h730/h730_10000.pl


echo "h830: domain size 1x"
time "$py" src/tester.py -f tests/experiments/domains_h830/h830_1.pl
echo "h830: domain size 10x"
time "$py" src/tester.py -f tests/experiments/domains_h830/h830_10.pl
echo "h830: domain size 100x"
time "$py" src/tester.py -f tests/experiments/domains_h830/h830_100.pl
echo "h830: domain size 1000x"
time "$py" src/tester.py -f tests/experiments/domains_h830/h830_1000.pl
echo "h830: domain size 10000x"
time "$py" src/tester.py -f tests/experiments/domains_h830/h830_10000.pl
# echo "h830: domain size 1000000x"
# time "$py" src/tester.py -f tests/experiments/domains_h830/h830_1000000.pl

echo "l305: domain size 1x"
time "$py" src/tester.py -f tests/experiments/domains_l305/l305_1.pl
echo "l305: domain size 10x"
time "$py" src/tester.py -f tests/experiments/domains_l305/l305_10.pl
echo "l305: domain size 100x"
time "$py" src/tester.py -f tests/experiments/domains_l305/l305_100.pl
echo "l305: domain size 1000x"
time "$py" src/tester.py -f tests/experiments/domains_l305/l305_1000.pl
echo "l305: domain size 10000x"
time "$py" src/tester.py -f tests/experiments/domains_l305/l305_10000.pl
# echo "l305: domain size 1000000x"
# time "$py" src/tester.py -f tests/experiments/domains_l305/l305_1000000.pl

echo "END TEST"
echo ""