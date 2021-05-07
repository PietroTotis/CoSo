#!/bin/bash

py="pyenv/bin/python3.8"
source pyenv/bin/activate
echo "h730: distinguishable positions"
"$py" src/tester.py -f tests/experiments/dist_h273.pl
echo "h730: indistinguishable positions"
"$py" src/tester.py -f tests/experiments/ind_h730.pl
echo "h830: distinguishable positions"
"$py" src/tester.py -f tests/experiments/dist_h830.pl
echo "h830: indistinguishable positions"
"$py" src/tester.py -f tests/experiments/ind_h830.pl
echo "h849: distinguishable positions"
"$py" src/tester.py -f tests/experiments/dist_h849.pl
echo "h849: indistinguishable positions"
"$py" src/tester.py -f tests/experiments/ind_h849.pl