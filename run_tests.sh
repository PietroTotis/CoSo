#!/bin/bash

if [[ ! (-d "minizinc" && -d "pyenv") ]]; then
    ./install_tools.sh
fi
py="pyenv/bin/python3.8"
source pyenv/bin/activate
"$py" src/tester.py --test-folder tests/examples 