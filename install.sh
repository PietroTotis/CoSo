#!/bin/bash

# install environment and python dependencies
echo "Installing python environment and packages"
python3.11 -m venv pyenv
source ./pyenv/bin/activate
pip3 install -q -r requirements.txt