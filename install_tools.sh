#!/bin/bash

minizinc="https://github.com/MiniZinc/MiniZincIDE/releases/download/2.4.3/MiniZincIDE-2.4.3-bundle-linux-x86_64.tgz"
wget $minizinc
tar zxvf "MiniZincIDE-2.4.3-bundle-linux-x86_64.tgz" 
mv "MiniZincIDE-2.4.3-bundle-linux-x86_64"  "minizinc"
python3.8 -m venv pyenv
source ./pyenv/bin/activate
pip3 install -r requirements.txt
git clone https://github.com/ML-KULeuven/problog.git