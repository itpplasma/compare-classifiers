#!/bin/bash

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

git clone https://github.com/itpplasma/SIMPLE
cd SIMPLE
git checkout compare-classifiers
./build.sh