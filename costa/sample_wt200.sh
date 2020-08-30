#!/bin/bash
pushd ~/code/generic_wip
make
popd
../generic_wip/glfunc 11 300 /home/madjp/data/g/mf/mf.200 sample_wt200.input 101 
