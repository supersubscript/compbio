#!/bin/bash
for ii in $(seq 1 1000000); do top -c -n 1 | grep 'exonerate' >> inspiration.txt; sleep 1800; done
