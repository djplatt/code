#!/bin/bash
awk '/1     / {print $5 " " $6 " " $7}' ${1} | arb_sum
awk '/2     / {print $5 " " $6 " " $7}' ${1} | arb_sum
awk '/kewes on file/ {print $13 " " $14 " " $15}' ${1} | arb_sum
