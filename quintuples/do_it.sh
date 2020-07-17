#!/bin/bash
./quintuples4 2 216666668 > quint_2.out &
./quintuples4 216666669 433333335 > quint_216666669.out &
./quintuples4 433333336 650000002 > quint_433333336.out &
./quintuples4 650000003 866666669 > quint_650000003.out &
./quintuples4 866666670 1083333336 > quint_866666670.out &
./quintuples4 1083333337 1300000000 > quint_1083333337.out &
wait
