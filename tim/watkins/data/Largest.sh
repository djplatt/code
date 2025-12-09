#!/bin/bash
grep -h Largest */*.out | sort -n -k 6 --key="8rn" | sort -s -k 6 -n -u
