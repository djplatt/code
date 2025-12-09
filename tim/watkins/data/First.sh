#!/bin/bash
grep -h First */*.out | sort -n -k 4 --key="9rn" | sort -s -k 4 -n -u
