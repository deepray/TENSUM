#!/bin/bash
#BSUB -W 4:00
#BSUB -n 16
tensum_run.py input.param > output.out
