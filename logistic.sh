#!/bin/sh

SEP=6.0

python logistic.py ../20_input/train_X.csv ../20_input/train_y_bin_${SEP}.csv logistic_${SEP}.pkl

