#!/bin/sh

SEP=6.0

python predict.py ../20_input/test_X.csv logistic_${SEP}.pkl test_y_predict_${SEP}.csv ../20_input/test_y_bin_${SEP}.csv

