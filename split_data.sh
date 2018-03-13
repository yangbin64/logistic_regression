#!/bin/sh

SEP=6.0

python split_data.py ../20_input/train_X.csv ../30_train/train_y_predict_${SEP}.csv ../20_input/train_y.csv train

