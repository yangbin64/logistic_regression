import sys
import numpy as np
import pandas as pd

def read_X(filename_X):
    df_X = pd.read_csv(filename_X, header=None, index_col=None)
    return df_X

def read_y_bin(filename_y_bin):
    df_y_bin = pd.read_csv(filename_y_bin, header=None, names=['y_bin'], index_col=None)
    return df_y_bin

def read_y(filename_y):
    df_y = pd.read_csv(filename_y, header=None, names=['y'], index_col=None)
    return df_y

def concat_data(df_y_bin, df_y, df_X):
    df_all = pd.concat([df_y_bin, df_y, df_X], axis=1)
    return df_all

def split_df(df_all):
    df_all_0 = df_all[df_all['y_bin']==0]
    df_all_1 = df_all[df_all['y_bin']==1]
    print('aaaaaaaaaa')
    print(df_all_0)
    print('bbbbbbbbbb')
    print(df_all_1)
    return df_all_0, df_all_1

def output_csv_1(filename_output, df):
    df_Xy = df.drop('y_bin', axis=1)
    df_Xy.to_csv(filename_output, index=False, header=False)

def output_csv_2(filename_output_X, filename_output_y, df):
    df['y'].to_csv(filename_output_y, index=False, header=False)
    df_X = df.drop('y_bin', axis=1).drop('y', axis=1)
    df_X.to_csv(filename_output_X, index=False, header=False)

def main():
    filename_X = sys.argv[1]
    filename_y_bin = sys.argv[2]
    filename_y = sys.argv[3]
    filename_output = sys.argv[4]

    df_X = read_X(filename_X)
    df_y_bin = read_y_bin(filename_y_bin)
    df_y = read_y(filename_y)

    df_all = concat_data(df_y_bin, df_y, df_X)

    df_all_0, df_all_1 = split_df(df_all)

    #output_csv_1(filename_output + '_0.csv', df_all_0)
    #output_csv_1(filename_output + '_1.csv', df_all_1)

    output_csv_2(filename_output + '_X_0.csv', filename_output + '_y_0.csv', df_all_0)
    output_csv_2(filename_output + '_X_1.csv', filename_output + '_y_1.csv', df_all_1)

if __name__=='__main__':
    main()

