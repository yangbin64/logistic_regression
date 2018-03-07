import sys
import pandas as pd
import numpy as np

COLUMN=147

def loadData(filename):
    df = pd.read_csv(filename, header=None, index_col=None)
    return df

def seperate(df, filename_Xy, sep):
    #df_new = pd.DataFrame(data=np.zeros((round(len(df.index)/10), COLUMN+1), dtype=int))
    array_new = np.zeros((round(len(df.index)/50), COLUMN+1), dtype=int)
    count=0
    for key,row in df.iterrows():
        count=count+1
        if count%100000==0:
            print(count)

        if key >= array_new.shape[0]:
            break

        if df.ix[key,0] < sep:
            array_new[key, 0] = 0
        else:
            array_new[key, 0] = 1

        for i in np.arange(5):
            col_x = int(row[i+1]) + 1
            array_new[key, col_x] = 1

    df_new = pd.DataFrame(data=array_new)
    print(df_new)
    return df_new

def main(filename, filename_Xy, sep):
    df = loadData(filename)
    df_new = seperate(df, filename_Xy, sep)
    df_new.to_csv(filename_Xy, header=False, index=False)

if __name__=='__main__':
    filename = sys.argv[1]
    filename_Xy = sys.argv[2]
    sep = float(sys.argv[3])
    main(filename, filename_Xy, sep)
