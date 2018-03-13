from sklearn.linear_model import Ridge
from sklearn.externals import joblib
import numpy as np
import time
import sys

def train(filename_train_X, filename_train_y, filename_model, alpha):
    start = time.time()
    print(start)

    X = np.loadtxt(filename_train_X, delimiter=',')
    y = np.loadtxt(filename_train_y, delimiter=',')

    clf = Ridge(alpha)
    clf.fit(X, y)

    print(clf.coef_)
    print(clf.intercept_)

    joblib.dump(clf, filename_model)

    end = time.time()
    print(end)

    print('training time :', str(end-start))

def main():
    filename_train_X = sys.argv[1]
    filename_train_y = sys.argv[2]
    filename_model = sys.argv[3]
    alpha = float(sys.argv[4])

    train(filename_train_X, filename_train_y, filename_model, alpha)

if __name__ == '__main__':
    main()

