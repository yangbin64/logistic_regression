import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from sklearn.metrics import classification_report
import time
import sys

def loadXy(filename_train):
    Xy = np.loadtxt(filename_train, dtype='int', delimiter=',')
    y = Xy[:,0]
    X = Xy[:,1:]
    #print(X)
    #print(y)
    return X, y

def train(filename_coef, X, y):
    start = time.time()
    print(start)

    clf = linear_model.LogisticRegression(C=0.01)
    clf.fit(X, y)
    print(clf.coef_)
    np.savetxt(filename_coef, clf.coef_, delimiter=',')

    y_pred = clf.predict(X)
    print(mean_squared_error(y, y_pred))
    print(classification_report(y, y_pred))

    np.savetxt('y.txt', y, delimiter=',')
    np.savetxt('y_pred.txt', y_pred, delimiter=',')

    end = time.time()
    print(end)

    print('training time :', str(end-start))

def main():
    filename_train = sys.argv[1]
    filename_coef = sys.argv[2]

    X, y = loadXy(filename_train)
    #train(filename_train, filename_coef, alpha, n_samples, n_features)

    train(filename_coef, X, y)

if __name__ == '__main__':
    main()

