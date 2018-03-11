import numpy as np
from sklearn import linear_model
from sklearn.externals import joblib
from sklearn.metrics import mean_squared_error
from sklearn.metrics import classification_report
import time
import sys

def loadX(filename_train_X):
    X = np.loadtxt(filename_train_X, dtype='int', delimiter=',')
    return X

def loady(filename_train_y):
    y = np.loadtxt(filename_train_y, dtype='int', delimiter=',')
    return y

def train(filename_model, X, y):
    start = time.time()
    print(start)

    clf = linear_model.LogisticRegression(C=0.01, class_weight='balanced')
    clf.fit(X, y)

    print(clf.coef_)
    print(clf.intercept_)
    #np.savetxt(filename_coef, clf.coef_, delimiter=',')
    #np.savetxt(filename_intercept, clf.intercept_, delimiter=',')

    joblib.dump(clf, filename_model)

    y_pred = clf.predict(X)
    print(mean_squared_error(y, y_pred))
    print(classification_report(y, y_pred))

    np.savetxt('y.txt', y, delimiter=',')
    np.savetxt('y_pred.txt', y_pred, delimiter=',')

    end = time.time()
    print(end)

    print('training time :', str(end-start))

def main():
    filename_train_X = sys.argv[1]
    filename_train_y = sys.argv[2]
    filename_model = sys.argv[3]

    X = loadX(filename_train_X)
    y = loady(filename_train_y)

    train(filename_model, X, y)

if __name__ == '__main__':
    main()

