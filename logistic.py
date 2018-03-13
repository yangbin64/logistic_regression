import numpy as np
from sklearn import linear_model
from sklearn.externals import joblib
from sklearn.metrics import mean_squared_error
from sklearn.metrics import classification_report
import time
import sys

def loadX(filename_train_X):
    train_X = np.loadtxt(filename_train_X, dtype='int', delimiter=',')
    return train_X

def loady(filename_train_y):
    train_y = np.loadtxt(filename_train_y, dtype='int', delimiter=',')
    return train_y

def train(filename_model, train_X, train_y, filename_train_y_predict):
    start = time.time()
    print(start)

    clf = linear_model.LogisticRegression(C=0.01, class_weight='balanced')
    clf.fit(train_X, train_y)

    print(clf.coef_)
    print(clf.intercept_)
    #np.savetxt(filename_coef, clf.coef_, delimiter=',')
    #np.savetxt(filename_intercept, clf.intercept_, delimiter=',')

    joblib.dump(clf, filename_model)

    train_y_predict = clf.predict(train_X)
    print(mean_squared_error(train_y, train_y_predict))
    print(classification_report(train_y, train_y_predict))

    np.savetxt(filename_train_y_predict, train_y_predict, delimiter=',', fmt="%.0f")

    end = time.time()
    print(end)

    print('training time :', str(end-start))

def main():
    filename_train_X = sys.argv[1]
    filename_train_y = sys.argv[2]
    filename_model = sys.argv[3]
    filename_train_y_predict = sys.argv[4]

    train_X = loadX(filename_train_X)
    train_y = loady(filename_train_y)

    train(filename_model, train_X, train_y, filename_train_y_predict)

if __name__ == '__main__':
    main()

