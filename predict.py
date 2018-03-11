import numpy as np
from sklearn import linear_model
from sklearn.externals import joblib
from sklearn.metrics import mean_squared_error
from sklearn.metrics import classification_report
import time
import sys

def predict(filename_testX, filename_model, filename_predict, filename_testy):
    start = time.time()
    print(start)

    clf = joblib.load(filename_model)

    print(clf.coef_)
    print(clf.intercept_)

    array_X = np.loadtxt(filename_testX, delimiter=',')
    array_y = np.loadtxt(filename_testy, delimiter=',')
    array_y_predict = clf.predict(array_X)

    print(mean_squared_error(array_y, array_y_predict))
    print(classification_report(array_y, array_y_predict))

    np.savetxt(filename_predict, array_y_predict, delimiter=',')

    end = time.time()
    print(end)

    print('training time :', str(end-start))

def main():
    filename_testX = sys.argv[1]
    filename_model = sys.argv[2]
    filename_predict = sys.argv[3]
    filename_testy = sys.argv[4]

    predict(filename_testX, filename_model, filename_predict, filename_testy)

if __name__ == '__main__':
    main()

