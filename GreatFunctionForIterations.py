import numpy
from scipy.optimize import minimize
import Utils

#Вся эта функция является одной итерацией разворотов вокруг всех трёх осей.
def func(SK1, SK2, rowsQuantity) :
    U = Utils.anglesBetweenPointsA(SK1, SK2, rowsQuantity)  # Вычисление углов между радиус-векторами точек.
    U = numpy.array(U)  # Примитивно перевёл 1D массив в 2D
    U = numpy.reshape(U, (-1, 3))

    XA = Utils.TX(SK1, SK2, rowsQuantity, U)[0]  # Вывод из процедуры значений смещений и угла доворота
    YA = Utils.TX(SK1, SK2, rowsQuantity, U)[1]
    ZA = Utils.TX(SK1, SK2, rowsQuantity, U)[2]
    BetA = Utils.TX(SK1, SK2, rowsQuantity, U)[3]

    def minimizeSKOA(x1, x2, x3, x4):  # Создание функции для последующей минимизации
        returnedValue = 0
        for i in range(0, rowsQuantity):
            returnedValue += \
                (SK2[i][0] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MA(x4), 0)) - x1) ** 2 + \
                (SK2[i][1] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MA(x4), 1)) - x2) ** 2 + \
                (SK2[i][2] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MA(x4), 2)) - x3) ** 2
        return numpy.sqrt(returnedValue / (rowsQuantity - 1))

    def argumentsForSKOAminimization(x):
        return minimizeSKOA(*x)

    minimizedSKOAvalues = minimize(argumentsForSKOAminimization, (XA, YA, ZA, BetA))

    # Вывод значений матриц разворота и смещений:
    MA_ = \
    Utils.MDA(minimizedSKOAvalues.x[0], minimizedSKOAvalues.x[1], minimizedSKOAvalues.x[2], minimizedSKOAvalues.x[3],
              rowsQuantity)[0]
    DA = \
    Utils.MDA(minimizedSKOAvalues.x[0], minimizedSKOAvalues.x[1], minimizedSKOAvalues.x[2], minimizedSKOAvalues.x[3],
              rowsQuantity)[1]
    SK1 = numpy.dot(SK1, MA_) + DA

    U = Utils.anglesBetweenPointsA(SK1, SK2, rowsQuantity)  # Вычисление углов между радиус-векторами точек.
    U = numpy.array(U)  # Примитивно перевёл 1D массив в 2D
    U = numpy.reshape(U, (-1, 3))

    XB = Utils.TY(SK1, SK2, rowsQuantity, U)[0]  # Вывод из процедуры значений смещений и угла доворота
    YB = Utils.TY(SK1, SK2, rowsQuantity, U)[1]
    ZB = Utils.TY(SK1, SK2, rowsQuantity, U)[2]
    BetB = Utils.TY(SK1, SK2, rowsQuantity, U)[3]

    def minimizeSKOB(x1, x2, x3, x4):  # Создание функции для последующей минимизации
        returnedValue = 0
        for i in range(0, rowsQuantity):
            returnedValue += \
                (SK2[i][0] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MB(x4), 0)) - x1) ** 2 + \
                (SK2[i][1] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MB(x4), 1)) - x2) ** 2 + \
                (SK2[i][2] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MB(x4), 2)) - x3) ** 2
        return numpy.sqrt(returnedValue / (rowsQuantity - 1))

    def argumentsForSKOBminimization(x):
        return minimizeSKOB(*x)

    minimizedSKOBvalues = minimize(argumentsForSKOBminimization, (XB, YB, ZB, BetB))

    MB_ = \
    Utils.MDB(minimizedSKOBvalues.x[0], minimizedSKOBvalues.x[1], minimizedSKOBvalues.x[2], minimizedSKOBvalues.x[3],
              rowsQuantity)[0]
    DB = \
    Utils.MDB(minimizedSKOBvalues.x[0], minimizedSKOBvalues.x[1], minimizedSKOBvalues.x[2], minimizedSKOBvalues.x[3],
              rowsQuantity)[1]
    SK1 = numpy.dot(SK1, MB_) + DB

    U = Utils.anglesBetweenPointsA(SK1, SK2, rowsQuantity)  # Вычисление углов между радиус-векторами точек.
    U = numpy.array(U)  # Примитивно перевёл 1D массив в 2D
    U = numpy.reshape(U, (-1, 3))

    XC = Utils.TZ(SK1, SK2, rowsQuantity, U)[0]  # Вывод из процедуры значений смещений и угла доворота
    YC = Utils.TZ(SK1, SK2, rowsQuantity, U)[1]
    ZC = Utils.TZ(SK1, SK2, rowsQuantity, U)[2]
    BetC = Utils.TZ(SK1, SK2, rowsQuantity, U)[3]

    def minimizeSKOC(x1, x2, x3, x4):  # Создание функции для последующей минимизации
        returnedValue = 0
        for i in range(0, rowsQuantity):
            returnedValue += \
                (SK2[i][0] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MC(x4), 0)) - x1) ** 2 + \
                (SK2[i][1] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MC(x4), 1)) - x2) ** 2 + \
                (SK2[i][2] - numpy.dot([SK1[i][0], SK1[i][1], SK1[i][2]], Utils.column(Utils.MC(x4), 2)) - x3) ** 2
        return numpy.sqrt(returnedValue / (rowsQuantity - 1))

    def argumentsForSKOCminimization(x):
        return minimizeSKOC(*x)

    minimizedSKOCvalues = minimize(argumentsForSKOCminimization, (XC, YC, ZC, BetC))

    MC_ = \
    Utils.MDC(minimizedSKOCvalues.x[0], minimizedSKOCvalues.x[1], minimizedSKOCvalues.x[2], minimizedSKOCvalues.x[3],
              rowsQuantity)[0]
    DC = \
    Utils.MDC(minimizedSKOCvalues.x[0], minimizedSKOCvalues.x[1], minimizedSKOCvalues.x[2], minimizedSKOCvalues.x[3],
              rowsQuantity)[1]
    SK1 = numpy.dot(SK1, MC_) + DC

    return SK1


def iterateFunc(SK1, SK2, rowsQuantity, numberOfIterations) :
    finalSK1 = SK1
    for i in range(0, numberOfIterations) :
        finalSK1 = func(finalSK1, SK2, rowsQuantity)
    return finalSK1