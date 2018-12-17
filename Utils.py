import numpy
from scipy.optimize import minimize

def coordinatesToDecart(ArraySC2, ArraySC1):
    for i in range(len(ArraySC2)):
        if i == 0 or i % 3 == 0:
            ArraySC2[i] = (ArraySC1[0 + i]
                                 * numpy.sin(ArraySC1[1 + i] / 57.2958)
                                 * numpy.cos(ArraySC1[2 + i] / 57.2958))
        elif i == 1 or i % 3 == 1:
            ArraySC2[i] = (ArraySC1[0 + i]
                                 * numpy.sin(ArraySC1[1 + i] / 57.2958)
                                 * numpy.sin(ArraySC1[2 + i] / 57.2958))
        elif i == 2 or i % 3 == 2:
            ArraySC2[i] = (ArraySC1[0 + i]
                                 * numpy.cos(ArraySC1[1 + i] / 57.2958))


def arrayToMatrix(Array):
    matrix = ([[], [], [], []])
    for i in range(0, len(Array)):
        if i >= 0 and i < 3:
            matrix[0].append(Array[i])
        elif i >= 3 and i < 6:
            matrix[1].append(Array[i])
        elif i >= 6 and i < 9:
            matrix[2].append(Array[i])
        elif i >= 9 and i < 12:
            matrix[3].append(Array[i])
        if i == 11:
            return matrix


def findAngle(a, b, c, d) :
    return numpy.arctan(b / a) - numpy.arctan(d / c) + 6.28318


def calculationOfTurningAngle(S1array, S2array, rowsQuantity):          #Вычисление угла разворота (это Ci)
    turningAnglesC = []
    for i in range(0, rowsQuantity ):
        turningAnglesC.append(findAngle((S1array[i + 1][0] - S1array[i][0]), (S1array[i + 1][1] - S1array[i][1]),
                                              (S2array[i + 1][0] - S2array[i][0]), (S2array[i + 1][1] - S2array[i][1])))
    return turningAnglesC


def calculationPreOffsetForDx(S1array, S2array, arrayC, rowsQuantity):
    PreOffsetForDx = []
    for i in range(0, rowsQuantity):
        PreOffsetForDx.append(S2array[i][0] - ((S1array[i][0] * numpy.cos(arrayC[i])) + (S1array[i][1] * numpy.sin(arrayC[i]))))
    return PreOffsetForDx


def calculationPreOffsetForDy(S1array, S2array, arrayC, rowsQuantity):
    PreOffsetForDy = []
    for i in range(0, rowsQuantity):
        PreOffsetForDy.append(S2array[i][1] - ((S1array[i][0] * (-numpy.sin(arrayC[i]))) + (S1array[i][1] * numpy.cos(arrayC[i]))))
    return PreOffsetForDy


def functionSKO(dX, dY, S1array, S2array, Bet, rowsQuantity):
    numeratorsSum = 0
    for i in range(0, rowsQuantity):
        numeratorsSum += (S2array[i][0] - (S1array[i][0] * numpy.cos(Bet) + (S1array[i][1] * numpy.sin(Bet))) - dX) ** 2 +\
                         (S2array[i][1] - (S1array[i][0] * (-numpy.sin(Bet)) + (S1array[i][1] * numpy.cos(Bet))) - dY) ** 2
    return numeratorsSum


def matrixD(x1, x2, Dz, rowsQuantity) :
    matrix = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        matrix[i] = [x1, x2, Dz]
    return matrix


def anglesBetweenPointsA(SK1array, SK2array, rowsQuantity) :
    U = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        U[i].append(findAngle(SK1array[i][1], SK1array[i][2], SK2array[i][1], SK2array[i][2]))
        U[i].append(findAngle(SK1array[i][0], SK1array[i][2], SK2array[i][0], SK2array[i][2]))
        U[i].append(findAngle(SK1array[i][0], SK1array[i][1], SK2array[i][0], SK2array[i][1]))
    return U


def MA(A) :
    return numpy.array([[1, 0, 0],
            [0, numpy.cos(A), -numpy.sin(A)],
            [0, numpy.sin(A), numpy.cos(A)]])

def MB(B) :
    return numpy.array([[numpy.cos(B), 0, numpy.sin(B)],
            [0, 1, 0],
            [-numpy.sin(B), 0, numpy.cos(B)]])

def MC(C) :
    return numpy.array([[numpy.cos(C), -numpy.sin(C), 0],
            [numpy.sin(C), numpy.cos(C), 0],
            [0, 0, 1]])


def column(matrix, i):                                                            #Взять столбец матрицы
    return numpy.array([row[i] for row in matrix])


def TX(SK1array, SK2array, rowsQuantity, U) :                          #Вычисление смещений и среднего угла между радиус-векторами точек
    Bet = numpy.sum(column(U, 0)) / len(U)
    arrayForDxDyDz = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        arrayForDxDyDz[i].append(SK2array[i][0] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MA(Bet), 0)))

        arrayForDxDyDz[i].append(SK2array[i][1] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MA(Bet), 1)))

        arrayForDxDyDz[i].append(SK2array[i][2] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MA(Bet), 2)))
    dX = numpy.sum(column(arrayForDxDyDz, 0)) / len(arrayForDxDyDz)
    dY = numpy.sum(column(arrayForDxDyDz, 1)) / len(arrayForDxDyDz)
    dZ = numpy.sum(column(arrayForDxDyDz, 2)) / len(arrayForDxDyDz)
    return [dX, dY, dZ, Bet]


def TY(SK1array, SK2array, rowsQuantity, U) :                          #Вычисление смещений и среднего угла между радиус-векторами точек
    Bet = numpy.sum(column(U, 1)) / len(U)
    arrayForDxDyDz = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        arrayForDxDyDz[i].append(SK2array[i][0] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MB(Bet), 0)))

        arrayForDxDyDz[i].append(SK2array[i][1] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MB(Bet), 1)))

        arrayForDxDyDz[i].append(SK2array[i][2] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MB(Bet), 2)))
    dX = numpy.sum(column(arrayForDxDyDz, 0)) / len(arrayForDxDyDz)
    dY = numpy.sum(column(arrayForDxDyDz, 1)) / len(arrayForDxDyDz)
    dZ = numpy.sum(column(arrayForDxDyDz, 2)) / len(arrayForDxDyDz)
    return [dX, dY, dZ, Bet]


def TZ(SK1array, SK2array, rowsQuantity, U) :                          #Вычисление смещений и среднего угла между радиус-векторами точек
    Bet = numpy.sum(column(U, 2)) / len(U)
    arrayForDxDyDz = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        arrayForDxDyDz[i].append(SK2array[i][0] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MC(Bet), 0)))

        arrayForDxDyDz[i].append(SK2array[i][1] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MC(Bet), 1)))

        arrayForDxDyDz[i].append(SK2array[i][2] - numpy.dot([SK1array[i][0], SK1array[i][1], SK1array[i][2]],
                                 column(MC(Bet), 2)))
    dX = numpy.sum(column(arrayForDxDyDz, 0)) / len(arrayForDxDyDz)
    dY = numpy.sum(column(arrayForDxDyDz, 1)) / len(arrayForDxDyDz)
    dZ = numpy.sum(column(arrayForDxDyDz, 2)) / len(arrayForDxDyDz)
    return [dX, dY, dZ, Bet]


def MDA(a, b, c, d, rowsQuantity) :                                           #Процедура формирования матриц разворота и смещений
    M = MA(d)
    D = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        D[i] = [a, b, c]
    return [M, D]


def MDB(a, b, c, d, rowsQuantity) :                                           #Процедура формирования матриц разворота и смещений
    M = MB(d)
    D = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        D[i] = [a, b, c]
    return [M, D]


def MDC(a, b, c, d, rowsQuantity) :                                           #Процедура формирования матриц разворота и смещений
    M = MC(d)
    D = [[] for x in range(rowsQuantity)]
    for i in range(0, rowsQuantity) :
        D[i] = [a, b, c]
    return [M, D]


def SKOu(SK1array, SK2array, rowsQuantity) :
    arrayForABC = [[] for x in range(rowsQuantity)]
    sumForSKOa = 0
    sumForSKOb = 0
    sumForSKOc = 0
    for i in range(0, rowsQuantity) :
        arrayForABC[i].append(findAngle(SK1array[i][1], SK1array[i][2], SK2array[i][1], SK2array[i][2]))
        arrayForABC[i].append(findAngle(SK1array[i][0], SK1array[i][2], SK2array[i][0], SK2array[i][2]))
        arrayForABC[i].append(findAngle(SK1array[i][0], SK1array[i][1], SK2array[i][0], SK2array[i][1]))

    for i in range (0, rowsQuantity) :
        sumForSKOa += (arrayForABC[i][0] - (numpy.sum(column(arrayForABC, 0)) / len(arrayForABC))) ** 2
    SKOa = numpy.sqrt(sumForSKOa / (rowsQuantity - 1))

    for i in range(0, rowsQuantity):
        sumForSKOb += (arrayForABC[i][1] - (numpy.sum(column(arrayForABC, 1)) / len(arrayForABC))) ** 2
    SKOb = numpy.sqrt(sumForSKOb / (rowsQuantity - 1))

    for i in range(0, rowsQuantity):
        sumForSKOc += (arrayForABC[i][2] - (numpy.sum(column(arrayForABC, 2)) / len(arrayForABC))) ** 2
    SKOc = numpy.sqrt(sumForSKOc / (rowsQuantity - 1))

    return [SKOa, SKOb, SKOc]


def SKOk(SK1array, SK2array, rowsQuantity) :
    arrayForDxDyDz = [[] for x in range(rowsQuantity)]
    sumForSKOx = 0
    sumForSKOy = 0
    sumForSKOz = 0
    for i in range(0, rowsQuantity):
        arrayForDxDyDz[i].append(SK2array[i][0] - SK1array[i][0])
        arrayForDxDyDz[i].append(SK2array[i][1] - SK1array[i][1])
        arrayForDxDyDz[i].append(SK2array[i][2] - SK1array[i][2])

    for i in range(0, rowsQuantity):
        sumForSKOx += (arrayForDxDyDz[i][0] - (numpy.sum(column(arrayForDxDyDz, 0)) / len(arrayForDxDyDz))) ** 2
    SKOx = numpy.sqrt(sumForSKOx / (rowsQuantity - 1))

    for i in range(0, rowsQuantity):
        sumForSKOy += (arrayForDxDyDz[i][1] - (numpy.sum(column(arrayForDxDyDz, 1)) / len(arrayForDxDyDz))) ** 2
    SKOy = numpy.sqrt(sumForSKOy / (rowsQuantity - 1))

    for i in range(0, rowsQuantity):
        sumForSKOz += (arrayForDxDyDz[i][2] - (numpy.sum(column(arrayForDxDyDz, 2)) / len(arrayForDxDyDz))) ** 2
    SKOz = numpy.sqrt(sumForSKOz / (rowsQuantity - 1))

    return [SKOx, SKOy, SKOz]
