import numpy
from scipy.optimize import minimize
import Utils
import GreatFunctionForIterations

CoordinatesSystemArray1 = numpy.fromfile("D:\Progs\Bestfit\LOPAST'\SK11.txt", dtype=float, sep=' ')  #Файл со сферическими координатами
CoordinatesSystemArray2 = numpy.fromfile("D:\Progs\Bestfit\LOPAST'\SK22.txt", dtype=float, sep=' ')  #Файл с данными для перевода в декартову систему

Utils.coordinatesToDecart(CoordinatesSystemArray2, CoordinatesSystemArray1)

#Разворот вокруг оси OZ и сдвиг одной СК.
SK1 = numpy.fromfile("D:\Progs\Bestfit\LOPAST'\SK1.txt", dtype=float, sep=' ')
SK2 = numpy.fromfile("D:\Progs\Bestfit\LOPAST'\SK2.txt", dtype=float, sep=' ')

rowsQuantity = SK1.size / 3                                                                                             #Количество точек (и строк)
rowsQuantity = int(rowsQuantity)

SK1firstRow = [SK1[0], SK1[1], SK1[2]]
S1 = numpy.concatenate((SK1, SK1firstRow), axis=0)                                                                      #Первая вспомогательная матрица координат

SK2firstRow = [SK2[0], SK2[1], SK2[2]]
S2 = numpy.concatenate((SK2, SK2firstRow), axis=0)                                                                      #Вторая вспомогательная матрица координат

S1 = Utils.arrayToMatrix(S1)                                                                                            #В матрицы перевёл одномерные массивы
S2 = Utils.arrayToMatrix(S2)

anglesC = Utils.calculationOfTurningAngle(S1, S2, rowsQuantity)                                                         #Углы разворота для разных направлений

PreOffsetForDx = Utils.calculationPreOffsetForDx(S1, S2, anglesC, rowsQuantity)                                         #Предварительные значения смещений для Dx
PreOffsetForDy = Utils.calculationPreOffsetForDy(S1, S2, anglesC, rowsQuantity)                                         #Предварительные значения смещений для Dy
PreOffsetForDz = ((S2[0][2] - S1[0][2]) + (S2[1][2] - S1[1][2]) + (S2[2][2] - S1[2][2])) / rowsQuantity
                                                                                                                        #Предварительное значение смещения для Dz

avgC = numpy.sum(anglesC) / len(anglesC)                                                                                #Средний угол разворота
Bet = avgC
#Средние значения угла разворота и предварительных смещений по X и Y:
dX = numpy.sum(PreOffsetForDx) / len(PreOffsetForDx)
dY = numpy.sum(PreOffsetForDy) / len(PreOffsetForDy)
Dz = PreOffsetForDz

SKO = numpy.sqrt(Utils.functionSKO(dX, dY, S1, S2, Bet, rowsQuantity) / (rowsQuantity - 1))


def minimizeFuncSKO(x1, x2) :                                                                                           #Создание функции для последующей минимизации
    returnedValue = 0
    for i in range(0, rowsQuantity):
        returnedValue += \
            ((S2[i][0] - (S1[i][0] * numpy.cos(Bet) + (S1[i][1] * numpy.sin(Bet))) - x1) ** 2 +\
             (S2[i][1] - (S1[i][0] * (-numpy.sin(Bet)) + (S1[i][1] * numpy.cos(Bet))) - x2) ** 2)
    return returnedValue / (rowsQuantity - 1)


def  argumentsForSKOminimization(x) :
    return  minimizeFuncSKO(*x)

minimizedSKOvalues = minimize(argumentsForSKOminimization, (dX, dY))
#Минимизация,новые значения dX и dY. Я не знаю, что за магия тут происходит,так что я не выносил минимизацию в utils

SKO = numpy.sqrt(Utils.functionSKO(minimizedSKOvalues.x[0], minimizedSKOvalues.x[1], S1, S2, Bet, rowsQuantity) / (rowsQuantity - 1)) #SKO после минимизации

D = Utils.matrixD(minimizedSKOvalues.x[0], minimizedSKOvalues.x[1], Dz, rowsQuantity)                                   #Формирование матрицы,где в каждой строке начения dX, dY, Dz

NewReverseM = [[numpy.cos(Bet), -numpy.sin(Bet), 0],                                                                    #Новая матрица разворота
               [numpy.sin(Bet), numpy.cos(Bet), 0],
               [0, 0, 1]]
SK1 = numpy.array(SK1)
SK1 = numpy.reshape(SK1, (-1, 3))
SK1 = numpy.dot(SK1, NewReverseM) + D                                                                                   #Пересчёт координат для SK1

SK2 = numpy.array(SK2)                                                                                                  #Примитивно перевёл 1D массив в 2D
SK2 = numpy.reshape(SK2, (-1, 3))

SK1 = GreatFunctionForIterations.iterateFunc(SK1, SK2, rowsQuantity, 15)

SKOuResult = Utils.SKOu(SK1, SK2, rowsQuantity)
SKOkResult = Utils.SKOk(SK1, SK2, rowsQuantity)

numpy.savetxt("D:\Progs\Bestfit\LOPAST'\SK1isp.txt", SK1, fmt='%10.15f', delimiter=" ", newline="\r\n")
numpy.savetxt("D:\Progs\Bestfit\LOPAST'\SKOu.txt", SKOuResult, fmt='%10.15f', delimiter=" ", newline="\r\n")
numpy.savetxt("D:\Progs\Bestfit\LOPAST'\SKOk.txt", SKOkResult, fmt='%10.15f', delimiter=" ", newline="\r\n")

