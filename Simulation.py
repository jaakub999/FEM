from Jacobian import *

class SOE:
    def __init__(self, data, element, node):
        HG = np.zeros((data.nN, data.nN))
        CG = np.zeros((data.nN, data.nN))
        PG = np.zeros(data.nN)
        HR = np.zeros((data.nN, data.nN))
        PR = np.zeros(data.nN)
        time = 0

        for ITERATION in range(int(data.simulationTime/data.stepTime)):
            jacobian = Jacobian(data.npc)
            for i in range(data.nN):
                for j in range(data.nN):
                    HG[i][j] = 0
                    CG[i][j] = 0
                    HR[i][j] = 0
                PG[i] = 0
                PR[i] = 0

            for smallITERATION in range(data.nE):
                for i in range(4):
                    for j in range(4):
                        element[smallITERATION].H[i][j] = 0
                        element[smallITERATION].C[i][j] = 0
                        element[smallITERATION].HBC[i][j] = 0
                    element[smallITERATION].P[i] = 0

                jacobian.calculate_H_and_C(element[smallITERATION], node, data)

                for side in range(4):
                    if node[int(element[smallITERATION].ID[side]) - 1].BC == 1 and node[int(element[smallITERATION].ID[(side + 1) % 4]) - 1].BC == 1:
                        jacobian.calculate_HBC_and_P(element[smallITERATION], node, data, side)

                for row in range(len(element[smallITERATION].H)):
                    for column in range(len(element[smallITERATION].H[row])):
                        element[smallITERATION].H[row][column] += element[smallITERATION].HBC[row][column]

            for i in range(data.nE):
                for j in range(4):
                    for k in range(4):
                        HG[int(element[i].ID[j] - 1)][int(element[i].ID[k] - 1)] += element[i].H[j][k]
                        CG[int(element[i].ID[j] - 1)][int(element[i].ID[k] - 1)] += element[i].C[j][k]
                    PG[int(element[i].ID[j] - 1)] += element[i].P[j]

            for i in range(data.nN):
                for j in range(data.nN):
                    HR[i][j] = HG[i][j] + (CG[i][j] / data.stepTime)

            for i in range(data.nN):
                PR[i] = -PG[i]

                for j in range(data.nN):
                    PR[i] += (CG[i][j] / data.stepTime) * node[j].t0

            t0 = elimination(data, HR, PR)
            max = 0
            min = 1000

            for i in range(data.nN):
                if t0[i] > max:
                    max = t0[i]

                if t0[i] < min:
                    min = t0[i]

                node[i].t0 = t0[i]
            time += data.stepTime

            print("{:>7}".format(time), end = "")
            print("{:>15}".format(round(min, 3)), end = "")
            print("{:>15}".format(round(max, 3)))


def elimination(data, HR, PR):
    vector = np.zeros(data.nN)
    matrix = np.zeros((data.nN, data.nN + 1))

    for i in range(data.nN):
        for j in range(data.nN + 1):
            if j < data.nN:
                matrix[i][j] = HR[i][j]

            else:
                matrix[i][j] = PR[i]

    for i in range(data.nN):
        for pivot in range(i+1, data.nN):
            if matrix[i][i] < matrix[pivot][i]:
                for j in range(data.nN + 1):
                    temp = matrix[pivot][j]
                    matrix[pivot][j] = matrix[i][j]
                    matrix[i][j] = temp

    for i in range(data.nN-1):
        for pivot in range(i+1, data.nN):
            t = matrix[pivot][i] / matrix[i][i]
            for j in range(data.nN+1):
                matrix[pivot][j] -= t * matrix[i][j]

    for i in range(data.nN-1, -1, -1):
        vector[i] = matrix[i][data.nN]

        for j in range(data.nN):
            if j != i:
                vector[i] -= matrix[i][j] * vector[j]

        vector[i] = vector[i] / matrix[i][i]

    return vector