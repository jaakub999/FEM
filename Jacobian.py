import math
import numpy as np

class Jacobian:
    def __init__(self, npc):
        self.npc = npc
        self.Nksi = np.zeros((npc**2, 4))
        self.Neta = np.zeros((npc**2, 4))
        self.N = np.zeros((npc**2, 4))

        if self.npc == 2:
            self.ksi = [-1/math.sqrt(3), 1/math.sqrt(3), 1/math.sqrt(3), -1/math.sqrt(3)]
            self.eta = [-1/math.sqrt(3), -1/math.sqrt(3), 1/math.sqrt(3), 1/math.sqrt(3)]
            self.ksiW = np.ones(4)
            self.etaW = np.ones(4)

        elif self.npc == 3:
            self.ksi = [-math.sqrt(3/5), 0, math.sqrt(3/5), -math.sqrt(3/5), 0, math.sqrt(3/5), -math.sqrt(3/5), 0, math.sqrt(3/5)]
            self.eta = [-math.sqrt(3/5), -math.sqrt(3/5), -math.sqrt(3/5), 0, 0, 0, math.sqrt(3/5), math.sqrt(3/5), math.sqrt(3/5)]
            self.ksiW = [5/9, 8/9, 5/9, 5/9, 8/9, 5/9, 5/9, 8/9, 5/9]
            self.etaW = [5/9, 5/9, 5/9, 8/9, 8/9, 8/9, 5/9, 5/9, 5/9]

        elif self.npc == 4:
            self.ksi = [-0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136, -0.861136, -0.339981, 0.339981, 0.861136]
            self.eta = [-0.861136, -0.861136, -0.861136, -0.861136, -0.339981, -0.339981, -0.339981, -0.339981, 0.339981, 0.339981, 0.339981, 0.339981, 0.861136, 0.861136, 0.861136, 0.861136]
            self.ksiW = [0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855]
            self.etaW = [0.347855, 0.347855, 0.347855, 0.347855, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.347855, 0.347855, 0.347855, 0.347855]

        for i in range(npc**2):
            self.Nksi[i][0] = -0.25 * (1 - self.eta[i])
            self.Nksi[i][1] = 0.25 * (1 - self.eta[i])
            self.Nksi[i][2] = 0.25 * (1 + self.eta[i])
            self.Nksi[i][3] = -0.25 * (1 + self.eta[i])

        for i in range(npc**2):
            self.Neta[i][0] = -0.25 * (1 - self.ksi[i])
            self.Neta[i][1] = -0.25 * (1 + self.ksi[i])
            self.Neta[i][2] = 0.25 * (1 + self.ksi[i])
            self.Neta[i][3] = 0.25 * (1 - self.ksi[i])

        for i in range(npc**2):
            self.N[i][0] = 0.25 * (1 - self.ksi[i]) * (1 - self.eta[i])
            self.N[i][1] = 0.25 * (1 + self.ksi[i]) * (1 - self.eta[i])
            self.N[i][2] = 0.25 * (1 + self.ksi[i]) * (1 + self.eta[i])
            self.N[i][3] = 0.25 * (1 - self.ksi[i]) * (1 + self.eta[i])

    def calculate_H_and_C(self, element, node, data):
        matrix = [[[0 for k in range(2)] for j in range(2)] for i in range(self.npc**2)]
        reversed_matrix = [[[0 for k in range(2)] for j in range(2)] for i in range(self.npc**2)]
        det = np.zeros(self.npc**2)
        dNdx = np.zeros((self.npc**2, 4))
        dNdy = np.zeros((self.npc**2, 4))

        for i in range(self.npc**2):
            for j in range(4):
                matrix[i][0][0] += self.Nksi[i][j] * node[int(element.ID[j] - 1)].x
                matrix[i][0][1] += self.Neta[i][j] * node[int(element.ID[j] - 1)].x
                matrix[i][1][0] += self.Nksi[i][j] * node[int(element.ID[j] - 1)].y
                matrix[i][1][1] += self.Neta[i][j] * node[int(element.ID[j] - 1)].y

        for i in range(len(matrix)):
            det[i] = (matrix[i][0][0] * matrix[i][1][1]) - (matrix[i][0][1] * matrix[i][1][0])

            for j in range(2):
                reversed_matrix[i][j][j] = matrix[i][1 - j][1 - j] / det[i]
                reversed_matrix[i][j][1 - j] = matrix[i][j][1 - j] / det[i]

        for i in range(len(matrix)):
            for j in range(4):
                dNdx[i][j] = reversed_matrix[i][0][0] * self.Nksi[i][j] + reversed_matrix[i][0][1] * self.Neta[i][j]
                dNdy[i][j] = reversed_matrix[i][1][0] * self.Nksi[i][j] + reversed_matrix[i][1][1] * self.Neta[i][j]

        for i in range(len(matrix)):
            for j in range(4):
                for k in  range(4):
                    HdNdx = dNdx[i][j] * dNdx[i][k]
                    HdNdy = dNdy[i][j] * dNdy[i][k]

                    element.H[j][k] += (HdNdx + HdNdy) * data.k * det[i] * self.ksiW[i] * self.etaW[i]
                    element.C[j][k] += self.N[i][j] * self.N[i][k] * data.ro * data.cp * det[i] * self.ksiW[i] * self.etaW[i]

    def calculate_HBC_and_P(self, element, node, data, side):
        N = np.zeros((self.npc, 4))

        pX = pow(node[int(element.ID[side]) - 1].x - node[int(element.ID[(side + 1) % 4]) - 1].x, 2)
        pY = pow(node[int(element.ID[side]) - 1].y - node[int(element.ID[(side + 1) % 4]) - 1].y, 2)
        det = math.sqrt(pX + pY) / 2

        if self.npc == 2:
            for i in range(self.npc):
                if side == 0:
                    _ksi = self.ksi[i]
                    _eta = -1

                elif side == 1:
                    _ksi = 1
                    _eta = self.eta[i + side]

                elif side == 2:
                    _ksi = self.ksi[i + side]
                    _eta = 1

                elif side == 3:
                    _ksi = -1
                    _eta = self.eta[(i + side) % 4]

                N[i][0] = 0.25 * ((1 - _ksi) * (1 - _eta))
                N[i][1] = 0.25 * ((1 + _ksi) * (1 - _eta))
                N[i][2] = 0.25 * ((1 + _ksi) * (1 + _eta))
                N[i][3] = 0.25 * ((1 - _ksi) * (1 + _eta))

        else:
            for i in range(self.npc):
                if side == 0:
                    _ksi = self.ksi[i]
                    _eta = -1

                elif side == 1:
                    _ksi = 1
                    _eta = self.ksi[i]

                elif side == 2:
                    _ksi = self.ksi[self.npc - i - 1]
                    _eta = 1

                elif side == 3:
                    _ksi = -1
                    _eta = self.ksi[self.npc - i - 1]

                N[i][0] = 0.25 * ((1 - _ksi) * (1 - _eta))
                N[i][1] = 0.25 * ((1 + _ksi) * (1 - _eta))
                N[i][2] = 0.25 * ((1 + _ksi) * (1 + _eta))
                N[i][3] = 0.25 * ((1 - _ksi) * (1 + _eta))

        for i in range(self.npc):
            for j in range(4):
                for k in range(4):
                    element.HBC[j][k] += det * data.alpha * N[i][j] * N[i][k] * self.ksiW[i]
                element.P[j] += -data.alpha * data.t * det * N[i][j] * self.ksiW[i]