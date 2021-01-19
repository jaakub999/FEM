from Element import Element
from Node import Node

class GlobalData:
    def __init__(self, element, node):
        file = open("data", "r")
        self.W = float(file.readline())
        self.H = float(file.readline())
        self.nW = int(file.readline())
        self.nH = int(file.readline())
        self.k = int(file.readline())
        self.cp = int(file.readline())
        self.ro = int(file.readline())
        self.alpha = int(file.readline())
        self.t = int(file.readline())
        self.t0 = int(file.readline())
        self.simulationTime = int(file.readline())
        self.stepTime = int(file.readline())
        self.npc = int(file.readline())
        file.close()

        self.nN = self.nH * self.nW
        self.nE = (self.nH - 1) * (self.nW - 1)

        j = 1
        for i in range(self.nE):
            element.append(Element())
            element[i].ID[0] = j
            element[i].ID[1] = j + self.nH
            element[i].ID[2] = j + self.nH + 1
            element[i].ID[3] = j + 1
            j += 1

            if (i + 1) % (self.nH - 1) == 0:
                j += 1

        deltaX = self.W / (self.nW - 1)
        deltaY = self.H / (self.nH - 1)

        for i in range(self.nN):
            x = int(i / self.nH) * deltaX
            y = int(i % self.nH) * deltaY

            if int(i / self.nH) == 0 or i % self.nH == 0 or int(i / self.nH) == self.nW - 1 or i % self.nH == self.nH - 1:
                BC = 1

            else:
                BC = 0

            node.append(Node(x, y, BC, self.t0))