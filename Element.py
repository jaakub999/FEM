import numpy as np

class Element:
    def __init__(self):
        self.ID = np.zeros(4)
        self.H = np.zeros((4, 4))
        self.C = np.zeros((4, 4))
        self.HBC = np.zeros((4, 4))
        self.P = np.zeros(4)