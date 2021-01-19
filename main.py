from GlobalData import GlobalData
from Simulation import SOE

if __name__ == "__main__":
    print("Time[s]", end = "")
    print("{:>15}".format("MinTemp"), end = "")
    print("{:>15}".format("MaxTemp"))

    node = []
    element = []

    data = GlobalData(element, node)
    Simulation = SOE(data, element, node)