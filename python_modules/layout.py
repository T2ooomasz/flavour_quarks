import numpy as np
import matplotlib.pyplot as plt

#
# Layout
#
def initialize_layout(site_elements, a):
    Layout = {"2Dneighbours": np.zeros((site_elements, site_elements), dtype=int),
              "realLocation": [[[]]]}
    k = 1
    for i in range(site_elements):
        for j in range(site_elements):
            Layout["2Dneighbours"][j][i] = k
            Layout["realLocation"].append([k,[a*i,a*j]])
            k += 1
    Layout["realLocation"] = Layout["realLocation"][1:]

    return Layout

def plot_layout(Layout):
    X = np.zeros(len(Layout["realLocation"]))
    Y = np.zeros(len(Layout["realLocation"]))
    for i in range(len(Layout["realLocation"])):
        elem = Layout["realLocation"][i]
        X[i] = elem[1][0]
        Y[i] = elem[1][1]
    plt.scatter(X, Y, c='blue', linewidths=1)
    plt.gca().set_aspect('equal', 'box')
    plt.title("grid")
    plt.show()   
