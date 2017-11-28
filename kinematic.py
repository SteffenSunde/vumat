def plot_results(file):
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.genfromtxt(file, delimiter=",")
    x_label = data[0,0]
    y_label = data[0,1]
    data = data[1:, :]

    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel(str(x_label))
    plt.ylabel(str(y_label))

    plt.show()

    return


def plotty():
    import numpy as np
    import matplotlib.pyplot as plt

    file = "debug.txt"

    data = np.genfromtxt(file, delimiter=",")
    x_label = data[0,0]
    y_label = data[0,1]
    data = data[1:, :]

    plt.plot(data[-3:, 0], data[-30:, 1])
    plt.plot(data[-30:-27, :])
    plt.xlabel(str(x_label))
    plt.ylabel(str(y_label))

    plt.show()

    return
