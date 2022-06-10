import argparse
from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams


def main():
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['IPAexGothic']
    rcParams['pdf.fonttype'] = 42  # TrueType
    rcParams['font.size'] = 9  # footnotesize
    markers = cycle(["o", "v", "D", "s", "*", "^", "d", "x"])


    parser = argparse.ArgumentParser()
    parser.add_argument("log_file")
    parser.add_argument("-o", "--output", default=argparse.SUPPRESS)
    args = parser.parse_args()
    log_file = args.log_file


    labels = ["fixed step", "Nesterov"]
    data = np.loadtxt(log_file)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r"iteration $k$")
    ax.set_ylabel(r"誤差 $f(w_k) - f(w^\star)$")
    ax.set_yscale("log")

    for log, label, marker in zip(data[:, 1:].T, labels, markers):
        ax.plot(data[:, 0], log, marker=marker, label=label)

    ax.legend()

    if "output" in args:
        plt.savefig(args.output, bbox_inches='tight')
    else:
        plt.show()


if __name__ == "__main__":
    main()
