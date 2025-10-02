from collections import Counter
import os
import re
import sys
import numpy as np
from ZZ_basicUtility import filePrefix, getAllDensityPaths, pathFixer
from reader import ParquetChunkReader, TmlReader
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def init(fullPath: str, MAX_SEEDS:int):
    files = os.listdir(fullPath)
    files = list(filter(lambda x: x.endswith(".tmlnp") or x.endswith(".tml"), files))
    assert len(files) > 0

    memorySizes = []
    for file in files:
        parts = file.split("_")
        for part in parts:
            if part.startswith("LC"):
                memorySizes.append(int(part[2:]))

    memorySizes = sorted(set(memorySizes))
    files = [f for f in files if any([f"_s{i}_" in f for i in range(0, MAX_SEEDS)])] # don't use more than MAX_SEEDS seeds

    memoryGroupedFileList = [[fullPath + file for file in files if f"_LC{m}_" in file] for m in memorySizes]
    memoryGroupedFileList = list(map(sorted, memoryGroupedFileList))

    # sort groups by memory size
    memoryGroupedFileList = sorted(memoryGroupedFileList, key=lambda memoryGroup: int(list(filter(lambda x: x.startswith("LC"), memoryGroup[0].split("_")))[0][2:]))
    return [memorySizes, memoryGroupedFileList]

def plt_formatPlot(fig: plt.Figure,
                   rows: int,
                   columns: int,
                   columnTitles: list[str]) -> None:
    global memorySizes
    axes = np.reshape(fig.axes, shape=(rows, columns))
    timeEvolutionAxes = axes[:, 0]
    densityAxes = axes[:, 1:3]

    for col in axes.T:
        ax1 = col[0]
        for ax in col[1:]:
            ax.sharex(ax1)

    for ax in axes[:-1, :].flatten():
        ax.set_xlabel("")

    for ax in timeEvolutionAxes:
        ax.set_ylim(-1.199,1.199)
        ax.set_ylabel("$\\mathcal{P}$")
        ax.set_xlabel("time / s")
        ax.yaxis.grid(which='major', color='black', linewidth=0.6, alpha=0.6)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.grid(True, "both", "y")

    for ax in densityAxes.flatten():
        ax.set_ylim(0,0.999)
        ax.set_xlim(-1.1,1.1)
        ax.set_ylabel("$\\mathcal{P}$")
        ax.yaxis.set_label_position("right")
        ax.yaxis.grid(which='major', color='gray', linewidth=0.6, alpha=0.8)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.yaxis.tick_right()
        ax.grid(True, "both", "y")

    for idx, ax in enumerate(densityAxes[:, -1]):
        ax.set_ylabel("$p$")
        ax.plot([],[], "k<", label=f"M={memorySizes[idx]}")
        ax.legend(loc="upper right")

    for idx, ax in enumerate(densityAxes[:-1, :].flatten()):
        ax.tick_params(labelbottom=False)

    for ax in densityAxes.flatten():
        ax.set_xlabel("$\\mathcal{P}$")

    for ax, title in zip(axes[0], columnTitles):
        ax.set_title(title)

    fig.subplots_adjust(hspace=0)

def doHist(ax, jdx, data):
    global metaData, memoryGroup
    counts = Counter(data)
    bins = np.array(list(counts.keys()))
    values = np.array(list(counts.values()))
    area = np.sum(values)
    dx = 2 / metaData["NN"]
    dx *= 0.8
    ax.bar(bins + dx * jdx / len(memoryGroup), np.array(values) / area, width=dx/len(memoryGroup))

def do2Dplots(pathToFolder: str, folder: str, MAX_SEEDS:int =  99) -> None:
    global memorySizes, metaData, memoryGroup

    fullPath = pathFixer(pathToFolder) + pathFixer(folder) + "timeline/"
    memorySizes, memoryGroupedFileList = init(fullPath, MAX_SEEDS)

    columns = 3
    rows = len(memoryGroupedFileList)
    fig, axs = plt.subplots(rows, columns, sharex=False, figsize=(7,7), squeeze=False, dpi=400)
    folder = re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", "")) # remove trailing 0s
    superTit = re.sub(r"([A-z])(\d)", r"\1 \2", folder) # space between folder name and ID
    fig.suptitle(superTit)

    for idx, memoryGroup in enumerate(memoryGroupedFileList):
        ax = axs[idx, :]
        for jdx, simulation in enumerate(memoryGroup):
            tml = TmlReader(simulation)
            data = tml.data
            metaData = tml.metaData
            T_simulation = metaData["T"]
            x = np.linspace(0, T_simulation, data.size)
            ax[0].plot(x, data)

            doHist(ax[1], jdx, data)
            doHist(ax[2], jdx, data[len(data)//2:])


    if "/!outData/" not in pathToFolder: raise RuntimeError("You are writing from a directory that isn't 'outData'. Adjustments neccissary")
    appendage = pathFixer(pathToFolder.split("outData/")[-1])
    # print(appendage)
    customGraphDataDir = pathFixer(f"graphs_v2/{appendage}histograms")
    try: os.makedirs(customGraphDataDir)
    except FileExistsError as e: print(e)

    plt_formatPlot(fig, rows, columns, ["time evolution", "density",  "reduced density"])
    plt.savefig(f"{customGraphDataDir}{filePrefix(customGraphDataDir)}histograms_{folder}_c{columns}_s{MAX_SEEDS}.png", bbox_inches='tight')
    # plt.show()
    fig.clear()
    plt.close()

try:
    base = pathFixer(sys.argv[1])
except:
    base = pathFixer(r".\outData\staticMem32")

allDensityPaths = [s for s in getAllDensityPaths(base)]
for path in allDensityPaths:
    singles = (path[:-1] if path.endswith("/") else path).split("/")
    pathToFolder = "/".join(singles[:-1])
    folder = singles[-1]
    do2Dplots(pathToFolder, folder)




