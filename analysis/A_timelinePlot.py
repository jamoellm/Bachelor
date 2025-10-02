import os
import re
import sys
import numpy as np
from ZZ_basicUtility import filePrefix, getAllDensityPaths, pathFixer
from reader import TmlReader
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def init(fullPath: str, MAX_SEEDS:int, tmlFolder: str):
    files = os.listdir(fullPath)
    files = list(filter(lambda x: x.endswith(".tmlnp") or x.endswith(".tml"), files))
    assert len(files) > 0
    LCType = lambda x: int(x)

    memorySizes = []
    for file in files:
        parts = file.split("_")
        for part in parts:
            if part.startswith("LC"):
                memorySizes.append(LCType(part[2:]))
                print(type(LCType(part[2:])))

    memorySizes = sorted(set(memorySizes))
    # memorySizes = [5]
    print(memorySizes)

    files = [f for f in files if any([f"_s{i}_" in f for i in range(0,MAX_SEEDS)])] # don't use more than MAX_SEEDS seeds
    # print(files)
    memoryGroupedFileList = [[fullPath + file for file in files if f"_LC{m}_" in file] for m in memorySizes]
    memoryGroupedFileList = list(map(sorted, memoryGroupedFileList))

    # sort groups by memory size
    memoryGroupedFileList = sorted(memoryGroupedFileList, key=lambda memoryGroup: LCType(list(filter(lambda x: x.startswith("LC"), memoryGroup[0].split("_")))[0][2:]))
    return [memorySizes, memoryGroupedFileList]

def plt_formatPlot(fig: plt.Figure,
                   rows: int,
                   columns: int,
                   columnTitles: list[str],
                   tmlFolder: str) -> None:
    global memorySizes
    axes = np.reshape(fig.axes, shape=(rows, columns))
    timeEvolutionAxes = axes[:, 0]

    for col in axes.T:
        ax1 = col[0]
        for ax in col[1:]:
            ax.sharex(ax1)

    for ax in axes[:-1, :].flatten():
        ax.set_xlabel("")

    for idx, ax in enumerate(timeEvolutionAxes):
        if tmlFolder == "influence":
            ax.set_ylim(-0.1, 2.1)
        else:
            ax.set_ylim(-1.199,1.199)
        ax.set_ylabel("$\\mathcal{P}$")
        ax.set_xlabel("time / s")
        ax.yaxis.grid(which='major', color='black', linewidth=0.6, alpha=0.6)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.grid(True, "both", "y")
        if tmlFolder == "influence":
            ax.plot([], [], "<k", label=f"$\\rho^*$={memorySizes[idx]}")
        else:
            ax.plot([], [], "<k", label=f"M={memorySizes[idx]}")
        ax.legend(loc="upper right")

    for ax, title in zip(axes[0], columnTitles):
        ax.set_title(title)

    for ax in axes.flatten():
        ax.label_outer()
    fig.subplots_adjust(hspace=0)


def doTimelinePlot(pathToFolder: str, folder: str, MAX_SEEDS: int = 99, tmlFolder: str = "timeline") -> None:
    global memorySizes, metaData, memoryGroup
    pathToFolder, folder = list(map(pathFixer, [pathToFolder, folder]))

    fullPath = pathFixer(pathToFolder) + pathFixer(folder) + pathFixer(tmlFolder)
    print(fullPath)
    if not os.path.isdir(fullPath): return
    memorySizes, memoryGroupedFileList = init(fullPath, MAX_SEEDS, tmlFolder)

    columns = 1
    rows = len(memoryGroupedFileList)
    fig, axs = plt.subplots(rows, columns, sharex=False, figsize=(7,9), squeeze=False)
    folder = re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", "")) # remove trailing 0s
    superTit = re.sub(r"([A-z])(\d)", r"\1 \2", folder) # space between folder name and ID
    fig.suptitle(superTit)

    for idx, memoryGroup in enumerate(memoryGroupedFileList):
        ax = axs[idx, :]
        for jdx, simulation in enumerate(memoryGroup):
            print(simulation)
            tml = TmlReader(simulation)
            data = tml.data
            if tmlFolder == "timeline":
                metaData = tml.metaData
                T_simulation = metaData["T"]
                # x = np.linspace(0, T_simulation, data.size)
                x = np.linspace(tml.trimmedZeros * metaData["dt_log"], T_simulation, data.size)
                ax[0].plot(x, data)
            elif tmlFolder == "influence":
                ax[0].plot(data) # , ls="", marker=".", markersize=1



    if "outData/" not in pathToFolder: raise RuntimeError("You are writing from a directory that isn't 'outData'. Adjustments neccissary")
    appendage = pathFixer(pathToFolder.split("outData/")[-1])
    print(appendage)
    customGraphDataDir = f"graphs_v2/{appendage}{tmlFolder}Plots/"
    try: os.makedirs(customGraphDataDir)
    except FileExistsError as e: print(e)

    plt_formatPlot(fig, rows, columns, ["time evolution", "density",  "reduced density"], tmlFolder)
    plt.savefig(f"{customGraphDataDir}{filePrefix(customGraphDataDir)}{tmlFolder}_{folder.removesuffix("/")}_c{columns}_s{len(memoryGroup)}.svg", bbox_inches='tight')
    # plt.show()
    fig.clear()
    plt.close()


try:
    base = pathFixer(sys.argv[1])
except:
    base = pathFixer(r".\outData\staticMem32")

tmlFolder = "timeline"

base= pathFixer(base)
allDensityPaths = [s for s in getAllDensityPaths(base)]
for path in allDensityPaths:
    singles = (path[:-1] if path.endswith("/") else path).split("/")
    pathToFolder = "/".join(singles[:-1])
    folder = singles[-1]

    doTimelinePlot(pathToFolder, folder, tmlFolder=tmlFolder)

