import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from ZZ_basicUtility import filePrefix, getAllDensityPaths, pathFixer, errstyle, titleTranslator
from reader import TmlReader
from C_polarizability import plt_formatPlot as C_plt_formatPlot
from C_polarizability import doErrorbarStuff as C_doErrorbarStuff
from C_polarizability import init as C_init
from C_polarizability import plt_formatSecondaryPlot

totalAbsMean = [0, 0]
totalMeanAbs = [0, 0]
totalAbsMean_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]
totalMeanAbs_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]

def init(fullPath: str):
    files = os.listdir(fullPath)
    files = list(filter(lambda x: x.endswith(".tmlnp") or x.endswith(".tml"), files))
    assert len(files) > 0
    return C_init(files, "LC")

def plt_formatPlot(fig: plt.Figure,
                   SEEDS: int,
                   rows: int,
                   columns: int,
                   columnTitles: list[str]) -> None:
    global memorySizes
    C_plt_formatPlot(fig, rows, columns, SEEDS, columnTitles, memorySizes, "LC")

def doErrorbarStuff(ax, jdx, data, errstyle, idxx):
    global totalAbsMean, totalMeanAbs, totalAbsMean_yerr, totalMeanAbs
    C_doErrorbarStuff(ax, jdx, data, errstyle,
                      idxx,
                      totalAbsMean, totalMeanAbs,
                      totalAbsMean_yerr, totalMeanAbs_yerr)

def doPolarizability(pathToFolder: str, folder: str, SEEDS: int):
    global memorySizes, errstyle, totalAbsMean, totalMeanAbs, totalMeanAbs_yerr, totalAbsMean_yerr
    fullPath = pathFixer(pathToFolder) + pathFixer(folder) + "timeline/"
    memorySizes, memoryGroupedFileList = init(fullPath)

    columns = 3
    rows = len(memoryGroupedFileList)
    fig, axs = plt.subplots(rows, columns, figsize=(7,9), sharex=False, squeeze=False)
    superTit = re.sub(r"([A-z])(\d)", r"\1 \2",
                        re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", "")) # remove trailing 0s
                      ) # space between folder name and ID
    fig.suptitle(superTit)

    fig2, axs2 = plt.subplots(1, 2, figsize=(7,3), squeeze=True)
    fig2.suptitle(superTit)

    for idx, memoryGroup in enumerate(memoryGroupedFileList):
        ax = axs[idx, :]
        totalAbsMean = [0, 0]
        totalMeanAbs = [0, 0]
        totalAbsMean_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]
        totalMeanAbs_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]

        for jdx, simulation in enumerate(memoryGroup):
            print(simulation)
            tml = TmlReader(pathFixer(pathToFolder) + pathFixer(folder) + "timeline/" + simulation)
            data = tml.data
            metaData = tml.metaData
            T_simulation = metaData["T"]
            # x = np.linspace(0, T_simulation, data.size)
            x = np.linspace(tml.trimmedZeros * metaData["dt_log"], T_simulation, data.size)
            ax[0].plot(x, data)

            # mean & std
            doErrorbarStuff(ax[1], jdx, data, errstyle, 0)

            # reduced mean & std
            data = data[len(data)//2:]
            doErrorbarStuff(ax[2], jdx, data, errstyle, 1)

        # print(len(memoryGroup) - 1/8, [totalAbsMean[0]], "<--")
        errstyle3 = errstyle.copy()
        errstyle3["fmt"] = "k"

        # abs(mean(data))
        ax[1].errorbar(len(memoryGroup) - 1/8, [totalAbsMean[0] / len(memoryGroup)], totalAbsMean_yerr[0] / len(memoryGroup), marker=".", **errstyle3)
        axs2[0].errorbar([memorySizes[idx]-1/4]  , [totalAbsMean[0] / len(memoryGroup)], totalAbsMean_yerr[0] / len(memoryGroup), marker=".", **errstyle3)
        # mean(abs(data))
        ax[1].errorbar(len(memoryGroup) + 1/8, [totalMeanAbs[0] / len(memoryGroup)], totalMeanAbs_yerr[0] / len(memoryGroup), marker="x", **errstyle3)
        axs2[0].errorbar([memorySizes[idx]+1/4]  , [totalMeanAbs[0] / len(memoryGroup)], totalMeanAbs_yerr[0] / len(memoryGroup), marker="x", **errstyle3)

        # reduced abs(mean(data))
        ax[2].errorbar(len(memoryGroup) - 1/8, [totalAbsMean[1] / len(memoryGroup)], totalAbsMean_yerr[1] / len(memoryGroup), marker=".", **errstyle3)
        axs2[1].errorbar([memorySizes[idx]-1/4]  , [totalAbsMean[1] / len(memoryGroup)], totalAbsMean_yerr[1] / len(memoryGroup), marker=".", **errstyle3)
        # reduced mean(abs(data))
        ax[2].errorbar(len(memoryGroup) + 1/8, [totalMeanAbs[1] / len(memoryGroup)], totalMeanAbs_yerr[1] / len(memoryGroup), marker="x", **errstyle3)
        axs2[1].errorbar([memorySizes[idx]+1/4]  , [totalMeanAbs[1] / len(memoryGroup)], totalMeanAbs_yerr[1] / len(memoryGroup), marker="x", **errstyle3)

    if "outData/" not in pathToFolder: raise RuntimeError("You are writing from a directory that isn't 'outData'. Adjustments neccissary")
    appendage = pathFixer(pathToFolder.split("outData/")[-1])
    print(appendage)
    customGraphDataDir = f"graphs_v2/{appendage}polarizability/"
    try: os.makedirs(customGraphDataDir)
    except FileExistsError as e: print(e)

    customGraphDataDir2 = f"graphs_v2/{appendage}polarizabilityScatter/"
    try: os.makedirs(customGraphDataDir2)
    except FileExistsError as e: print(e)

    plt_formatSecondaryPlot(fig2, "LC")

    plt_formatPlot(fig, SEEDS, rows, columns, ["time evolution", "mean & uncertainty",  "reduced\nmean & uncertainty"])
    fig.savefig(f"{customGraphDataDir}{filePrefix(customGraphDataDir)}polarizability_{re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", ""))}_c{columns}.svg", bbox_inches='tight')
    fig2.savefig(f"{customGraphDataDir2}{filePrefix(customGraphDataDir)}polarizabilityScatter_{re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", ""))}_c{columns}.svg", bbox_inches='tight')

    # plt.show()


    fig.clear()
    fig2.clear()
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
    # print(path)
    # print(pathToFolder, folder)
    doPolarizability(pathToFolder, folder, 8)
