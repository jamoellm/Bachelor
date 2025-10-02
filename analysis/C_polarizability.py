import os
import sys
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
import matplotlib.pyplot as plt
from ZZ_basicUtility import filePrefix, getAllDensityPaths, labelTranslator, pathFixer, errstyle, titleTranslator
from reader import TmlReader
import re
import pandas as pd

totalAbsMean = [0, 0]
totalMeanAbs = [0, 0]
totalAbsMean_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]
totalMeanAbs_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]
allAbsMeanPolarizabilities = []
allAbsMeanPolarizabilities_reduced = []
allMeanAbsPolarizabilities = []
allMeanAbsPolarizabilities_reduced = []

def plt_formatPlot(fig: plt.Figure,
                   rows: int,
                   columns: int,
                   SEEDS: int,
                   columnTitles: list[str],
                   groupIDs: list[float | int],
                   subGroupKey: str) -> None:
    axes = np.reshape(fig.axes, shape=(rows, columns))
    timeEvolutionAxes = axes[:, :1]
    polAxes = axes[:, 1:3]

    # share x axes in each column
    for col in axes.T:
        ax1 = col[0]
        for ax in col[1:]:
            ax.sharex(ax1)

    # place xlabels only at the bottom
    for ax in axes[:-1, :].flatten():
        ax.set_xlabel("")

    for ax in timeEvolutionAxes[:, 0]:
        ax.set_ylabel("$\\mathcal{P}$")

    for ax in timeEvolutionAxes.flatten():
        ax.set_ylim(-1.15,1.15)
        ax.set_xlabel("time / s")
        ax.yaxis.grid(which='major', color='black', linewidth=0.6, alpha=0.6)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.grid(True, "both", "y")

    for idx, ax in enumerate(timeEvolutionAxes.flatten()):
        ax.plot([], [], ">k", label=f"{titleTranslator(subGroupKey)}={groupIDs[idx]}")
        ax.legend(loc="upper left", fontsize="small", framealpha=0.5)

    for ax in polAxes.flatten():
        ax.set_ylim(-0.1,1.1)
        ax.set_xlabel("seed")
        ax.yaxis.set_label_position("right")
        ax.yaxis.grid(which='major', color='black', linewidth=0.6, alpha=0.6)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.yaxis.tick_right()
        ax.grid(True, "both", "y")
        ax.set_ylabel("$\\mathscr{P}$ & $\\widetilde\\mathscr{P}$")

    for ax in polAxes[-1:, :].flatten():
        ax.set_xticks(range(SEEDS+1))
        ax.set_xticklabels([f"{i}" for i in range(SEEDS)] + ["avg"])

    for ax, title in zip(axes[0], columnTitles):
        ax.set_title(title)

    for ax in axes.flatten():
        ax.label_outer()

    fig.subplots_adjust(hspace=0)

def plt_formatSecondaryPlot(fig2, innerKey: str):
    axs2 = fig2.axes
    axs2[0].set_title("average abs(mean) & mean(abs)")
    axs2[1].set_title("reduced\naverage abs(mean) & mean(abs)")
    for ax in axs2:
        ax.set_ylim(-0.1, 1.1)
        for value in [0, 0.5, 1]:
            ax.axhline(value, color="gray", ls="--")
        ax.set_xlabel(labelTranslator(innerKey))
    axs2[0].set_ylabel("(absolute) polarizability")
    fig2.tight_layout()

def doErrorbarStuff(ax, jdx, data, errstyle, idxx, totalAbsMean, totalMeanAbs, totalAbsMean_yerr, totalMeanAbs_yerr):
    mean = np.mean(data)

    stddev = np.std(data)
    yerr = [[stddev], [stddev]]

    ax.errorbar([jdx-1/8], [abs(mean)], yerr, **errstyle, color=f"C{jdx}")
    totalAbsMean[idxx] += abs(mean)
    totalAbsMean_yerr[idxx] += np.array(yerr)

    absData = np.abs(data)
    meanAbs = np.mean(absData)
    stddev = np.std(data)
    yerr = [[stddev], [stddev]]
    errstyle2 = errstyle.copy()
    errstyle2["fmt"] = "x"
    ax.errorbar([jdx+1/8], [meanAbs], yerr, **errstyle2, color=f"C{jdx}")
    totalMeanAbs[idxx] += meanAbs
    totalMeanAbs_yerr[idxx] += np.array(yerr)

def doOuterGrouping(allFiles: list[str], groupingKeyword: str, sortingKeyword: str):
    """groups the files by figure"""
    groupIDs = []
    for file in allFiles:
        parts = file.split("_")
        for part in parts:
            if part.startswith(groupingKeyword):
                groupIDs.append(float(part[len(groupingKeyword):]))

    groupIDs = sorted(list(set(groupIDs)))
    print(groupIDs)
    if len(groupIDs) == 0:
        raise ValueError("No groups for `groupingKeyword` were found. Check for typos.")

    if all(map(lambda x: x == int(x), groupIDs)):
        groupIDs = list(map(int, groupIDs))


    groupedFileList = [[file for file in allFiles if bool(re.search(f"_{groupingKeyword}{str(g).replace(".", r"\.")}0*_", file)) ] for g in groupIDs]
    sortKey = lambda file: float(list(filter(lambda part: part.startswith(sortingKeyword), file.split("_")))[0][len(sortingKeyword):])

    groupedFileList = list(map(lambda x: sorted(x, key=sortKey), groupedFileList))
    return [groupIDs, groupedFileList]

def init(files: list[str], groupingKeyword: str) -> tuple[list[int | float], list[list[str]]]:
    groupIDs = []
    for file in files:
        parts = file.split("_")
        for part in parts:
            if part.startswith(groupingKeyword):
                groupIDs.append(float(part[len(groupingKeyword):]))

    groupIDs = sorted(list(set(groupIDs)))

    if all(map(lambda x: x == int(x), groupIDs)):
        groupIDs = list(map(int, groupIDs))

    groupedFileList = [[file for file in files if bool(re.search(f"_{groupingKeyword}{str(g).replace(".", r"\.")}0*_", file)) ] for g in groupIDs]
    return [groupIDs, groupedFileList]

def doPolarizability(outerID: float, fileList: list[str], outerKey: str, innerKey: str):
    """takes a list of all files that will be plotted in a single figure"""
    global allAbsMeanPolarizabilities, allAbsMeanPolarizabilities_reduced, allMeanAbsPolarizabilities, allMeanAbsPolarizabilities_reduced

    innerIDs, groupedFileList = init(fileList, innerKey)
    pathToFolder = pathFixer(fileList[0].split("/density")[0])
    folder = outerKey + str(np.round(outerID, 2))
    SEEDS = len(groupedFileList[0])
    print(SEEDS)

    columns = 3
    rows = len(innerIDs)
    fig, axs = plt.subplots(rows, columns, figsize=(7,9), sharex=False, squeeze=False)
    fig2, axs2 = plt.subplots(1, 2, figsize=(7,3), squeeze=True)

    graphTitle = titleTranslator(outerKey) + " " + str(np.round(outerID, 2))
    fig.suptitle(graphTitle)
    fig2.suptitle(graphTitle)

    minDx = 999
    for i in range(1, len(innerIDs)):
        dx = innerIDs[i] - innerIDs[i-1]
        if dx < minDx: minDx = dx

    for idx, innerGroup in enumerate(groupedFileList):
        ax = axs[idx, :]
        totalAbsMean = [0, 0]
        totalMeanAbs = [0, 0]
        totalAbsMean_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]
        totalMeanAbs_yerr = [np.array([[0], [0]], dtype=float), np.array([[0], [0]], dtype=float)]

        for jdx, simulation in enumerate(innerGroup):
            print(simulation)
            tml = TmlReader(simulation)
            data = tml.data
            metaData = tml.metaData
            T_simulation = metaData["T"]
            # x = np.linspace(0, T_simulation, data.size)
            x = np.linspace(tml.trimmedZeros * metaData["dt_log"], T_simulation, data.size)
            ax[0].plot(x, data)

            # mean & stddev
            doErrorbarStuff(ax[1], jdx, data, errstyle, 0, totalAbsMean, totalMeanAbs, totalAbsMean_yerr, totalMeanAbs_yerr)

            # reduced mean & stddev
            data = data[len(data)//2:]
            doErrorbarStuff(ax[2], jdx, data, errstyle, 1, totalAbsMean, totalMeanAbs, totalAbsMean_yerr, totalMeanAbs_yerr)

        # print(len(innerGroup) - 1/8, [totalAbsMean[0]], "<--")
        errstyle3 = errstyle.copy()
        errstyle3["fmt"] = "k"

        # abs(mean(data))
        ax[1].errorbar(len(innerGroup) - 1/8, [totalAbsMean[0] / len(innerGroup)], totalAbsMean_yerr[0] / len(innerGroup), marker=".", **errstyle3)
        axs2[0].errorbar([innerIDs[idx]-1/8 * minDx]  , [totalAbsMean[0] / len(innerGroup)], totalAbsMean_yerr[0] / len(innerGroup), marker=".", **errstyle3)
        # mean(abs(data))
        ax[1].errorbar(len(innerGroup) + 1/8, [totalMeanAbs[0] / len(innerGroup)], totalMeanAbs_yerr[0] / len(innerGroup), marker="x", **errstyle3)
        axs2[0].errorbar([innerIDs[idx]+1/8 * minDx]  , [totalMeanAbs[0] / len(innerGroup)], totalMeanAbs_yerr[0] / len(innerGroup), marker="x", **errstyle3)
        allAbsMeanPolarizabilities.append(np.array([outerID,
                                                innerIDs[idx],
                                                totalAbsMean[0] / len(innerGroup),
                                                totalAbsMean_yerr[0][0][0] / len(innerGroup)]))
        allMeanAbsPolarizabilities.append(np.array([outerID,
                                                innerIDs[idx],
                                                totalMeanAbs[0] / len(innerGroup),
                                                totalMeanAbs_yerr[0][0][0] / len(innerGroup)]))

        # reduced abs(mean(data))
        ax[2].errorbar(len(innerGroup) - 1/8, [totalAbsMean[1] / len(innerGroup)], totalAbsMean_yerr[1] / len(innerGroup), marker=".", **errstyle3)
        axs2[1].errorbar([innerIDs[idx]-1/8 * minDx]  , [totalAbsMean[1] / len(innerGroup)], totalAbsMean_yerr[1] / len(innerGroup), marker=".", **errstyle3)
        # reduced mean(abs(data))
        ax[2].errorbar(len(innerGroup) + 1/8, [totalMeanAbs[1] / len(innerGroup)], totalMeanAbs_yerr[1] / len(innerGroup), marker="x", **errstyle3)
        axs2[1].errorbar([innerIDs[idx]+1/8 * minDx]  , [totalMeanAbs[1] / len(innerGroup)], totalMeanAbs_yerr[1] / len(innerGroup), marker="x", **errstyle3)
        allAbsMeanPolarizabilities_reduced.append(np.array([outerID,
                                                        innerIDs[idx],
                                                        totalAbsMean[1] / len(innerGroup),
                                                        totalAbsMean_yerr[1][0][0] / len(innerGroup)]))
        allMeanAbsPolarizabilities_reduced.append(np.array([outerID,
                                                        innerIDs[idx],
                                                        totalMeanAbs[1] / len(innerGroup),
                                                        totalMeanAbs_yerr[1][0][0] / len(innerGroup)]))

    if "outData/" not in pathToFolder: raise RuntimeError("You are writing from a directory that isn't 'outData'. Adjustments neccissary.")
    appendage = pathFixer(pathToFolder.split("outData/")[-1])
    print(appendage)
    customGraphDataDir = f"graphs_v2/{appendage}polarizability/"
    try: os.makedirs(customGraphDataDir)
    except FileExistsError as e: print(e)

    customGraphDataDir2 = f"graphs_v2/{appendage}polarizabilityScatter/"
    try: os.makedirs(customGraphDataDir2)
    except FileExistsError as e: print(e)

    plt_formatSecondaryPlot(fig2, innerKey)

    titles = ["time evolution", "mean & uncertainty",  "reduced\nmean & uncertainty"]
    plt_formatPlot(fig, rows, columns, SEEDS, titles, innerIDs, innerKey)
    fig.savefig(f"{customGraphDataDir}{filePrefix(customGraphDataDir)}polarizability_{folder}_c{columns}.svg", bbox_inches='tight')
    fig2.savefig(f"{customGraphDataDir2}{filePrefix(customGraphDataDir)}polarizabilityScatter_{folder}_c{columns}.svg", bbox_inches='tight')

    # plt.show()

    fig.clear()
    fig2.clear()


if __name__ == "__main__":
    try:
        base = pathFixer(sys.argv[1])
        outerGrouping = sys.argv[2]
        innerGrouping = sys.argv[3]
    except:
        base = pathFixer(r".\outData\staticMem")
        outerGrouping = "LC"
        innerGrouping = "density"


    allDensityPaths = [s for s in getAllDensityPaths(base) if "density0.3" not in s]
    print(allDensityPaths)
    allFiles = []
    for path in allDensityPaths:
        allFiles.extend([path + "timeline/" + file for file in os.listdir(path + "timeline/")])

    groupIDs, groupedFileList = doOuterGrouping(allFiles, outerGrouping, innerGrouping)
    for id, fileList in zip(groupIDs, groupedFileList):
        doPolarizability(id, fileList, outerGrouping, innerGrouping)


    # see `statiSim1_p2.py`
    dataPath = pathFixer(f"data/polarizability/{base[:-1].split("/")[-1]}")
    try: os.makedirs(dataPath)
    except FileExistsError as e: print(e)

    allAbsMeanPolarizabilities = np.array(allAbsMeanPolarizabilities)
    allAbsMeanPolarizabilities_reduced = np.array(allAbsMeanPolarizabilities_reduced)

    allAbsMeanPolarizabilitiesPD = pd.DataFrame({
        outerGrouping: allAbsMeanPolarizabilities[:, 0],
        innerGrouping: allAbsMeanPolarizabilities[:, 1],
        "polarizability": allAbsMeanPolarizabilities[:, 2],
        "dpolarizability": allAbsMeanPolarizabilities[:, 3]})
    allAbsMeanPolarizabilitiesPD_reduced = pd.DataFrame({
        outerGrouping: allAbsMeanPolarizabilities_reduced[:, 0],
        innerGrouping: allAbsMeanPolarizabilities_reduced[:, 1],
        "polarizability": allAbsMeanPolarizabilities_reduced[:, 2],
        "dpolarizability": allAbsMeanPolarizabilities_reduced[:, 3]})

    allAbsMeanPolarizabilitiesPD.to_csv(dataPath + "allAbsMeanPolarizabilities.txt", index=False)
    allAbsMeanPolarizabilitiesPD_reduced.to_csv(dataPath + "allAbsMeanPolarizabilities_reduced.txt", index=False)




    allMeanAbsPolarizabilities = np.array(allMeanAbsPolarizabilities)
    allMeanAbsPolarizabilities_reduced = np.array(allMeanAbsPolarizabilities_reduced)

    allMeanAbsPolarizabilitiesPD = pd.DataFrame({
        outerGrouping: allMeanAbsPolarizabilities[:, 0],
        innerGrouping: allMeanAbsPolarizabilities[:, 1],
        "polarizability": allMeanAbsPolarizabilities[:, 2],
        "dpolarizability": allMeanAbsPolarizabilities[:, 3]})
    allMeanAbsPolarizabilitiesPD_reduced = pd.DataFrame({
        outerGrouping: allMeanAbsPolarizabilities_reduced[:, 0],
        innerGrouping: allMeanAbsPolarizabilities_reduced[:, 1],
        "polarizability": allMeanAbsPolarizabilities_reduced[:, 2],
        "dpolarizability": allMeanAbsPolarizabilities_reduced[:, 3]})

    allMeanAbsPolarizabilitiesPD.to_csv(dataPath + "allMeanAbsPolarizabilities.txt", index=False)
    allMeanAbsPolarizabilitiesPD_reduced.to_csv(dataPath + "allMeanAbsPolarizabilities_reduced.txt", index=False)
