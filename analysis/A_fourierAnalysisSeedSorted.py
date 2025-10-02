import ast
import re
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from scipy.fft import fft, fftfreq, ifft
import numpy as np
import pandas as pd
import os
from ZZ_basicUtility import filePrefix, getAllDensityPaths, pathFixer
from reader import TmlReader


def init(fullPath: str):
    files = os.listdir(fullPath)
    files = list(filter(lambda x: x.endswith(".tmlnp") or x.endswith(".tml"), files))
    assert len(files) > 0

    memorySizes = []
    seeds = []
    for file in files:
        parts = file.split("_")
        for part in parts:
            if part.startswith("LC"):
                memorySizes.append(int(part[2:]))
            elif part.startswith("s"):
                try:
                    seeds.append(int(part[1:]))
                except:
                    pass

    memorySizes = sorted(set(memorySizes))
    # print(memorySizes)

    memoryGroupedFileList = [[fullPath + file for file in files if f"_LC{m}_" in file] for m in memorySizes]
    memoryGroupedFileList = list(map(sorted, memoryGroupedFileList))

    # sort groups by memory size
    memoryGroupedFileList = sorted(memoryGroupedFileList, key=lambda memoryGroup: int(list(filter(lambda x: x.startswith("LC"), memoryGroup[0].split("_")))[0][2:]))
    return [seeds, list(zip(*memoryGroupedFileList))]

def plt_formatPlot(fig: plt.Figure,
                   rows: int,
                   columns: int,
                   columnTitles: list[str]) -> None:
    global seeds
    axes = np.reshape(fig.axes, shape=(rows, columns))
    if columns == 5:
        timeEvolutionAxes = axes[:, (0,1)]
        frequencyAxes = axes[:, 2:]
    elif columns == 2 or columns == 3:
        timeEvolutionAxes = axes[:, :1]
        frequencyAxes = axes[:, 1:]

    for col in axes.T:
        ax1 = col[0]
        for ax in col[1:]:
            ax.sharex(ax1)

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

    for ax in frequencyAxes.flatten():
        ax.set_xlim(-0.001,0.02)
        ax.set_xlabel("$f$ / Hz")
        ax.yaxis.set_label_position("right")
        ax.yaxis.grid(which='major', color='gray', linewidth=0.6, alpha=0.8)
        ax.yaxis.grid(which='minor', color='gray', linewidth=0.5, alpha=0.4)
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.yaxis.tick_right()
        ax.grid(True, "both", "y")

    for idx, ax in enumerate(frequencyAxes[:, -1]):
        ax.set_ylabel("magnitude")
        ax.plot([],[], "k<", label=f"s={seeds[idx]}")
        ax.legend(loc="upper right", fontsize="small", framealpha=0.5)

    for idx, ax in enumerate(frequencyAxes[:-1, :].flatten()):
        ax.tick_params(labelbottom=False)

    for ax, title in zip(axes[0], columnTitles):
        ax.set_title(title)

    fig.subplots_adjust(hspace=0, wspace=.27)

def makeFourierPlotForSingleFolder(pathToFolder: str, folder: str):
    global seeds
    fullPath = pathFixer(pathToFolder) + pathFixer(folder) + "timeline/"

    seeds, seedGroupedFileList = init(fullPath)
    seedGroupedFileList = seedGroupedFileList
    columns = 3
    rows = len(seedGroupedFileList)
    fig, axs = plt.subplots(rows, columns, sharex=False, figsize=((7-2)*columns/3 + 2, 9), squeeze=False)
    axs2 = np.zeros(shape=(rows, 5), dtype=object)
    axs2[:, 0] = axs[:, 0]
    axs2[:, 2] = axs[:, 1]
    axs2[:, 4] = axs[:, 2]
    axs = axs2

    folder = re.sub(r"(\w*\d\.([^0]+|0))0*", r"\1", folder.replace("/", "")) # remove trailing 0s
    superTit = re.sub(r"([A-z])(\d)", r"\1 \2", folder) # space between folder name and ID
    fig.suptitle(superTit)
    cutoffA0 = True

    for idx, seedGroup in enumerate(seedGroupedFileList):
        ax = axs[idx, :]
        fourier_datas = []
        fourier_ys = []

        lens = []
        for simulation in seedGroup:
            tml = TmlReader(simulation)
            lens.append(tml.data.size)

        if abs(max(lens) - min(lens)) < 10:
            usePadding = True
        else:
            usePadding = False

        for simulation in seedGroup:
            print(simulation)
            tml = TmlReader(simulation)
            data = tml.data

            if usePadding:
                data = np.pad(data, (max(lens) - data.size, 0), mode='constant')


            metaData = tml.metaData
            T_simulation = metaData["T"]
            # x = np.linspace(0, T_simulation, data.size)
            x = np.linspace(tml.trimmedZeros * metaData["dt_log"], T_simulation, data.size)
            if len(x) == 0:
                continue

            ax[0].plot(x, data)
            if onlyLaterHalf:
                data = data[data.size//2:]

            ### fft ###
            N_samplePoints = len(x)
            dt = T_simulation / data.size
            if "trueDynamic" in simulation:
                if np.trim_zeros(data, "f").size == 0:
                    continue
                fourier_data = fft(np.trim_zeros(data, "f"))
                fourier_data = np.pad(fourier_data, (0, data.size - fourier_data.size), mode='constant', constant_values=0)
            else:
                fourier_data = fft(data)
            fourier_x = fftfreq(N_samplePoints, dt)[:N_samplePoints//2]
            fourier_y = 2/N_samplePoints * np.abs(fourier_data[0:N_samplePoints//2])
            ax[2].plot(fourier_x[int(cutoffA0):], fourier_y[int(cutoffA0):])
            fourier_datas.append(fourier_data)
            fourier_ys.append(fourier_y)

        if len(fourier_datas) == 0:
            continue
        fourier_datas = np.array(fourier_datas)

        fourier_ys = np.array(fourier_ys)
        sum_fourier_ys = np.sum(fourier_ys, axis=0)

        ax[4].plot(fourier_x[int(cutoffA0):], sum_fourier_ys[int(cutoffA0):])


    if "outData/" not in pathToFolder: raise RuntimeError("You are writing from a directory that isn't 'outData'. Adjustments neccissary")
    appendage = pathFixer(pathToFolder.split("outData/")[-1])
    # print(appendage)
    customGraphDataDir = f"graphs_v2/{appendage}fourierTransforms/"
    try: os.makedirs(customGraphDataDir)
    except FileExistsError as e: print(e)

    columnTitles = ["time evolution", "abs$($FT$)$", "$\\sum \\text{abs}\\left(FT\\right)$"]

    plt_formatPlot(fig, rows, columns, columnTitles)

    plt.savefig(f"{customGraphDataDir}{filePrefix(customGraphDataDir)}fourierSeedGrouped_{folder}_c{columns}{"_laterHalf" if onlyLaterHalf else ""}.svg", bbox_inches='tight')
    # plt.show()
    fig.clear()
    plt.close()

try:
    base = pathFixer(sys.argv[1])
except:
    base = pathFixer(r".\outData\staticMem32")

onlyLaterHalf = True

allDensityPaths = [s for s in getAllDensityPaths(base)]
for path in allDensityPaths:
    singles = (path[:-1] if path.endswith("/") else path).split("/")
    pathToFolder = "/".join(singles[:-1])
    folder = singles[-1]
    # print(path)
    # print(pathToFolder, folder)
    makeFourierPlotForSingleFolder(pathToFolder, folder)
