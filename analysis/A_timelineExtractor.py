import os
import numpy as np
import matplotlib.pyplot as plt
from ZZ_basicUtility import getAllDensityPaths, pathFixer
from reader import DatChunkReader, ParquetChunkReader

FILE_EXTENSION = ".parquet"

def init(pathToFolder: str, folderName: str):
    global FILE_EXTENSION
    basePath = pathFixer(pathToFolder + folderName)

    allFiles = os.listdir(basePath)
    allFiles = list(filter(lambda x: not x.startswith("INCOMPLETE_"), allFiles))
    if len(list(filter(lambda x: x.endswith(FILE_EXTENSION), allFiles))) == 0: FILE_EXTENSION = ".dat"
    allFiles = list(filter(lambda x: x.endswith(FILE_EXTENSION), allFiles))


    try: assert len(allFiles) > 0
    except:
        print(basePath)
        print(*allFiles, sep="\n")
        assert len(allFiles) > 0

    memorySizes = []
    for file in allFiles:
        parts = file.split("_")
        for part in parts:
            if part.startswith("LC"):
                memorySizes.append(int(part[2:]))

    memorySizes = sorted(set(memorySizes))
    print(memorySizes)

    starterFiles = list(filter(lambda x: x.endswith("_0" + FILE_EXTENSION), allFiles))
    memoryGroupedFileList = [[basePath + starterFile for starterFile in starterFiles if f"_LC{m}_" in starterFile] for m in memorySizes]
    memoryGroupedFileList = list(map(sorted, memoryGroupedFileList))

    # sort groups by memory size
    memoryGroupedFileList = sorted(memoryGroupedFileList, key=lambda memoryGroup: int(list(filter(lambda x: x.startswith("LC"), memoryGroup[0].split("_")))[0][2:]))
    return [memorySizes, memoryGroupedFileList]


def extractAvgTimeline(pathToFolder, folder):
    global FILE_EXTENSION
    pathToFolder = pathFixer(pathToFolder)
    folder = pathFixer(folder)
    try: os.makedirs(pathToFolder + folder + "timeline")
    except FileExistsError as e: print (e)

    _, memoryGroupedFileList = init(pathToFolder, folder)

    allTimelines = []
    for memoryGroup in memoryGroupedFileList:
        memoryGroupTimelines = []
        for simulation in memoryGroup:
            # get average timeline
            simulationTimeline = []
            reader = ParquetChunkReader(simulation) if FILE_EXTENSION == ".parquet" else DatChunkReader(simulation)
            layout = np.zeros(reader.shape)
            layout[reader.metaData["gapPosition1"], :reader.metaData["NN"]] = 1
            reader.mask = np.array(layout, dtype=bool)

            for dataChunk in reader:
                sign0 = dataChunk
                avg = np.mean(sign0)
                simulationTimeline.append(avg)

            # save average timeline
            ## insert the metaData
            newPath = pathToFolder + folder + "timeline/" + simulation.split("/")[-1].replace(FILE_EXTENSION, "") + ".tml"
            print(str(reader.metaData))
            with open(newPath, "w") as f: f.write(str(reader.metaData) + "\n")
            ## append with numpy array
            with open(newPath, "a") as f: np.savetxt(f, simulationTimeline)

            memoryGroupTimelines.append(simulationTimeline)
        allTimelines.append(memoryGroupTimelines)


base = pathFixer(r".\outData\dynamicMem32_allneg")
allDensityPaths = [s for s in getAllDensityPaths(base)]
for path in allDensityPaths:
    singles = (path[:-1] if path.endswith("/") else path).split("/")
    pathToFolder = "/".join(singles[:-1])
    folder = singles[-1]
    extractAvgTimeline(pathToFolder, folder)
