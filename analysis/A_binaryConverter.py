import pandas as pd
import os
from ZZ_basicUtility import getAllDensityPaths, getMetaData, pathFixer


def doBinaryConverting(pathToFolder: str, folder: str):
    fullPath = pathFixer(pathToFolder) + pathFixer(folder)
    allFiles = list(filter(lambda x: x.endswith(".dat"), os.listdir(fullPath)))
    print(allFiles)
    print(fullPath + list(filter(lambda x: ".dat" in x, allFiles))[0])

    metaData = getMetaData(fullPath + list(filter(lambda x: "_0.dat" in x, allFiles))[0])
    allFiles = sorted(list(filter(lambda x: x.endswith(".dat"), allFiles)))

    print(metaData)
    for f in allFiles:
        file = fullPath + f
        skip = 1 if "_0.dat" in f else 0
        df = pd.read_csv(file, sep=" ", skiprows=skip, names=list( range(metaData["maxLineLength"] + 1) ), dtype=float)
        df.dropna(axis=1, how="all", inplace=True)
        df.reset_index(drop=True, inplace=True)
        print(file)
        if "_0.dat" in f:
            with open(file) as F: metaDataLine = F.readline().strip()
            metaData = getMetaData(file)
            print(f",chunkSize {metaData["chunkSize"]}" )
            # ! always needed
            extraMetaData = f",chunkSize {metaData["chunkSize"]}" + f",gapPosition1 {metaData["gapPosition1"]}"

            metaDataLine = metaDataLine + extraMetaData
            print(metaDataLine)
            df.attrs["metaDataLine"] = metaDataLine

        newfile = fullPath + f.replace(".dat", ".parquet")
        df.to_parquet(newfile)
        # print(newfile)



base = r".\outData\staticMem32"
allDensityPaths = [s for s in getAllDensityPaths(base)]
for path in allDensityPaths:
    singles = (path[:-1] if path.endswith("/") else path).split("/")
    pathToFolder = "/".join(singles[:-1])
    folder = singles[-1]
    # print(path)
    # print(pathToFolder, folder)
    doBinaryConverting(pathToFolder, folder)
