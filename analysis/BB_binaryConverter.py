import sys
import pandas as pd
import os
from ZZ_basicUtility import getMetaData, pathFixer
# from time import time

base = sys.argv[1]
base = pathFixer(base)
print(base)

allFiles = list(filter(lambda x: x.endswith(".dat"), os.listdir(base)))
print(base + list(filter(lambda x: ".dat" in x, allFiles))[0])

metaData = getMetaData(base + list(filter(lambda x: "_0.dat" in x, allFiles))[0])
allFiles = sorted(list(filter(lambda x: x.endswith(".dat"), allFiles)))

print(metaData)


for f in allFiles:
    ""
    file = base + f
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
    # print(df)

    newfile = base + ".".join(f.split(".")[:-1]) + ".parquet"
    df.to_parquet(newfile)
    print(newfile)

