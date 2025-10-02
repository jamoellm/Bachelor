import ast
import sys
import pandas as pd
from ZZ_basicUtility import getMetaData, pathFixer



filePath = sys.argv[1]
metaData = ast.literal_eval(sys.argv[2])


skip = 1 if "_0.dat" in filePath else 0
df = pd.read_csv(filePath, sep=" ", skiprows=skip, names=list( range(metaData["maxLineLength"] + 1) ), dtype=float)
df.dropna(axis=1, how="all", inplace=True)
df.reset_index(drop=True, inplace=True)
print(filePath)
if "_0.dat" in filePath:
    with open(filePath) as F: metaDataLine = F.readline().strip()
    metaData = getMetaData(filePath)
    print(f",chunkSize {metaData["chunkSize"]}" )
    # ! always needed
    extraMetaData = f",chunkSize {metaData["chunkSize"]}" + f",gapPosition1 {metaData["gapPosition1"]}"

    metaDataLine = metaDataLine + extraMetaData
    print(metaDataLine)
    df.attrs["metaDataLine"] = metaDataLine

newfile = filePath.removesuffix(".dat") + ".parquet"
df.to_parquet(newfile)
print(newfile)

