import numpy as np
import pandas as pd
class ParquetChunkReader:
    def __init__(self,
                 filePath: str,
                 *,
                 skipFrames: int = 0,
                 preFormatFunction = None,
                 ):

        self.path = filePath
        self.skipFrames = skipFrames
        self.preFormatFunction = preFormatFunction
        self.mask = None

        self.metaData = {}
        self.partIdx = 0
        self.i = 0
        self.globalTime = 0

        df = pd.read_parquet(self.path)
        self.metaData = self.get_metaData(df)
        self.data = df.to_numpy()

        self.shape = (self.metaData["chunkSize"], self.data.shape[1])
        self.size = self.shape[0] * self.shape[1]

    def __iter__(self):
        assert self.mask is not None
        return self

    def __next__(self):
        try:
            self.data[self.i * self.metaData["chunkSize"]]
        except:
            self._openNextFile()

        chunk = self.data[self.i * self.metaData["chunkSize"]: (self.i+1) * self.metaData["chunkSize"]]
        self.i += 1 + self.skipFrames
        self.globalTime += 1 + self.skipFrames

        if self.preFormatFunction is not None:
            try:
                fetch = chunk[self.mask]
                # fetch = fetch[~np.isnan(fetch)]
            except:
                raise RuntimeError("Mask could not properly applied.")
            try:
                return self.preFormatFunction(fetch)
            except:
                raise RuntimeError("'preFormatFunction' isn't doing its thing properly.")
        # print("------", chunk)
        # fetch = chunk[self.mask]
        # fetch = fetch[~np.isnan(fetch)]
        return chunk[self.mask]

    def __getitem__(self, chunkIndex: int):
        chunksPerFile = self.data.shape[0] / self.metaData["chunkSize"]
        if chunkIndex >= 0:
            if chunkIndex > chunksPerFile:
                k = chunkIndex // chunksPerFile
                self.partIdx += int(k - 1)
                self.globalTime
                self._openNextFile()
                return self.__getitem__(int(chunkIndex - k * chunksPerFile))
            else:
                self.i = chunkIndex
                self.globalTime = self.partIdx * chunksPerFile + chunkIndex
                return self.__next__()
        elif chunkIndex < 0:
            maxIndex = self.metaData["n_files"] * self.data.shape[0] / self.metaData["chunkSize"]
            i2 = int(maxIndex - abs(chunkIndex))
            return self.__getitem__(i2)[0] # ! could be bad

    def _openNextFile(self):
        try:
            self.partIdx += 1
            self.path = ".".join(self.path.split(".")[:-1])[:-1] + str(self.partIdx) + ".parquet"
            print(self.partIdx, self.path)
            self.data = pd.read_parquet(self.path).to_numpy()
            self.i = 0
        except:
            print("Done")
            raise StopIteration

    def _openPreviousFile(self):
        self.partIdx -= 1
        self.path = ".".join(self.path.split(".")[:-1])[:-1] + str(self.partIdx) + ".parquet"
        print(self.partIdx, self.path)
        self.data = pd.read_parquet(self.path).to_numpy()
        self.i = 0

    def get_metaData(self, df):
        metaData = {}
        metaDataLine = df.attrs["metaDataLine"]

        data = metaDataLine.split(":")[1] # "T 1800,dt 1e-4"
        Ldata = data.split(",") # ["T 1800", "dt 1e-4"]
        for qty in Ldata:
            s = qty.split(" ")
            if float(s[1].strip()) == float(s[1])//1: metaData[s[0]] = int(float(s[1])//1)
            else: metaData[s[0]] = float(s[1])

        metaData["maxLineLength"] = max(metaData["NN"], metaData["LC"]) # + 1 # ! caution
        return metaData

    def make_zeroMask(self):
        self.mask = np.zeros(shape=(self.metaData["chunkSize"], self.data.shape[1]))

import ast
class TmlReader:
    def __init__(self, filePath: str):
        # assert filePath.endswith(".tml")

        with open(filePath) as f: line0 = f.readline()
        try:
            self.metaData = ast.literal_eval(line0)
        except:
            self.metaData = {}
            data = line0.split(":")[1] # "T 1800,dt 1e-4"
            Ldata = data.split(",") # ["T 1800", "dt 1e-4"]
            for qty in Ldata:
                s = qty.split(" ")
                if float(s[1].strip()) == float(s[1])//1: self.metaData[s[0]] = int(float(s[1])//1)
                else: self.metaData[s[0]] = float(s[1])

        if filePath.endswith(".tml"):
            with open(filePath) as f:
                f.readline()
                self.data = np.loadtxt(f)

        elif filePath.endswith(".tmlnp"):
            with open(filePath, "rb") as f:
                f.seek(len(line0) + 1)
                self.data = np.load(f)

        self.trimmedZeros = 0

from ZZ_basicUtility import getMetaData, pathFixer
class DatChunkReader:
    def __init__(self,
                 filePath: str,
                 *,
                 preFormatFunction = None
                 ):

        self.path = filePath
        self.mask = None

        self.metaData = getMetaData(self.path)
        self.partIdx = 0

        self.file = open(self.path)
        self.file.readline() # remove metaDataLine

        self.shape = (self.metaData["chunkSize"], self.metaData["maxLineLength"])
        self.size = self.shape[0] * self.shape[1]

        self.preFormatFunction = preFormatFunction

    def __iter__(self):
        assert self.mask is not None
        return self

    def __next__(self):
        try:
            i = 0
            cc = 0
            data = []
            while i < self.shape[0]:
                l = self.file.readline().strip()
                if l == "":
                    cc += 1
                    if cc >= 3:
                        raise EOFError()
                    continue
                else:
                    cc = 0
                    l = list(map(float, l.split(" ")))
                    data.append(l)
                    i += 1
            # print(*data, sep="\n")
        except:
            self._openNextFile()
            return self.__next__()

        arr = np.full((self.metaData["chunkSize"], self.metaData["maxLineLength"]), np.nan)
        for i, row in enumerate(data):
            arr[i, :len(row)] = row


        if self.preFormatFunction is not None:
            return self.preFormatFunction(arr[self.mask])
        else:
            return arr[self.mask]

    def _openNextFile(self):
        try:
            self.partIdx += 1
            self.path = ".".join(self.path.split(".")[:-1])[:-1] + str(self.partIdx) + ".dat"
            print(self.partIdx, self.path)
            self.file = open(self.path)
        except:
            print("Done")
            raise StopIteration

def ChunkReader(filePath: str,
                *,
                skipFrames: int = 0,
                preFormatFunction = None,
                ):
    if filePath.endswith(".parquet"):
        return ParquetChunkReader(filePath, skipFrames=skipFrames, preFormatFunction=preFormatFunction)
    elif filePath.endswith(".dat"):
        if skipFrames != 0 or preFormatFunction is not None:
            print("Warning: .dat files does not support `skipFrames`.")
        return DatChunkReader(filePath, preFormatFunction=preFormatFunction)


if __name__ == "__main__":
    path = pathFixer(r".\outData\staticMem\density1.0")
    file = path + "spinner_T1800_LC21_density0.100000_s3_0.dat"
    reader = ChunkReader(file)
    print(type(reader))
    mask = np.zeros(reader.shape)
    mask[reader.metaData["gapPosition1"], :reader.metaData["NN"]] = 1
    mask = mask.astype(bool)
    print(mask)

    reader.mask = mask
    for idx, chunk in enumerate(reader):
        if idx >= 5: break
        print(chunk)


