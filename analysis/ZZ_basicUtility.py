import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

errstyle = {"fmt":".","ms":7,"ecolor":"k", "elinewidth":1, "capsize":3, "capthick":1}

def makeZeroChunk(filepath: str)->list[list]:
    ""
    lines = []
    c = 0
    with open(filepath, "r") as file:
        while True:
            l = file.readline().strip()
            if l.startswith("meta"):
                continue
            elif l == "":
                c += 1
            else:
                c = 0
            lines.append(l.split(" "))
            if c == 2:
                lines.pop()
                lines.pop()
                break

    res = [[] for _ in lines]
    for i in range(len(lines)):
        res[i] = [0 for _ in lines[i] if _ != ""]


    return [r for r in res if r != []]

def make1At(M: list[list], positions: list[tuple | int])-> list[list]:
    """takes a list of `tuples`"""
    for pos in positions:
        M[pos[0]][pos[1]] = 1
    return M

def layout2mask(layout, maxLineLength, returnType: str) -> list[tuple]:
    mask = np.array([arr + [0] * (maxLineLength - len(arr)) for arr in layout])
    if returnType == "positions":
        return np.argwhere(mask == 1)
    elif returnType == "mask":
        return np.array(mask, dtype=bool)
    else:
        raise ValueError("returnType must be 'postitions' or 'mask'")

def getMetaData(path: str):
    with open(path, "r") as file:
        metaData = {}
        linebreaks = 0
        lineCounter = 2
        maxLineLength = 0
        while True:
            l = file.readline().strip() # todo: changed - added strip()
            if l.startswith("meta") or l.startswith("#meta"):
                meta = l.split(":")[1]
                Ldata = meta.split(",") # ["T 1800", "dt 1e-4"]
                for qty in Ldata:
                    s = qty.split(" ")
                    value = float(s[1])
                    if value == int(value): metaData[s[0]] = int(value)
                    else: metaData[s[0]] = value
            else:
                maxLineLength = max(maxLineLength, len(l.split(" ")))

                if l.strip() == "":
                    linebreaks += 1
                    if "gapPosition1" not in metaData.keys():
                        metaData["gapPosition1"] = lineCounter - 2
                        print(lineCounter - 2, "nils")
                else:
                    lineCounter += 1
                    linebreaks = 0
                if linebreaks == 2:
                    lineCounter = lineCounter - 2
                    break

        # chunkSize = lineCounter - 2
        metaData["chunkSize"] = lineCounter
        metaData["maxLineLength"] = maxLineLength
    if "T" in path:
        fileT = float([s for s in path.split("_") if "T" in s][0][1:])
        if metaData["T"] != fileT:
            print(f"Warning: ({metaData["T"]} = metaData['T']) != (fileT = {fileT})")
            metaData["T"] = fileT
    else:
        print("Warning: File does not contain time in its title:", path)
    return metaData

def getSignIndex(path: str):
    f = open(path, "r")
    ll = 0
    cc = 0
    while True:
        line = f.readline()
        if line.startswith("meta"): continue
        if line.strip() == "":
            f.close()
            return ll
        else:
            ll += 1

def getChunkSize(path: str) -> int:
    """counts non-empty lines"""
    with open(path, "r") as file:
        ll = 0
        cc = 0
        while True:
            line = file.readline()
            if line.startswith("meta"): continue
            if line.strip() == "":
                cc += 1
            if cc == 2:
                return ll - 1
            ll += 1

def pathFixer(path: str) -> str:
    if path == "": return ""
    path = path.replace("\\", "/")
    return path if path[-1] == "/" else path + "/"

def _getAllDensityPaths(basePath: str, result: list):
    basePath = pathFixer(basePath)
    allFolders = [pathFixer(s) for s in os.listdir(basePath) if os.path.isdir(basePath + s)]
    for folder in allFolders:
        if folder.startswith("density"):
            result.append(basePath + folder)
        else:
            _getAllDensityPaths(basePath + folder, result)

    return result

def getAllDensityPaths(basePath: str):
    return _getAllDensityPaths(basePath, [])

def titleTranslator(title: str) -> str:
    match title:
        case "LC": return "M"
        case _: return title

def labelTranslator(label: str) -> str:
    match label:
        case "LC": return "memory size $M$"
        case _: return label

def chi_2(x, y, dy, func, para):
    return sum([((y[i] - func(x[i],*para)) / dy[i])**2 for i in range(len(x))]) / (len(x)-len(para))

def out(para, pcov, chi2):
    print("para:", para)
    perr = np.sqrt(np.diag(pcov))
    print("perr:", perr)
    print("pcov:", perr**2)
    print("chi2:", chi2)
    for p, e in zip(para, perr):
        print(np.round(p,3), "+-", np.round(e, 3))
    if chi2 is not None:
        print(np.round(chi2, 3))
    else:
        print(chi2)

def fitter(ax, fit, x, y, dy, dx=None, p0=None, errstyle=errstyle):
    x, y, dy = list(map(np.array, [x, y, dy]))
    dx = None if dx is None else np.array(dx)

    para, pcov = curve_fit(fit, x, y, p0=p0, sigma=dy, absolute_sigma=True)
    chi2 = chi_2(x, y, dy, fit, para)
    out(para, pcov, chi2)

    x_fit = np.linspace(min(x), max(x), 400)

    line, = ax.plot(x_fit, fit(x_fit, *para), label="fit", ls="--")
    ax.plot([], [], marker=">", color=line.get_color(), ls="", label=f"$\\chi^2_{{ndof}}=${np.round(chi2, 2)}")

    return para, pcov


if __name__ == "__main__":
    getAllDensityPaths(r".\outData\staticMem")
