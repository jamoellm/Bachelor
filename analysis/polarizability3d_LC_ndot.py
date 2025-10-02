import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.optimize import least_squares, curve_fit
from ZZ_basicUtility import labelTranslator, pathFixer

print(plt.rcParams["font.size"])
plt.rcParams.update({'font.size': 15})


xKey = "density"
yKey = "LC"

scenario = "allneg"


try:
    scenario = sys.argv[1]
    polarizability = sys.argv[2]
except:
    scenario = "allneg"
    polarizability = "MeanAbs" # "AbsMean"

print(scenario, polarizability)


# load all data
basePath = f"data/polarizability/dynamicMem32_{scenario}" #?
df = pd.read_csv(f"{basePath}/all{polarizability}Polarizabilities_reduced.txt")
xValues = np.array(sorted(list(set(list(df[xKey])))))
yValues = np.array(sorted(list(set(list(df[yKey])))))
if (xValues == xValues.astype(int)).all(): xValues = xValues.astype(int)
if (yValues == yValues.astype(int)).all(): yValues = yValues.astype(int)


X, Y = np.meshgrid(xValues, yValues)
Z = np.array([df[df[xKey] == xValue]["polarizability"] for xValue in xValues]).T

# ------------------------------------------------------------------------
fig2, ax2 = plt.subplots(1, 1, figsize=(6,5))
tit = "$\\widetilde\\mathscr{P}$" if "MeanAbs" in polarizability else "$\\mathscr{P}$" if "AbsMean" in polarizability else "WeeWoo"
plt.imshow(Z, origin='lower',
           extent=[X.min(), X.max(), Y.min(), Y.max()],
           aspect='auto', cmap='viridis',
           vmin=0, vmax=1)
plt.colorbar(label=tit)
ax2.set_xlabel(labelTranslator(xKey))
ax2.set_ylabel(labelTranslator(yKey))
ax2.set_title(f"reduced {tit}")


dxTickPos = (max(xValues) - min(xValues)) / len(xValues)
xTickPos = np.arange(min(xValues), max(xValues), dxTickPos) + dxTickPos / 2
ax2.set_xticks(xTickPos)
ax2.set_xticklabels(xValues)

dyTickPos = (max(yValues) - min(yValues)) / len(yValues)
yTickPos = np.arange(min(yValues), max(yValues), dyTickPos) + dyTickPos / 2
ax2.set_yticks(yTickPos)
ax2.set_yticklabels(yValues)

base = "trueDynamicMem/" + pathFixer(basePath.split("/polarizability/")[-1]) + "dIV0.300000/" #?
print(base)
os.makedirs(f"./graphs_v2/{base}", exist_ok=True)
plt.tight_layout()
fig2.savefig(f"./graphs_v2/{base}3D_{polarizability}Pol({xKey}, {yKey}){"_reduced" if True else ""}.svg")

# plt.show()