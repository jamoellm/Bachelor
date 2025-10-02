import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ZZ_basicUtility import labelTranslator

plt.rcParams.update({'font.size': 15})

fixKey = "density"
fixKeyValue = 0.5

notFixKey = "LC"

scenario = "allneg"

curmudgeons = np.array([0, 2, 4, 6, 8])

try:
    fixKeyValue = float(sys.argv[1])
    scenario = sys.argv[2]
except:
    fixKeyValue = 1.0
    scenario = "allneg"



# load all data
basePath = f"data/polarizability/staticMem32" #?
curmSuffix = lambda i: f"_C{int(i)}" if i != 0 else "" #?

polarizability = "AbsMean" # "MeanAbs"

dfs = [pd.read_csv(f"{basePath}{curmSuffix(i)}/all{polarizability}Polarizabilities_reduced.txt") for i in curmudgeons]


# # keep only relevant data
for i, df in enumerate(dfs):
    dfs[i] = df[df[fixKey] == fixKeyValue]


# fetch axis values
notFixKeyValues = np.array(dfs[0][notFixKey])

X, Y = np.meshgrid(curmudgeons, notFixKeyValues)
Z = np.array([df["polarizability"] for df in dfs])

# ------------------------------------------------------------------------
fig2, ax2 = plt.subplots(1, 1, figsize=(6,5))
tit = "$\\widetilde\\mathscr{P}$" if "MeanAbs" in polarizability else "$\\mathscr{P}$" if "AbsMean" in polarizability else "WeeWoo"
# extent legt Achsenbeschriftung anhand curmudgeons,y fest
plt.imshow(Z.T, origin='lower',
           extent=[X.min(), X.max(), Y.min(), Y.max()],
           aspect='auto', cmap='viridis',
           vmin=0, vmax=1)
plt.colorbar(label=f"reduced {tit}")
ax2.set_xlabel("#curmudgeon")
ax2.set_ylabel(labelTranslator(notFixKey))
ax2.set_title(f"{fixKey} {fixKeyValue}")


dxTickPos = (max(curmudgeons) - min(curmudgeons)) / len(curmudgeons)
xTickPos = np.arange(min(curmudgeons), max(curmudgeons), dxTickPos) + dxTickPos / 2
ax2.set_xticks(xTickPos)
ax2.set_xticklabels(curmudgeons.astype(int))

ax2.set_yticks((notFixKeyValues[:-1] + 1).tolist() + [20] ) #! designed for M / LC
ax2.set_yticklabels(notFixKeyValues.astype(int))

plt.tight_layout()

base = "staticMem32/" #?
print(base)
fig2.savefig(f"./graphs_v2/{base}3D_{polarizability}Pol(Curm, {notFixKey})_{fixKey}{fixKeyValue}{"_reduced" if True else ""}.svg")

# plt.show()