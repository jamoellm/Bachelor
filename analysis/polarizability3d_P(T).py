import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.optimize import least_squares, curve_fit
from ZZ_basicUtility import out

def chi_2(x, y, dy, func, para):
    return sum([((y[i] - func(x[i],*para)) / dy[i])**2 for i in range(len(x))]) / (len(x)-len(para))

def chi_2(x, z, dz, func, para):
    diff = z - func(x, *para)
    chi = (diff / dz)**2
    ndof = x.size - len(para)
    return np.sum(chi) / ndof

_df = pd.read_csv("data/polarizability/staticMem32/allAbsMeanPolarizabilities.txt")
_df_reduced = pd.read_csv("data/polarizability/staticMem32/allAbsMeanPolarizabilities_reduced.txt")
reduced = True
if not reduced:
    df = _df
if reduced:
    df = _df_reduced


#enter collision rate
df["ndot"] = np.zeros_like(df["density"])
densities = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])# * 18.8
totalCollisions = [7.0, 21.125, 369.875, 2827.25, 3835.375, 5309.4375, 6474.125, 8369.625, 10072.875, 11782.625, 13434.171875] # averaged over 64 simulations
totalCollisionsPerMinutePerSpinner = np.array(totalCollisions) / 1800*60/32
for d, ndot in zip(densities, totalCollisionsPerMinutePerSpinner):
    df.loc[df["density"] == d, "ndot"] = ndot


print("--------")
print(df)
print("--------")



fig, axs = plt.subplots(1, 3, subplot_kw={"projection": "3d"}, dpi=250, figsize=(6,2.5))

k1, k2 = "ndot", "LC"
matrix = df.pivot(index=k1, columns=k2, values="polarizability").values.T
dmatrix = df.pivot(index=k1, columns=k2, values="dpolarizability").values.T
ndots = np.array(sorted(list(set(df["ndot"].to_list()))))
LC = np.array(sorted(list(set(df[df.keys()[1]].to_list()))))


###### ZIEL: fit a * (mesured polarizability)^b = LC/ndot
def P_T(T, a, b):
    return a * np.pow(T, b)

measuredPolarizability = df.pivot(index=k1, columns=k2, values="polarizability")
dmeasuredPolarizability = df.pivot(index=k1, columns=k2, values="dpolarizability")

df_mess = df[["ndot", "LC", "polarizability", "dpolarizability"]]
df_mess = df_mess.rename(columns={"polarizability": "P", "dpolarizability": "dP"})

df_theo = pd.DataFrame({ "ndot": df["ndot"], "LC": df["LC"], "T": np.zeros_like(df["LC"]) })
df_theo["T"] = df_theo["LC"] / df_theo["ndot"]

para, pcov = curve_fit(P_T, df_theo["T"], df_mess["P"])
df_fit = pd.DataFrame({ "ndot": df["ndot"], "LC": df["LC"], "T'": np.zeros_like(df["LC"]) })
df_fit["P_fit"] = P_T(df_theo["T"], *para)
chi2 = chi_2(df_theo["T"], df_mess["P"], df_mess["dP"], P_T, para)

out(para, pcov, chi2)
print("______")


#### plotting ####
dx = min(np.diff(ndots))
dy = min(np.diff(LC))

## measured values ##
P_measure = df_mess["P"]
axs[0].bar3d(df_mess["ndot"], df_mess["LC"], np.zeros(df_mess["P"].size), dx, dy, P_measure)
axs[0].set_zlabel("$\\mathscr{P}$")
axs[0].set_title("measured $\\mathscr{P}$")

## naiive calculations (T) plotted ##
P_naiiveFit = P_T(df_theo["T"], 1,1)
axs[1].bar3d(df_theo["ndot"], df_theo["LC"], np.zeros(df_theo["T"].size), dx, dy, P_naiiveFit)
axs[1].set_zlabel("$\\mathscr{T}$")
axs[1].set_title("$\\mathscr{T}$")

## best fit of a * (M/N')^b onto df_mess["P"]
P_bestFit = P_T(df_theo["T"], *para)
axs[2].bar3d(df_fit["ndot"], df_fit["LC"], np.zeros(df_fit["P_fit"].size), dx, dy, P_bestFit)
axs[2].set_zlabel("$\\mathscr{P}$")
axs[2].set_title("fitted $\\mathscr{P}(\\mathscr{T})$")

for ax in axs.flatten():
    ax.set_xlabel("$\\dot{\\overline{N}}$")
    ax.set_ylabel("M")

fig.subplots_adjust(wspace=0.5)
fig.savefig(f"graphs_v2/staticMem32/polarizabilityComparison{"_reduced" if reduced else ""}")
# plt.show()

