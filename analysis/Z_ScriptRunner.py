import subprocess

programms = ""
programms +=       "A_2Dplots.py"
programms += "," + "A_fourierAnalysis.py"
programms += "," + "A_fourierAnalysisSeedSorted.py"
programms += "," + "A_polarizability.py"
programms += "," + "A_timelinePlot.py"

programms = programms.split(",")

for p in programms:
    assert p.startswith("A_")
    assert p.endswith(".py")
    print(f"Starting {p} ...")
    subprocess.run(["python", p, r".\outData\staticMem32"])
