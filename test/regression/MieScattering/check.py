"""
DESCRIPTION:
A python script to compare the results of a Mie Scattering test to that
expected analytically

IMPORTANT:
I (Jesse) don't know how to do a Mie scattering myself
analytically yet, so I just took a result off web calculator.
That means it is currently limited to a 20nm radius sphere
"""
import glob
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

# current directory
here = os.path.abspath(os.path.dirname(__file__))


def Compare(comparePairs):
	for label, sim, true in comparePairs:
		relMaxDev = np.nanmax(np.abs(sim - true) / (0.5 * (sim + true))) * 100
		relMeanDev = np.nanmean(np.abs(sim - true) / (0.5 * (sim + true))) * 100
		maxDev = np.nanmax(np.abs(sim - true))
		meanDev = np.nanmean(np.abs(sim - true))
		print(f"Max Relative Deviation in {label}:\t{relMaxDev:.2f}%")
		print(f"Mean Relative Deviation in {label}:\t{relMeanDev:.2f}%")
		print(f"Max Absolute Deviation in {label}:\t{maxDev:.6f}")
		print(f"Mean Absolute Deviation in {label}:\t{meanDev:.6f}\n")
	return


def latest(out, pattern):
	"""the file of the last step matching `pattern` in <out>/reports"""
	files = glob.glob(os.path.join(out, "reports", pattern))
	if not files:
		sys.exit(f"no {pattern} in {out}/reports -- usage: check.py <output directory of the run>")
	return max(files, key=lambda f: int(re.findall(r"(\d+)\.txt$", f)[0]) if re.findall(r"(\d+)\.txt$", f) else 0)


def main():
	sphereRadius = 20e-9
	mieFile = f"{here}/MieTheoryResult.csv"
	oldFile = f"{here}/OldResults.txt"
	simFile = latest(sys.argv[1] if len(sys.argv) > 1 else here + "/OUTPUT", "Qcoeff_double_90pol.txt")

	# initialize figure:
	qFig, qAx = plt.subplots()

	# Load in and plot Analytical Result
	mieWL, cs_Ex, cs_Sc, cs_Ab = np.loadtxt(mieFile, unpack=True, skiprows=1, delimiter=',')
	cs_sphere = np.pi * sphereRadius**2
	for tmpData, tmpLabel, tmpFMT in zip([cs_Ex, cs_Sc, cs_Ab], ["Mie A", "Mie T", "Mie R"], ['k-', 'b-', 'g-']):
		qAx.plot(mieWL, tmpData / cs_sphere, tmpFMT, label=tmpLabel)
	
	# Load in and plot Old Results
	oldWL, oldT, oldR, oldA = np.loadtxt(oldFile, unpack=True)
	for tmpData, tmpLabel, tmpFMT in zip([oldT, oldR, oldA], ["Old T", "Old R", "Old A"], ['b--', 'g--', 'k--']):
		qAx.plot(oldWL, tmpData, tmpFMT, label=tmpLabel)
	
	# Load in and plot simulation Result
	simWL, T, R, A = np.loadtxt(simFile, unpack=True, skiprows=1)
	for tmpData, tmpLabel, tmpFMT in zip([T, R, A], ["Sim T", "Sim R", "Sim A"], ['bo', 'go', 'ko']):
		qAx.plot(simWL, tmpData, tmpFMT, label=tmpLabel)

	qAx.set_xlabel("Wavelength (nm)")
	qAx.set_ylabel("Q")
	qAx.set_title("Mie Theory Comparison of T, R, A")
	qAx.legend()
	qAx.grid()
	qFig.tight_layout()

	# calculate some statistics
	print("\n\nCOMPARISON TO THEORY:")
	interpT = np.interp(simWL, mieWL, cs_Sc / cs_sphere)
	interpR = np.interp(simWL, mieWL, cs_Ab / cs_sphere)
	interpA = np.interp(simWL, mieWL, cs_Ex / cs_sphere)
	comparePairs = [("T", T, interpT), ("R", R, interpR), ("A", A, interpA)]
	Compare(comparePairs)


	print("\n\nCOMPARISON TO PREVIOUS RESULTS:")
	comparePairs = [("T", T, oldT), ("R", R, oldR), ("A", A, oldA)]
	Compare(comparePairs)


	# show results
	plt.savefig('mie.png')

	return


if __name__ == "__main__":
	main()
