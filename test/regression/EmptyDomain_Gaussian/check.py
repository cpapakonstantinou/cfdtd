import glob
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

# current directory
here = os.path.abspath(os.path.dirname(__file__))

"""
DESCRIPTION:
A python script to compare the results of a Thin Film test to that
expected analytically

IMPORTANT:
Assumes that the background media is vacuum. Can be changed by editing
'n1' variable in the ThinFilmCalculator function, however should be
expanded to read from the matini file in the future
Assumes no absorption
"""
import os
import numpy as np
import matplotlib.pyplot as plt

# current directory
here = os.path.abspath(os.path.dirname(__file__))








def compareResults(r1, t1, r2, t2):
	# calculate some statistics
	relMaxDevT = np.nanmax(np.abs(t1 - t2) / (0.5 * (t1 + t2))) * 100
	relMeanDevT = np.nanmean(np.abs(t1 - t2) / (0.5 * (t1 + t2))) * 100
	maxDevT = np.nanmax(np.abs(t1 - t2))
	meanDevT = np.nanmean(np.abs(t1 - t2))

	
	maxDevR = np.nanmax(np.abs(r1 - r2))
	meanDevR = np.nanmean(np.abs(r1 - r2))

	print(f"Max Relative Deviation in T:\t{relMaxDevT:.2f}%")
	print(f"Mean Relative Deviation in T:\t{relMeanDevT:.2f}%")
	print(f"Max Absolute Deviation in T:\t{maxDevT:.6f}")
	print(f"Mean Absolute Deviation in T:\t{meanDevT:.6f}\n")

	
	print(f"Max Absolute Deviation in R:\t{maxDevR:.6f}")
	print(f"Mean Absolute Deviation in R:\t{meanDevR:.6f}\n")


def latest(out, pattern):
	"""the file of the last step matching `pattern` in <out>/reports"""
	files = glob.glob(os.path.join(out, "reports", pattern))
	if not files:
		sys.exit(f"no {pattern} in {out}/reports -- usage: check.py <output directory of the run>")
	return max(files, key=lambda f: int(re.findall(r"(\d+)\.txt$", f)[0]) if re.findall(r"(\d+)\.txt$", f) else 0)


def main():

	# load in previous results:
	oldFile = f"{here}/OldResults.txt"
	oldWL, oldT, oldR, _ = np.loadtxt(oldFile, unpack=True, skiprows=1)

	# Load in the simulation results
	simFile = latest(sys.argv[1] if len(sys.argv) > 1 else here + "/OUTPUT", "ART-TFSF_double_90pol_*.txt")
	simWL, simT, simR, _ = np.loadtxt(simFile, unpack=True, skiprows=1)
	simWL *= 1e-9

	plotDict = {
		'wl': [simWL, simWL, oldWL * 1e-9, oldWL * 1e-9],
		'coeff': [simT, simR, oldT, oldR],
		'fmt': ['bo', 'go', 'b--.', 'g--.'],
		'label': ['Sim T', 'Sim R', 'Old T', 'Old R']
	}

	

	# perform comparison to previous results:
	print("\nComparison to Previous:")
	compareResults(simR, simT, oldR, oldT)

	# initialize figure:
	qFig, qAx = plt.subplots()

	nPlots = len(plotDict['wl'])
	for i in range(nPlots):
		qAx.plot(plotDict['wl'][i] / 1e-9, plotDict['coeff'][i],
				 plotDict['fmt'][i], label=plotDict['label'][i])

	qAx.set_xlabel("Wavelength (nm)")
	qAx.set_title("Gaussian Test")
	qAx.legend()
	qAx.grid()
	qFig.tight_layout()

	# show results
	plt.show()

	return






if __name__ == "__main__":
	main()
