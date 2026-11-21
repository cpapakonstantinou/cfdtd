import glob
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

# current directory
here = os.path.abspath(os.path.dirname(__file__))


def latest(out, pattern):
	"""the file of the last step matching `pattern` in <out>/reports"""
	files = glob.glob(os.path.join(out, "reports", pattern))
	if not files:
		sys.exit(f"no {pattern} in {out}/reports -- usage: check.py <output directory of the run>")
	return max(files, key=lambda f: int(re.findall(r"(\d+)\.txt$", f)[0]) if re.findall(r"(\d+)\.txt$", f) else 0)


def main():
	simFile = latest(sys.argv[1] if len(sys.argv) > 1 else here + "/OUTPUT", "art_tfsf_[0-9]*.txt")

	# initialize figure:
	qFig, qAx = plt.subplots()

	# Load in and plot simulation Result
	simWL, _, R, _ = np.loadtxt(simFile, unpack=True)

	qAx.plot(simWL, R, 'o', label="R")

	qAx.set_xlabel("Wavelength (nm)")
	qAx.set_title("TFSF Reflection Coefficient")
	qAx.legend()
	qAx.grid()
	qFig.tight_layout()

	maxDevR = np.nanmax(np.abs(R))
	meanDevR = np.nanmean(np.abs(R))

	print(f"Max Absolute Deviation in R:\t{maxDevR:.8f}")
	print(f"Mean Absolute Deviation in R:\t{meanDevR:.8f}")

	# show results
	plt.show()

	return


if __name__ == "__main__":
	main()
