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
import configparser
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# current directory
here = os.path.abspath(os.path.dirname(__file__))


def ThinFilmCalculator(n2, wl, matLen):
	'''
	DESCRIPTION
		An analytical thin film calculator
		(built using formulas from Hecht 2016)
	INPUTS
		n2:         Thin film refractive index
		wl:         Wavelengths to evaluate at
		matLen:     Thin film thickness
	'''
	# formulas from Hect 2016
	n1 = 1
	r = (n1 - n2) / (n1 + n2)
	delta = 4 * np.pi * n2 * matLen / wl
	F = (2 * r / (1 - r**2))**2
	It_Ii = 1 / (1 + F * np.sin(delta / 2)**2)
	return It_Ii


def GetAnalyticalFromConfig(iniFile, wl):
	'''
	DESCRIPTION
		Analytical T and R of the slab described in run.ini
		([geometry] middle_medium / middle_thickness, [materials] n<i>)
	INPUTS
		iniFile:    path to run.ini
		wl:         wavelengths to evaluate T and R at
	'''
	ini = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))
	if not ini.read(iniFile):
		raise FileNotFoundError(iniFile)
	medium = ini.getint("geometry", "middle_medium")
	matLen = ini.getfloat("geometry", "middle_thickness")
	matIndex = ini.getfloat("materials", f"n{medium}")

	trueT = ThinFilmCalculator(matIndex, wl, matLen)

	return trueT, 1 - trueT


def compareResults(r1, t1, r2, t2):
	# calculate some statistics
	relMaxDevT = np.nanmax(np.abs(t1 - t2) / (0.5 * (t1 + t2))) * 100
	relMeanDevT = np.nanmean(np.abs(t1 - t2) / (0.5 * (t1 + t2))) * 100
	maxDevT = np.nanmax(np.abs(t1 - t2))
	meanDevT = np.nanmean(np.abs(t1 - t2))

	relMaxDevR = np.nanmax(np.abs(r1 - r2) / (0.5 * (r1 + r2))) * 100
	relMeanDevR = np.nanmean(np.abs(r1 - r2) / (0.5 * (r1 + r2))) * 100
	maxDevR = np.nanmax(np.abs(r1 - r2))
	meanDevR = np.nanmean(np.abs(r1 - r2))

	print(f"Max Relative Deviation in T:\t{relMaxDevT:.2f}%")
	print(f"Mean Relative Deviation in T:\t{relMeanDevT:.2f}%")
	print(f"Max Absolute Deviation in T:\t{maxDevT:.6f}")
	print(f"Mean Absolute Deviation in T:\t{meanDevT:.6f}\n")

	print(f"Max Relative Deviation in R:\t{relMaxDevR:.2f}%")
	print(f"Mean Relative Deviation in R:\t{relMeanDevR:.2f}%")
	print(f"Max Absolute Deviation in R:\t{maxDevR:.6f}")
	print(f"Mean Absolute Deviation in R:\t{meanDevR:.6f}\n")


def main():

	# load in previous results:
	oldFile = f"{here}/OldResults.txt"
	oldWL, oldT, oldR, _ = np.loadtxt(oldFile, unpack=True, skiprows=1)

	# Load in the simulation results
	simFile = f"{sys.argv[1] if len(sys.argv) > 1 else here + '/OUTPUT'}/reports/art_4000.txt"
	simWL, simT, simR, _ = np.loadtxt(simFile, unpack=True)
	simWL *= 1e-9

	plotDict = {
		'wl': [simWL, simWL, oldWL * 1e-9, oldWL * 1e-9],
		'coeff': [simT, simR, oldT, oldR],
		'fmt': ['bo', 'go', 'b--.', 'g--.'],
		'label': ['Sim T', 'Sim R', 'Old T', 'Old R']
	}

	try:
		# calculate expected results from analytical formula
		refWL = simWL
		trueT, trueR = GetAnalyticalFromConfig(f"{here}/run.ini", simWL)
		plotDict['wl'].extend([refWL, refWL])
		plotDict['coeff'].extend([trueT, trueR])
		plotDict['fmt'].extend(['b--', 'g--'])
		plotDict['label'].extend(['Ref T', 'Ref R'])

		# perform comparison to analytical results:
		print("\nComparison to Theory:")
		compareResults(simR, simT, trueR, trueT)
	except FileNotFoundError:
		print("Cannot find 'run.ini', skipping analytical")

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
	qAx.set_title("Thin Film Test")
	qAx.legend()
	qAx.grid()
	qFig.tight_layout()

	# show results
	plt.show()

	return


if __name__ == "__main__":
	main()
