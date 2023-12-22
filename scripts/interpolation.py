import math
from scipy.stats import linregress
from scipy.optimize import curve_fit
import re
import numpy

# Fit is the main function to interpolate! (uses MCfit)
# data = yValues that are not interpolated
# X2 = xValues belonging to yValues
# FurtherPoints = additionalValues that should also be interpolated
def fit(data, X2=None, FurtherPoints=None, pre=""):
	l = len(data)
	X = list(range(l))
	if X2 is not None:
		X = X2

    # searches for the best a and b for a given (linear) function and the data
	popt, pcov = curve_fit(MCfit, X, data)

	a, b = popt[0], popt[1]

	I = X
	if FurtherPoints is not None:
		I = I + FurtherPoints

	#data = [(b + (a * x)) for x in I]
	data = MCfit(I, a, b)
	print(pre, "a * x + b ", a, b)
	return data, a, b

# Calculates for a given value set of n (here called x) the corresponding interpolated values using
# linear vars a and b that resulted from the fit function.
def MCfit(x, a, b):
	return [a*xx/(math.log(xx, 2)*5) + b for xx in x]