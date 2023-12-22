from typing import Any
from matplotlib import rc
import matplotlib.pyplot as plt
import math
import Complexity as cmpl
from interpolation import *
from Instance import *

INTERPOINTS = [2129, 2265, 2472, 2681, 2893, 3108]
# Dont have benched prange, must interpolate that point... adjust these if necc
INTERPOINTSSP = [2129, 2197, 2265, 2472, 2681, 2893, 3108]

# Helper class for writeTimeTableToTexFile function
class TableRow:
	# valid data_types are: tu (time_unit), l (literal in logarithm)
	def __init__(self, label : str, data: list, data_type: str):
		self.label = label
		self.data = data
		self.data_type = data_type

# Helper function to join two lists while pairing their respective items
# eg a=[1, 2] b = [a, b] -> return [ (1, a), (2, b) ]
def mergeListsToTupleLists(a, b):
	return list(zip(a,b))

# Helper function to filter lines starting with '#' when reading benchmark data
def expectedFiler(line):
	if len(line) == 0:
		return False
	if line[0] == "#":
		return False
	return True

# Returns list of [n loops time_in_microseconds]
# file must consists of lines following the format above (example in Measurements folder)
def readExectedFile(filename: str):
	f = open(filename, 'r')
	lines = list(filter(expectedFiler, f.readlines()[1:]))
	#print(lines)
	s = [x.strip().split(" ") for x in lines]
	return [[int(x[0]), float(x[1]), float(x[2])] for x in s]

# Returns list of [n loops time_in_microseconds l p]
# file must consists of lines following the format above (example in Measurements folder)
def readDumerFile(filename: str):
	f = open(filename, 'r')
	lines = list(filter(expectedFiler, f.readlines()))
	s = [x.strip().split(" ") for x in lines]
	return [[int(x[0]), float(x[1]), float(x[2]), int(x[3]), int(x[4])] for x in s]

# Returns 3 lists [n],[time_in_seconds log2], [(n, time_in_seconds log2)]
# the first 2 lists are required to plot the data using matplotlib, the last list is a
# convenience list. This function is for Prange data and the second parameter decides if it
# is data for StandardPrange or TemplatePrange (important for complexity estimation)
def prepData(fileData, templateP = True):
	x = []
	y = []
	t = []
	for p in fileData:
		x.append(p[0])
		t_per_loop = p[2] / p[1] #in microseconds
		t_per_loop = t_per_loop / (1000 * 1000) # in seconds
		exp_perm = cmpl.prangeComplexitySpecialPerm(MapFromNToTemplateInstance[p[0]], False) if templateP else cmpl.prangeComplexityRandomPerm(MapFromNToStandardInstance[p[0]], False)
		print("{}: {}".format(p[0], exp_perm))
		#print(exp_perm)
		y.append(math.log(t_per_loop*exp_perm,2))
		t.append( (p[0], math.log(t_per_loop*exp_perm,2)) )

	return x,y,t

# Basically the same function as above, but just for dumer.
def prepDumerData(fileData):
	x = []
	y = []
	t = []
	for p in fileData:
		x.append(p[0])
		t_per_loop = p[2] / p[1] # in microseconds
		t_per_loop = t_per_loop / (1000*1000) # in seconds
		exp_perm = cmpl.dumerComplextityRandomPerm(MapFromNToTemplateInstance[p[0]], p[3], p[4], False)
		print(p[0], exp_perm, t_per_loop)
		y.append(math.log(t_per_loop*exp_perm,2))
		t.append( (p[0], math.log(t_per_loop*exp_perm,2)) )

	return x,y,t

# Prints the data prepared by preData or prepDumerData. Converts the estimated time
# to more or less every time unit.
def printPreparedData(n ,t, message = ""):
	print(message)
	nr_threads = 200.0
	for nn, tt in zip(n,t):
		s = "2^{}".format(round(tt,2))
		h = round((2**tt) / (60*60),2)
		d = round(h / 24, 2)
		y = round(d / 365,2)
		d_sc = round(d / nr_threads, 2)
		y_sc = round(y/ nr_threads, 2)
		print("{}: {}s = {}h = {}d = {}y (Scaled to {} threads: {}d, {}y)".format(nn, s, h, d, y, nr_threads, d_sc, y_sc))

# Computes the factor data1/data2 for each entry in these lists. Both lists must have
# the same size and must contain tuples in the form (n, time_in_seconds log2).
def cmpData(data1: list, data2: list, should_print = False):
	assert len(data1) == len(data2)
	data = []
	for t1, t2 in zip(data1, data2):
		factor = round((2**t1[1])/(2**t2[1]),2)
		factor_log = round(math.log(factor, 2),1)
		if should_print:
			print("{}: {} (2^{})".format(t1[0], factor, factor_log))
		data.append( (t1[0], factor_log) )
	return data

# Writes a latex table to "table.tex". The first parameter specifies which instances you
# want to be displayed. Then you can give any number of TableRows as arguments (not as a list)
# containg a label as a string, the data as list in the form [ (n, some_value)], and the type of
# some_value, e.g. if it is a time unit (runtime of n in sec log2) or a literal (speedup in log2).
# Additionally a caption can be given with a keyword argument caption="A cool caption".
def writeTimeTableToTexFile(instances: list, *args, **kwargs):
	def convert_time(t_s):
		flat_t = 2**t_s
		if flat_t < 60:
			return str(round(flat_t,2))+"s"
		flat_t /= 60.0
		if flat_t < 60:
			return str(round(flat_t, 2))+"m"
		flat_t /= 60.0
		if flat_t < 24:
			return str(round(flat_t, 2))+"h"
		flat_t /= 24.0
		if flat_t < 365:
			return str(round(flat_t, 2))+"d"
		flat_t /= 365.0
		if flat_t < 512:
			return str(round(flat_t, 2))+"y"
		else:
			return "$2^{" +str(round(math.log2(flat_t),2)) + "}$y"
		
	col_data_length = len(instances)
	row_length = len(args)
	# add all the begin environments and such
	latex_code = "\\begin{table}[h!]\n"
	latex_code += "\t\\centering\n"
	latex_code += "\t\\setlength\\extrarowheight{7pt}\n"
	latex_code += "\t\\begin{tabular}{|l||*{" + str(col_data_length) + "}{c|}}\hline\n"
	# add the table header
	latex_code += "\t\t\\backslashbox{Category}{n}"
	for n in instances:
		latex_code += " & " + str(n)
	latex_code += "\\\\\\hline\\hline\n"

	# add the table data (rows)
	for row in args:
		latex_code += "\t\t" + row.label
		for n in instances: # search the row_data if specified instance is available
			data_string = "-"
			for t in row.data:
				if t[0] == n: # found data, check if tu type or not
					if row.data_type == "tu":
						data_string = convert_time(t[1])
					else:
						data_string = "$2^{" + str(round(t[1],2))+"}$"
					break
			latex_code += " & " + data_string
		latex_code += "\\\\\\hline\n"
	
	# add the ending
	latex_code += "\t\\end{tabular}\n"
	latex_code += "\t\\caption{ " + kwargs.get("caption", "Here could be your caption!") + " }\n"
	latex_code += "\\end{table}\n"

	with open("table.tex", "w") as f:
		f.write(latex_code)

# IMPORTANT: Benched instances must be the same for all input data, e.g. if you benched the n=1995 instance for prange but
# not for dumer you are going to have a bad time 
if __name__ == "__main__":
	tStandardPrange = readExectedFile("Measurements/StandardPrangeTime.txt")
	tTemplatePrange = readExectedFile("Measurements/TemplatePrangeTime.txt")
	tTemplateDumerV1 = readDumerFile("Measurements/TemplateDumerV1Time.txt")
	print("--------------Prep Template Data------------------")
	x,yTP, tTP = prepData(tTemplatePrange)
	print("--------------Prep Standard Data------------------")
	x2, ySP, tSP = prepData(tStandardPrange, False)
	print("--------------Prep Template Dumer Data------------------")
	x3, yTD, tTD = prepDumerData(tTemplateDumerV1)
	#xInterpolationP = x + INTERPOINTS
	interpolated, a1, b1 = fit(yTP, x)
	interpolated2, a2, b2 = fit(ySP, x2)
	interpolated3, a3, b3 = fit(yTD, x3)
	missingInterpolTemplate = MCfit(INTERPOINTS, a1, b1)
	missingInterpolStandard = MCfit(INTERPOINTSSP, a2, b2)
	missingInterpolTemplateDumer = MCfit(INTERPOINTS, a3, b3)
	printPreparedData(x,yTP, "TemplatePrangeTime Samples")
	printPreparedData(INTERPOINTS, missingInterpolTemplate, "TemplatePrangeTime Interpol")
	printPreparedData(x2,ySP, "StandardPrangeTime Samples")
	printPreparedData(INTERPOINTSSP, missingInterpolStandard, "StandardPrangeTime Interpol")
	printPreparedData(x3,yTD, "TemplateDumerTime Samples")
	printPreparedData(INTERPOINTS, missingInterpolTemplateDumer, "TemplateDumerTime Interpol")

	mSP = mergeListsToTupleLists(INTERPOINTSSP, missingInterpolStandard)
	mTP = mergeListsToTupleLists(INTERPOINTS, missingInterpolTemplate)
	mTD = mergeListsToTupleLists(INTERPOINTS, missingInterpolTemplateDumer)

	allPoints = x + INTERPOINTS # in x is also 2197 wich is interpolated for now...
	finalDataTP = tTP + mTP
	finalDataSP = tSP + mSP
	finalDataTD = tTD + mTD

	# can sort like so
	finalDataTP.sort(key=lambda x: x[0])
	finalDataSP.sort(key=lambda x: x[0])
	finalDataTD.sort(key=lambda x: x[0])
	allPoints.sort()

	cmpSPTP = cmpData(finalDataSP, finalDataTP)
	cmpTPTD = cmpData(finalDataTP, finalDataTD)

	writeTimeTableToTexFile([640, 1101, 1473, 1995, 2472, 2893, 3108, 3488],
						 TableRow("Standard Prange", finalDataSP, "tu"),
						 TableRow("Template Prange", finalDataTP, "tu"),
						 TableRow("Template Dumer", finalDataTD, "tu"),
						 TableRow("Speedup (SP->TP)", cmpSPTP, "l"),
						 TableRow("Speedup (TP-> TD)", cmpTPTD, "l"),
						 caption="Example caption!")

	# ------ EVERYTHING FROM HERE ON IS CURRENTLY NOT TESTED, CODE TO CREATE GRAPHS AND PRINT FACTORS---------

	# # print("-------Factor (standard/template) Measured+Expected+Interpol---------")
	# # cmpData(x, interpolated2, interpolated)

	# # print("-------Factor (standard/template) Measured+Expected+Actual---------")
	# # cmpData(x, ySP, yTP)

	# print("-------Factor (standard/template) Missing+Expected+Interpol-----------")
	# cmpData(INTERPOINTS, missingInterpolStandard, missingInterpolTemplate)

	# print("-------Factor (TemplatePrange/TemplateDumer) Missing+Expected+Actual----------")
	# cmpData(x, yTP, yTD)

	# print("-------Factor (TemplatePrange/TemplateDumer) Missing+Expected+Interpol----------")
	# cmpData(INTERPOINTS, missingInterpolTemplate, missingInterpolTemplateDumer)

	# print("--------Direct Expected Permutation Comparison (Standard/Template)----------------")
	# for n, inst in MapFromNToStandardInstance.items():
	# 	standard_perm = cmpl.prangeComplexityRandomPerm(inst, False)
	# 	special_perm = cmpl.prangeComplexitySpecialPerm(MapFromNToTemplateInstance[n], False)
	# 	factor = standard_perm/special_perm
	# 	factor_log = round(math.log(factor,2), 1)
	# 	print("{}: 2^{} / 2^{} = 2^{}".format(n,round(math.log(standard_perm,2), 1), round(math.log(special_perm,2), 1), factor_log))


	# font = {'size': 9}
	# rc('font', **font)
	# plt.plot(x, y, 'b^', mfc='none', label = "TemplatePrange Samples")
	# plt.plot(x, interpolated, 'b--', mfc='none', label = "TemplatePrange Interpolation")
	# plt.plot(x2, y2, 'r^', mfc='none', label = "StandardPrange Samples")
	# plt.plot(x, interpolated2, 'r--', mfc='none', label = "StandardPrange Interpolation")
	# plt.plot(x3, y3, 'g^', mfc='none', label = "TemplateDumerV1 Samples")
	# plt.plot(x, interpolated3, 'g--', mfc='none', label = "TemplateDumerV1 Interpolation")
	# plt.xticks(x, rotation = 70)

	# plt.xlabel(r"$n$")
	# plt.ylabel(r'Estimated runtime log$ _2$ s')
	# #plt.legend(loc="upper left")
	# plt.legend()

	# plt.savefig("times.pdf", format="pdf", bbox_inches="tight")
	# plt.show()