import os, time, math, pickle, subprocess, sys, decimal
import numpy as np
basedir = os.getcwd().replace("\\", "/") + "/"
ms0_path = basedir + "output/BGR.ms0.csv"
ms1_path = basedir + "output/BGR.ms1.csv"
ma0_path = basedir + "output/BGR.ma0.csv"
libpath = basedir + "result/library.csv"
respath = basedir + "result/result.csv"
sumpath = basedir + "result/summary.csv"
def log_write(string):
	output_str= str(int(time.time())) + ",;\t,;" + string
	print(output_str[:output_str.index(",")] + ":\t" + string)
	logfile = open(basedir + "run.log", "at")
	logfile.write("\n" + output_str)
	logfile.close()
def sub_call(string):
	log_write("calling \'" + str(string) + "\' on PowerShell 7")
	subprocess.call(str(string), shell=True)
def sub_popen(string, popen_in=None):
	log_write("calling \'" + str(string) + "\' on PowerShell 7, returning output to program")
	if popen_in is not None:
		popen_in = str(popen_in)
		log_write("\'"+popen_in+"\' will be automatically inputted in the program call above")
	return subprocess.Popen(str(string), shell=True, stdin=popen_in, stdout=subprocess.PIPE).communicate()[0]
def avg(inputlist):
	return sum(inputlist) / float(len(inputlist))
def std(inputlist, avg_init=None):
	if avg_init is None:
		avg_ = avg(inputlist)
	else:
		avg_ = avg_init
	return pow(avg([pow(float(x - avg_), 2) for x in inputlist]), 0.5)
def measform3(x):
	if isinstance(x, float):
		if x < 0:
			return "-" + measform3(-x)
		elif x == 0:
			return "0e0"
		else:
			exponent = int(math.floor(math.log10(x)))
			return str((decimal.Decimal(x) * decimal.Decimal(10 ** (0 - exponent)), decimal.Decimal(x) / decimal.Decimal(10 ** exponent))[exponent > 0]) + "e" + str(exponent)
	else:
		return str(x)
def stats(inputlist):
	return (lambda x: x + [std(inputlist, x[-1])])([min(inputlist), max(inputlist), avg(inputlist)])
def extract(output_file, col_list, monte=False):
	col_index = {}
	content = {}
	passed_title = False
	with open(output_file, "rt") as output_csv:
		for line in output_csv:
			splitline = line.rstrip().split(",")
			if "alter#" in line:
				passed_title = True
				for col in col_list + ["alter#"]:
					col_index[col] = splitline.index(col)
			elif passed_title:
				try:
					alter_index = int(splitline[col_index["alter#"]])
					if alter_index in content:
						for col in col_list:
							try:
								value = float(splitline[col_index[col]])
								# value = int(round(value * 2 ** 27)) / 1024.0
								value = (value, int(value))[value == int(value)]
								if monte:
									content[alter_index][col].append(value)
								else:
									content[alter_index][col] = value
							except Exception:
								log_write(str(sys.exc_info()[0]) + " @ extract, existing alter : " + str(alter_index) + " > " + str(sys.exc_info()[2]))
					else:
						content[alter_index] = {}
						for col in col_list:
							try:
								value = float(splitline[col_index[col]])
								# value = int(round(value * 1024)) / 1024.0
								value = (value, int(value))[value == int(value)]
								content[alter_index][col] = (value, [value])[monte]
							except Exception:
								log_write(str(sys.exc_info()[0]) + " @ extract, new alter : " + str(alter_index) + " > " + str(sys.exc_info()[2]))
								content[alter_index][col] = (None, [])[monte]
					if False in [bool(content[alter_index][col]) or content[alter_index][col] == 0 for col in col_list]:
						content.pop(alter_index)
				except ValueError:
					log_write(str(sys.exc_info()[0]) + " @ extract, alter determination : " + str(alter_index) + " > " + str(line.rstrip().lstrip()))
	# if monte:	content[alter_value][col in col_list] = [monte0, monte1, ..., monte_last]
	# else:		content[alter_value][col in col_list] = monte_last
	# alter_value starts from 1 to len(content)
	return content

# step 0: making step-by-step flowchart (this)
# step 1: extract
# 1.1. 		Extract *.ms0.csv
# 1.1.1.		Obtain index of 'creation_time', 'alter#', and everything on var_list on the line that has 'alter#' in it
# 1.1.2.		Obtain index of everything in outvar_list_dc on the same line mentioned in 1.1.1
# 1.1.3.		Create a (alter#, creation_time) dictionary for both library and result
# 1.1.4.		Create a (alter#, (var_list items, var_list value)) dictionary for library
# 1.1.5.		Create a (alter#, (outvar_list_dc items, outvar_list_dc value)) for result
# 1.2.		Extract *.ma0.csv
# 1.2.1.		Obtain index of 'alter#' and everything on outvar_list_ac on the line that has 'alter#' in it
# 1.2.2.		Create a (alter#, (outvar_list_ac items, outvar_list_ac value)) for result
# 1.3.		Extract *.ms1.csv
# 1.3.1.		Obtain index of 'alter#' and everything on outvar_list_ls on the line that has 'alter#' in it
# 1.3.2.		Create a (alter#, (outvar_list_ls items, outvar_list_ls value)) for result
# step 2: analysis
# 2.1.		Meta analysis and processing
# 2.1.1.		Calculate FoM and create (alter#, FoM) dictionary using values obtained in step 1
# 2.1.2.		Generate and pickle [(FoM, (hyperparameter dict))] parents by comparing FoM of each alter and adding its respective hyperparameters
# 2.1.3.		Calculate secondary FoM
# 2.1.4.		Calculate FoM from prior art: https://www.mdpi.com/2072-666X/14/7/1420
# step 3: documentation
# 3.0.		Record current time
# 3.1.		Add time from 3.0, creation_time, and hyperparameters to library
# 3.2.		Add time from 3.0, creation_time, FoM(s), and everything on outvar_list_* to result
decision_value = 0
if len(sys.argv) > 1:
	try:
		decision_value = int(float(sys.argv[1]) / 20) # switches FoM every 20 cycles
		decision_value = (0, decision_value)[decision_value > 0] % 2
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ decision_value : " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		decision_value = 0
outvar_list_dc = ["stage0_vref_maxat", "stage0_vref_minat", "stage0_vref_tc", "stage0_out_maxat", "stage0_out_minat", "stage0_out_tc","output_avg",'output_max','output_maxat','output_min','output_minat','output_tc','output_curve','output_mintemp','output_maxtemp','output_vout0','output_temp0','output_vout1','output_temp1','output_vout2','output_temp2','power_all']
outvar_list_ac = ["vref_psrr", "vref_khz_psrr_max", "vref_khz_max_at", "vref_mhz_psrr_max", "vref_mhz_max_at"]
var_list = ["vref_op_pch3_w", "vref_op_pch3_l", "vref_op_pch3_m", "vref_op_nch_w", "vref_op_nch_l", "vref_op_nch_m", "vref_op_stg2_w", "vref_op_stg2_m", "vref_op_stg2_l", "vref_op_join_l", "vref_op_outn_l", "vref_op_out3_l", "stage0_pch_w", "stage0_pch_l", "stage0_pch_m", "op1_stack0_w", "op1_stack0_l", "op1_stack0_m", "op1_stack1_w", "op1_stack1_l", "op1_stack1_m", "op1_join_w", "op1_join_l", "op1_join_m", "op2_input_w", "op2_input_l", "op2_input_m", "op2_base_w", "op2_base_l", "op2_base_m", "res_rr0_w", "res_rr0_l", "res_rr0_m", "res_rr1_w", "res_rr1_l", "res_rr1_m", "res_rr3_w", "res_rr3_l", "res_rr3_m", "res_rr5_w", "res_rr5_l", "res_rr5_m", "psrr_stabilizer_w", "psrr_stabilizer_l", "psrr_stabilizer_m"]
outvar_list_ls = ['ls_vref_max','ls_vsrc30_maxat','ls_vsrc18_maxat','ls_vref_min','ls_vsrc30_minat','ls_vsrc18_minat','ls_vref_avg']
ms0_csv = extract(ms0_path, outvar_list_dc, True)
ma0_csv = extract(ma0_path, outvar_list_ac, True)
ms1_csv = extract(ms1_path, outvar_list_ls + var_list + ["creation_time"], True)

ms0_tc_stats = {}
# TC stats holds min, max, avg, std for TC of stage 0's voltage reference, stage 0's output, and final output, in that order
ms0_tc_linefit = {}
# linefit holds c, m, R^2 for line fitting out0->out1, out1->out2, and out0->out2, in that order
# both tc_stats and tc_linefit takes polarity into account (argmax < argmin -> TC < 0)
ms0_out_stats = {}
# out_stats holds min, max, avg, std for average value of final output
ms0_pow_stats = {}
# pow_stats holds min, max, avg, std for power consumption
ms0_out_curve = {}
# out_curve holds 5-7 (temp, vout) points of interest
ms0_out_curveavg = {}
for alter in ms0_csv:
	ms0_tc_stats[alter] = []
	out0 = [x * (1, -1)[ms0_csv[alter]["stage0_vref_maxat"][z] < ms0_csv[alter]["stage0_vref_minat"][z]] for z, x in enumerate(ms0_csv[alter]["stage0_vref_tc"])]
	out1 = [x * (1, -1)[ms0_csv[alter]["stage0_out_maxat"][z] < ms0_csv[alter]["stage0_out_minat"][z]] for z, x in enumerate(ms0_csv[alter]["stage0_out_tc"])]
	out2 = [x * (1, -1)[ms0_csv[alter]["output_maxat"][z] < ms0_csv[alter]["output_minat"][z]] for z, x in enumerate(ms0_csv[alter]["output_tc"])]
	# stats = lambda list0: (lambda x: x + [std(list0, x[-1])])([min(list0), max(list0), avg(list0)])
	ms0_tc_stats[alter] = ms0_tc_stats[alter] + stats(out0)
	ms0_tc_stats[alter] = ms0_tc_stats[alter] + stats(out1)
	ms0_tc_stats[alter] = ms0_tc_stats[alter] + stats(out2)[0:2]
	ms0_tc_stats[alter] = ms0_tc_stats[alter] + stats(ms0_csv[alter]["output_tc"])[2:4]
	ms0_tc_linefit[alter] = []
	ms0_out_stats[alter] = stats(ms0_csv[alter]["output_avg"])
	ms0_pow_stats[alter] = stats(ms0_csv[alter]["power_all"])
	linefit = lambda list0, list1: (lambda x: list(x[0]) + [1 - float(x[1][0]) / (len(list1) * math.pow(std(list1), 2))])(np.polynomial.polynomial.polyfit(list0, list1, 1, full=True))
	ms0_tc_linefit[alter] = ms0_tc_linefit[alter] + linefit(out0, out1)
	ms0_tc_linefit[alter] = ms0_tc_linefit[alter] + linefit(out1, out2)
	ms0_tc_linefit[alter] = ms0_tc_linefit[alter] + linefit(out0, out2)
	ms0_out_curve[alter] = []
	# for z, offset in enumerate(ms0_csv[alter]["output_curve"]):
		# points = []
		# min_point = (ms0_csv[alter]["output_minat"][z], ms0_csv[alter]["output_min"][z])
		# max_point = (ms0_csv[alter]["output_maxat"][z], ms0_csv[alter]["output_max"][z])
		# poi = [(ms0_csv[alter]["output_temp" + str(z0)][z], ms0_csv[alter]["output_vout" + str(z0)][z]) for z0 in range(3)]
		# points.append((-40, ms0_csv[alter]["output_mintemp"][z]))
		# points.append((125, ms0_csv[alter]["output_maxtemp"][z]))
		# points.append(min_point)
		# points.append(max_point)
		# points = sorted(list(set(points + poi)), key:lambda x: x[0])
		# minmax_line = [min_point[1] - max_point[1], max_point[0] - min_point[0]]
		# if minmax_line[1] < 0:
			# minmax_line = [-x for x in minmax_line]
		# minmax_line.append(0 - minmax_line[0] * min_point[0] - minmax_line[1] * min_point[1])
		# poi_offset = [(minmax_line[0] * point[0] + minmax_line[1] * point[1] + minmax_line[2]) / (2 * minmax_line[1]) for point in poi]
			# # |poi_offset * minmax_range| < |area of curve with monotonical slope above the line|
			# # offset * minmax_range = area of actual curve above the line, indicates how concave/convex is the overall curve
		# inflection_in_minmax_curve = True in [halfheight * offset < 0 or abs(halfheight) > abs(offset) for halfheight in poi_offset]
		# # if any of |poi_offset| > |offset| or poi_offset * offset < 0 is true, there must exist a convex and concave part
		# # all of poi_offset <= offset doesn't mean that all of it is convex/concave, it just means no meaningful convex/concave part
		# ms0_out_curve[alter].append([offset, inflection_in_minmax_curve, points])
	ms0_out_curve[alter].append((-40, avg(ms0_csv[alter]["output_mintemp"])))
	ms0_out_curve[alter].append((125, avg(ms0_csv[alter]["output_maxtemp"])))
	ms0_out_curve[alter].append((avg(ms0_csv[alter]["output_maxat"]), avg(ms0_csv[alter]["output_max"])))
	ms0_out_curve[alter].append((avg(ms0_csv[alter]["output_minat"]), avg(ms0_csv[alter]["output_min"])))
	temp_line = np.polynomial.polynomial.polyfit(ms0_csv[alter]["output_temp0"], ms0_csv[alter]["output_vout0"], 1)
	vout_avg = avg(ms0_csv[alter]["output_vout0"])
	ms0_out_curve[alter].append(((vout_avg - temp_line[0])/temp_line[1], vout_avg))
	#ms0_out_curve[alter].append((avg(ms0_csv[alter]["output_temp0"]), vout_avg))
	temp_line = np.polynomial.polynomial.polyfit(ms0_csv[alter]["output_temp1"], ms0_csv[alter]["output_vout1"], 1)
	vout_avg = avg(ms0_csv[alter]["output_vout1"])
	ms0_out_curve[alter].append(((vout_avg - temp_line[0])/temp_line[1], vout_avg))
	#ms0_out_curve[alter].append((avg(ms0_csv[alter]["output_temp1"]), vout_avg))
	temp_line = np.polynomial.polynomial.polyfit(ms0_csv[alter]["output_temp2"], ms0_csv[alter]["output_vout2"], 1)
	vout_avg = avg(ms0_csv[alter]["output_vout2"])
	ms0_out_curve[alter].append(((vout_avg - temp_line[0])/temp_line[1], vout_avg))
	#ms0_out_curve[alter].append((avg(ms0_csv[alter]["output_temp2"]), vout_avg))
	ms0_out_curve[alter] = sorted([point for point in set(ms0_out_curve[alter]) if point[0] >= -40 and point[0] <= 125], key=lambda x: x[0])
	ms0_out_curveavg[alter] = stats(ms0_csv[alter]["output_curve"])

ma0_psrr_stats = dict([(alter, stats(ma0_csv[alter]["vref_psrr"]))for alter in ma0_csv])
# psrr_stats holds min, max, avg, std for power supply rejection ratio (AC gain from 1.8 V voltage source) @ 100 Hz

ms1_ls_stats = {}
for alter in ms1_csv:
	current_list = []
	ext = min(ms1_csv[alter]["ls_vref_min"])
	ext_30, ext_18 = (lambda z: (ms1_csv[alter]["ls_vsrc30_minat"][z], ms1_csv[alter]["ls_vsrc18_minat"][z]))(ms1_csv[alter]["ls_vref_min"].index(ext))
	current_list = current_list + [ext, ext_30, ext_18]
	ext = max(ms1_csv[alter]["ls_vref_max"])
	ext_30, ext_18 = (lambda z: (ms1_csv[alter]["ls_vsrc30_maxat"][z], ms1_csv[alter]["ls_vsrc18_maxat"][z]))(ms1_csv[alter]["ls_vref_max"].index(ext))
	current_list = current_list + [ext, ext_30, ext_18]
	average = avg(ms1_csv[alter]["ls_vref_avg"])
	line_swing = abs(current_list[3] - current_list[0]) / average
	ms1_ls_stats[alter] = current_list + [average, line_swing]
# ls_stats holds min, VDD_30 @ min, VDD_18 @ min, max, VDD_30 @ max, VDD_18 @ max, avg, and (max - min)/avg of final output, sweeping 1.8 <= VDD_18 <= 5 and 2.6 <= VDD_30 <= 5

fom = {}
for alter in ms1_csv:
	try:
		temperature_coefficient = max(abs(ms0_tc_stats[alter][8]), abs(ms0_tc_stats[alter][9]))
		current_fom = (0 - ma0_psrr_stats[alter][2]) * 165 * 3.0 / (temperature_coefficient * ms0_pow_stats[alter][2])
		# log_write("DEBUG: " + str(int(ms1_csv[alter]["creation_time"][-1])) + ".sp > " + str((temperature_coefficient, current_fom)))
		# FoM from equation (9)
		current_fom = (1 - math.pow(1 + abs(ms0_out_stats[alter][2] - 3), -1)) * (1 - math.pow(1 + abs(ms0_out_stats[alter][2] - 1.8), -1)) * current_fom / (ms1_ls_stats[alter][-1] * temperature_coefficient)
		fom[alter] = current_fom
	except ZeroDivisionError:
		log_write(str(sys.exc_info()[0]) + " @ FoM calculation for sim_"+ str(int(ms1_csv[alter]["creation_time"][-1]))+ ".sp: " + ''.join([str(x) + ", " for x in [ms0_tc_stats[alter][8], ms0_tc_stats[alter][9], ms0_pow_stats[alter][2], ms1_ls_stats[alter][-1]]])[:-2])
		#log_write("DEBUG: " + str(ms0_pow_stats[alter]))
		#log_write("DEBUG: " + str(ms0_csv[alter]["power_all"]))
	except KeyError:
		log_write(str(sys.exc_info()[0]) + " @ FoM calculation for sim_"+ str(int(ms1_csv[alter]["creation_time"][-1]))+ ".sp: " + str(alter) + " > " + ''.join([str(alter in x) + ", " for x in [ma0_psrr_stats, ms0_tc_stats, ms0_pow_stats, ms0_out_stats, ms1_ls_stats]])[:-2])
parents = [(fom[alter], dict([(key, float(ms1_csv[alter][key][-1])) for key in var_list] + [("title", "sim_"+str(int(ms1_csv[alter]["creation_time"][-1]))+".sp")])) for alter in fom if int(ms1_csv[alter]["creation_time"][-1]) != 1e7]
if len(parents) > 0:
	with open(basedir + "selectedparents.pickle", "wb") as parent_pickle:
		pickle.dump(parents, parent_pickle, 0)

logtime = int(time.time() * 1e6)
try:
	open(respath, "rt").close()
except OSError:
	try:
		with open(respath, "wt") as rescsv:
			title_list = ["time","identifier"]
			title_list.append("fom - pso")
			title_list.append("fom - (9)")
			title_list.append("sep_0")
			title_list = title_list + ["final output - vout(T, 3.0, 1.8) - " + x for x in ["min", "max", "avg", "std"]]
			title_list = title_list + ["power consumption - " + x for x in ["min", "max", "avg", "std"]]
			title_list.append("sep_1")
			title_list = title_list + ["vout(T, 3.0, 1.8) TC - " + x for x in ["min", "max", "avg", "std"]]
			title_list = title_list + ["stage 0 output TC - " + x for x in ["min", "max", "avg", "std"]]
			title_list = title_list + ["stage 0 vref TC - " + x for x in ["min", "max", "avg", "std"]]
			title_list.append("sep_2")
			title_list = title_list + ["linefit vref-output TC - " + x for x in ["const", "slope", "R^2"]]
			title_list = title_list + ["linefit output-vout TC - " + x for x in ["const", "slope", "R^2"]]
			title_list = title_list + ["linefit vref-vout TC - " + x for x in ["const", "slope", "R^2"]]
			title_list.append("sep_3")
			title_list = title_list + ["vout(25, VDD_18, VDD_30) min", "VDD_30 @ vout min", "VDD_18 @ vout min"]
			title_list = title_list + ["vout(25, VDD_18, VDD_30) max", "VDD_30 @ vout max", "VDD_18 @ vout max"]
			title_list = title_list + ["vout(25, VDD_18, VDD_30) avg", "line sensitivity"]
			title_list = title_list + ["1.8 V PSRR (AC gain @ 100 Hz) - " + x for x in ["min", "max", "avg", "std"]]
			title_list.append("sep_4")
			title_list = title_list + ["vout(T, 3.0, 1.8) curvature identifier (avg - (min + max) / 2) - " + x for x in ["min", "max", "avg", "std"]]
			title_list.append("sep_5")
			title_list = title_list + ["point " + str(z) + " - temp, point " + str(z) + " - output" for z in range(7)]
			rescsv.write(''.join([str(x) + ", " for x in title_list])[:-2] + "\n")
	except OSError:
		log_write(str(sys.exc_info()[0]) + " @ CSV result init write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
try:
	with open(respath, "at") as rescsv:
		for alter in sorted(list(fom), key=lambda x: fom[x], reverse=True):
			if int(ms1_csv[alter]["creation_time"][-1]) != 1e7:
				use_measform3 = [2, 3, 8, 9, 10, 11, 12] + [z + 14 for z in range(12)] + [z + 27 for z in range(9) if z % 3 != 2] + [z + 50 for z in range(4)]
				try:
					content_list = [logtime, int(ms1_csv[alter]["creation_time"][-1])]
					content_list.append(fom[alter])
					content_list.append((0 - ma0_psrr_stats[alter][2]) * 165 * 3.0 / (temperature_coefficient * ms0_pow_stats[alter][2]))
					content_list.append("000000")	# zero_sep_0, 4
					content_list = content_list + ms0_out_stats[alter]
					content_list = content_list + ms0_pow_stats[alter]
					content_list.append(111111)		# zero_sep_1, 13
					content_list = content_list + ms0_tc_stats[alter][8:12]
					content_list = content_list + ms0_tc_stats[alter][4:8]
					content_list = content_list + ms0_tc_stats[alter][0:4]
					content_list.append(222222)		# zero_sep_2, 26
					content_list = content_list + ms0_tc_linefit[alter]
					content_list.append(333333)		# zero_sep_3, 36
					content_list = content_list + ms1_ls_stats[alter]
					content_list = content_list + ma0_psrr_stats[alter]
					content_list.append(444444)		# zero_sep_4, 49
					content_list = content_list + ms0_out_curveavg[alter]
					content_list.append(555555)		# zero_sep_5, 54
					for points in ms0_out_curve[alter]:
						content_list = content_list + list(points)
					rescsv.write(''.join([(str(x), measform3(x))[z in use_measform3] + ", " for z, x in enumerate(content_list)])[:-2] + "\n")
				except KeyError:
					log_write(str(sys.exc_info()[0]) + " @ CSV result content write: sim_" + str(int(ms1_csv[alter]["creation_time"][-1])) + ".sp > " + str(len(content_list)))
except OSError:
	log_write(str(sys.exc_info()[0]) + " @ CSV result content write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
	raise
	
try:
	open(libpath, "rt").close()
except OSError:
	try:
		with open(libpath, "wt") as libcsv:
			libcsv.write(''.join([str(x) + ", " for x in ["time","identifier"] + var_list])[:-2] + "\n")
	except OSError:
		log_write(str(sys.exc_info()[0]) + " @ CSV library init write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
try:
	with open(libpath, "at") as libcsv:
		for alter in sorted(list(fom), key=lambda x: fom[x], reverse=True):
			if int(ms1_csv[alter]["creation_time"][-1]) != 1e7:
				libcsv.write(''.join([str(x) + ", " for x in [logtime, int(ms1_csv[alter]["creation_time"][-1])] + [int(round(ms1_csv[alter][col][-1] * 1000)) / 1000.0 for col in var_list]])[:-2] + "\n")
except OSError:
	log_write(str(sys.exc_info()[0]) + " @ CSV library content write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
	raise

ctime_index = None
fom_index = None
time_index = None
past_time = []
text_dict = {}
first_line = ""
try:
	with open(respath, "rt") as rescsv:
		for line in rescsv:
			splitline = line.rstrip().split(", ")
			if None in [time_index, ctime_index, fom_index]:
				time_index = splitline.index("time")
				ctime_index = splitline.index("identifier")
				fom_index = splitline.index("fom - pso")
				first_line = line.rstrip()
			else:
				past_time.append(int(splitline[time_index]))
				creation_time = int(splitline[ctime_index])
				if (text_dict[creation_time][0] if creation_time in text_dict else -float("inf")) < float(splitline[fom_index]):
					text_dict[creation_time] = (float(splitline[fom_index]), line.rstrip())
	duration = (lambda x: [float(x[z + 1] - x[z]) / 1e6 for z in range(len(x) - 1)])(sorted(list(set(past_time))))
	with open(sumpath, "wt") as sumcsv:
		sumcsv.write(first_line + "\n")
		for creation_time in sorted(text_dict, key=lambda x: text_dict[x][0], reverse=True):
			sumcsv.write(text_dict[creation_time][1] + "\n")
		sumcsv.write("\ncycle time - min, cycle time - max, cycle time - avg, cycle time - std\n")
		stats_duration = stats(duration)
		sumcsv.write(''.join([str(x) + ", " for x in stats_duration])[:-2] + "\n")
		stats_duration = stats([x for x in duration if abs((stats_duration[2] - x) / stats_duration[3]) < 2])
		sumcsv.write(''.join([str(x) + ", " for x in stats_duration])[:-2] + "\n")
except OSError:
	log_write(str(sys.exc_info()[0]) + " @ CSV result content summary: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
	raise