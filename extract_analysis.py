import os, time, math, pickle, subprocess, sys
basedir = os.getcwd() + "/"
def log_write(string):
	output_str= str(int(time.time())) + ",;\t,;" + string
	print(output_str[:output_str.index(",")] + ":\t" + string)
	logfile = open(basedir + "run.log", "at")
	logfile.write("\n" + output_str)
	logfile.close()
def sub_call(string):
	log_write("calling \'" + str(string) + "\' on shell")
	subprocess.call(str(string), shell=True)
def sub_popen(string, popen_in=None):
	log_write("calling \'" + str(string) + "\' on shell, returning output to program")
	if popen_in is not None:
		popen_in = str(popen_in)
		log_write("\'"+popen_in+"\' will be automatically inputted in the program call above")
	return subprocess.Popen(str(string), shell=True, stdin=popen_in, stdout=subprocess.PIPE).communicate()[0]
def avg(inputlist):
	return sum(inputlist) / float(len(inputlist))
def std(inputlist):
	avg_ = avg(inputlist)
	return pow(avg([pow(float(x - avg_), 2) for x in inputlist]), 0.5)
def measform3(x):
	return str(x * pow(10, -int(math.floor(math.log10(x)))))+"e"+str(int(math.floor(math.log10(x))))

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
# 1.4.		Extract *.mt0.csv
# 1.4.1.		Obtain index of 'alter#' and everything on outvar_list_tr on the line that has 'alter#' in it
# 1.4.2.		Create a (alter#, (outvar_list_ls items, [outvar_list_ls values])), note, plural values (Monte results)
# 1.4.3.		Create a (alter#, ([vref_avg_proc, vref_std_proc, pow_avg_proc, pow_std_proc])) for result
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
var_list = ["op1_stage1_pmos_w","op1_stage1_pmos_l","op1_stage1_pmos_m","op1_stage1_nmos_w","op1_stage1_nmos_l","op1_stage1_nmos_m","op1_stage2_nmos_w","op1_stage2_nmos_l","op1_stage2_nmos_m","op1_vref_nmos_w","op1_vref_nmos_l","op1_vref_nmos_m","op1_vref_pmos_w","op1_vref_pmos_l","op1_vref_pmos_m","op1_stage2_pmos_w","op1_stage2_pmos_l","op1_stage2_pmos_m","psrr_stabilizer_w","psrr_stabilizer_l","psrr_stabilizer_m","res_rr0_w","res_rr0_l","res_rr0_m","res_rr2_w","res_rr2_l","res_rr2_m"]
outvar_list_dc = ["vref_tc", "temp_avg_vref", "temp_vref_max", "temp_max_vref", "temp_vref_min", "temp_min_vref"]
outvar_list_ac = ["vref_psrr", "vref_khz_psrr_max", "vref_khz_max_at", "vref_mhz_psrr_max", "vref_mhz_max_at"]
outvar_list_ls = ["ls_vref_avg","ls"]
outvar_list_tr = ["monte_avg_vref", "monte_avg_pow"]
if True:
	# step 1
	sub_call("".join(["cd ",basedir,"; mkdir temp"]))
	sub_call("".join(["cd ",basedir,"; scp selectedparents.pickle ./temp/"]))
	sub_call("".join(["cd ",basedir,"output/; scp *.csv ",basedir,"temp/"]))
	
	# step 1.1
	dict_113 = {}
	dict_114 = {}
	dict_115 = {}
	try:
		inp_ind = {}
		odc_ind = {}
		alt_ind = 0
		crt_ind = 0
		src = open(basedir+"temp/BGRN.ms0.csv", "rt")
		for line in src:
			splitline = [x.rstrip().lstrip() for x in line.split(",")]
			if len(inp_ind) == 0 and len(odc_ind) == 0:
				if "alter#" in splitline:
					# step 1.1.1
					alt_ind = splitline.index("alter#")
					crt_ind = splitline.index("creation_time")
					for key in var_list:
						inp_ind[key] = splitline.index(key)
					# step 1.1.2
					for key in outvar_list_dc:
						odc_ind[key] = splitline.index(key)	
			else:
				creation_time = int(float(splitline[crt_ind]))
				alter_11 = int(float(splitline[alt_ind]))
				if not creation_time == 1e7 or True:
					dict_113[alter_11] = creation_time
					dict_114[alter_11] = dict([(key, float(splitline[inp_ind[key]])) for key in var_list])
					dict_115[alter_11] = dict([(key, float(splitline[odc_ind[key]])) for key in outvar_list_dc])
		src.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ read *.ms0.csv: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	
	# step 1.2
	dict_122 = {}
	try:
		oac_ind = {}
		alt_ind = 0
		src = open(basedir+"temp/BGRN.ma0.csv", "rt")
		for line in src:
			splitline = [x.rstrip().lstrip() for x in line.split(",")]
			if len(oac_ind) == 0:
				if "alter#" in splitline:
					alt_ind = splitline.index("alter#")
					for key in outvar_list_ac:
						oac_ind[key] = splitline.index(key)	
			else:
				alter_12 = int(float(splitline[alt_ind]))
				if True:
					dict_122[alter_12] = dict([(key, float(splitline[oac_ind[key]])) for key in outvar_list_ac])
		src.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ read *.ma0.csv: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	
	# step 1.3
	dict_132 = {}
	try:
		ols_ind = {}
		alt_ind = 0
		src = open(basedir+"temp/BGRN.ms1.csv", "rt")
		for line in src:
			splitline = [x.rstrip().lstrip() for x in line.split(",")]
			if len(ols_ind) == 0:
				if "alter#" in splitline:
					alt_ind = splitline.index("alter#")
					for key in outvar_list_ls:
						ols_ind[key] = splitline.index(key)	
			else:
				alter_13 = int(float(splitline[alt_ind]))
				if True:
					dict_132[alter_13] = dict([(key, float(splitline[ols_ind[key]])) for key in outvar_list_ls])
		src.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ read *.ms1.csv: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise

	# step 1.4
	dict_142 = {}
	dict_143 = {}
	try:
		otr_ind = {}
		alt_ind = 0
		src = open(basedir+"temp/BGRN.mt0.csv", "rt")
		for line in src:
			splitline = [x.rstrip().lstrip() for x in line.split(",")]
			if len(otr_ind) == 0:
				if "alter#" in splitline:
					alt_ind = splitline.index("alter#")
					for key in outvar_list_tr:
						otr_ind[key] = splitline.index(key)	
			else:
				alter_14 = int(float(splitline[alt_ind]))
				if not alter_14 in dict_142:
					dict_142[alter_14] = {}
					for key in outvar_list_tr:
						dict_142[alter_14][key] = [float(splitline[otr_ind[key]])]
				else:
					for key in outvar_list_tr:
						dict_142[alter_14][key].append(float(splitline[otr_ind[key]]))
		for alter in dict_142:
			dict_143[alter] = {}
			dict_143[alter]["vref_avg_proc"] = avg(dict_142[alter]["monte_avg_vref"])
			dict_143[alter]["vref_std_proc"] = std(dict_142[alter]["monte_avg_vref"])
			dict_143[alter]["pow_avg_proc"] = avg(dict_142[alter]["monte_avg_pow"])
			dict_143[alter]["pow_std_proc"] = std(dict_142[alter]["monte_avg_pow"])
		src.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ read *.mt0.csv: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
		
	# step 2
	dict_211 = {}
	dict_213 = {}
	dict_214 = {}
	# step 2.1
	for alter in dict_113:
		try:
			psrr_factor = 0 - dict_122[alter]["vref_psrr"]
			proc_factor = 10 * math.log10(dict_143[alter]["vref_avg_proc"] / dict_143[alter]["vref_std_proc"])
			tc_factor = 10 * math.log10(165.0 / dict_115[alter]["vref_tc"])
			ls_factor = 0 - 10 * math.log10(dict_132[alter]["ls"])
			stability_factor = 0 + 2.5 * math.log10(3.0/(pow(dict_143[alter]["vref_avg_proc"] - dict_115[alter]["temp_avg_vref"], 4) + pow(dict_115[alter]["temp_avg_vref"] - dict_132[alter]["ls_vref_avg"], 4) + pow(dict_132[alter]["ls_vref_avg"] - dict_143[alter]["vref_avg_proc"], 4)))
			# primary FoM is decided by Joshua Adiel Wijaya to be -psrr * 10 * log10(v_avg/v_std * (temp_range)/tc * 1/(ls) * 1/(4th power avg of delta V)), no known citation
			# note, power consumption do not factor in to FoM
			dict_211[alter] = proc_factor + tc_factor + ls_factor + stability_factor
			# alternate FoM is decided by Joshua Adiel Wijaya to be e^(-4/3 * psrr/tc - 1), exaggerates circuits with TC < 80 ppm/C
			dict_213[alter] = math.exp((4.0/3.0) * psrr_factor / dict_115[alter]["vref_tc"] - 1)
			# FoM from prior art is |PSRR| * temp_range * VDD / (TC * power)
			dict_214[alter] = psrr_factor * 165 * 1.8 / (dict_115[alter]["vref_tc"] * dict_143[alter]["pow_avg_proc"])
		except Exception:
			log_write("Step 2 FoM calculation error - skipping " + str(dict_113[alter]))
			log_write(str(sys.exc_info()[0]) + " @ FoM calc: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
	# step 2.2
	try:
		fom_couple_list = []		
		for alter in dict_113:
			used_fom = (dict_213[alter], dict_213[alter])[decision_value]
			if dict_113[alter] != int(1e7):
				fom_couple_list.append((used_fom, dict([(key, dict_114[alter][key]) for key in var_list] + [("title", "sim_"+str(dict_113[alter])+".sp")])))
		output_dest = open(basedir + "selectedparents.pickle", "wb")
		pickle.dump(fom_couple_list, output_dest, 0)
		output_dest.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ parent pickling: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	
	# step 3
	sub_call("".join(["cd ", basedir,"; mkdir result"]))
	var_list_min1 = list(set([key[:-1] for key in var_list if key[:-1].endswith("_")]))
	for alter in dict_114:
		for partialkey in var_list_min1:
			if False and partialkey + "w" in dict_114[alter] and partialkey + "m" in dict_114[alter]:
				conv_m = 2 ** int(round(math.log(float(dict_114[alter][partialkey + "m"]), 2)))
				conv_w = int(float(dict_114[alter][partialkey + "w"]) * float(dict_114[alter][partialkey + "m"]) / (conv_m * 0.025)) * 0.025
				dict_114[alter][partialkey + "m"] = str(conv_m)
				dict_114[alter][partialkey + "w"] = str(conv_w)
			if False and partialkey + "width" in dict_114[alter] and partialkey + "ratio" in dict_114[alter]:
				conv_l = int(float(dict_114[alter][partialkey + "width"]) * float(dict_114[alter][partialkey + "ratio"]) / 0.025) * 0.025
				dict_114[alter][partialkey + "w"] = dict_114[alter][partialkey + "width"]
				dict_114[alter][partialkey + "l"] = str(conv_l)
	var_list = ["op1_stage1_pmos_w","op1_stage1_pmos_l","op1_stage1_pmos_m","op1_stage1_nmos_w","op1_stage1_nmos_l","op1_stage1_nmos_m","op1_stage2_nmos_w","op1_stage2_nmos_l","op1_stage2_nmos_m","op1_vref_nmos_w","op1_vref_nmos_l","op1_vref_nmos_m","op1_vref_pmos_w","op1_vref_pmos_l","op1_vref_pmos_m","op1_stage2_pmos_w","op1_stage2_pmos_l","op1_stage2_pmos_m","psrr_stabilizer_w","psrr_stabilizer_l","psrr_stabilizer_m","res_rr0_w","res_rr0_l","res_rr0_m","res_rr2_w","res_rr2_l","res_rr2_m"]
	include_original = False
	try:
		open(basedir + "result/library.csv","rt").close()
		log_write("Found " + basedir + "result/library.csv")
	except Exception:
		log_write("First creation of " + basedir + "result/library.csv")
		try:
			initial = open(basedir + "result/library.csv","wt")
			initial.write("time,identifier,"+"".join([x + "," for x in var_list])+"\n")
			initial.close()
		except Exception:
			log_write(str(sys.exc_info()[0]) + " @ CSV library init write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
			raise
	try:
		open(basedir + "result/result.csv","rt").close()
		log_write("Found " + basedir + "result/result.csv")
	except Exception:
		log_write("First creation of " + basedir + "result/result.csv")
		try:
			initial = open(basedir + "result/result.csv","wt")
			initial.write("time,identifier,FoM_0,FoM_1,FoM_2,decision_value,DC - temperature coefficient (ppm/C),DC - average vref (temp -40C to 125C),DC - max vref,DC - temp @ max vref,DC - min vref,DC - temp @ min vref,AC - PSRR @ 100 Hz (dB),max PSRR @ kHz range,frequency @ kHz PSRR max,max PSRR @ MHz range,frequency @ MHz PSRR max,LS - average vref (vdd 1.8V to 3.5V),LS - line sensitivity,MC - avg/std vref,MC - average vref,MC - stdev vref,MC - average power,MC - stdev power"+"\n")
			initial.close()
			include_original = True
		except Exception:
			log_write(str(sys.exc_info()[0]) + " @ CSV result init write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
			raise
	# step 3.0
	logtime = int(time.time() * 1e6)
	dict_113_sorted = sorted(list(dict_113), key=lambda x: dict_213[x], reverse=True)
	if not include_original:
		dict_113_sorted = [alter for alter in dict_113_sorted if dict_113[alter] != 10000000]
	# step 3.1
	try:
		log_write("Writing to " + basedir + "result/library.csv")
		csv_file = open(basedir + "result/library.csv","at")
		for line in ["".join([str(x) + ", " for x in [logtime, dict_113[alter]] + [dict_114[alter][key] for key in var_list]]) for alter in dict_113_sorted]:
			csv_file.write(line + "\n")
		csv_file.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ CSV library content write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	# step 3.2
	try:
		log_write("Writing to " + basedir + "result/result.csv")
		csv_file = open(basedir + "result/result.csv","at")
		for line in ["".join([str(x) + ", " for x in [logtime, dict_113[alter], dict_211[alter], measform3(dict_213[alter]), measform3(dict_214[alter]), decision_value] + [dict_115[alter][key] for key in outvar_list_dc] + [dict_122[alter][key] for key in outvar_list_ac] + [measform3(dict_132[alter][key]) for key in outvar_list_ls] + [measform3(dict_143[alter]["vref_avg_proc"] / dict_143[alter]["vref_std_proc"])] + [measform3(dict_143[alter][key]) for key in ["vref_avg_proc", "vref_std_proc", "pow_avg_proc", "pow_std_proc"]] ]) for alter in dict_113_sorted]:
			csv_file.write(line + "\n")
		csv_file.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ CSV result content write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	sub_call("".join(["cd ",basedir,"; rm -rf temp"]))
	sub_call("".join(["cd ",basedir,"; rm -rf output/"]))
