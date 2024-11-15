import sys, pickle, math, dict_creator_lib, csv_data_sifter
# import random
def measform3(x):
	try:
		return str(x * pow(10, -int(math.floor(math.log10(abs(x))))))+"e"+str(int(math.floor(math.log10(abs(x)))))
	except Exception:
		return str(x)
def measform3_short(x):
	try:
		return (lambda l, m, r: l[:5 + (0, 1)[int(l.startswith("-"))]] + m + r)(*measform3(x).partition("e")).replace("\'", "")
	except Exception:
		return str(x).replace("\'", "")
basedir = dict_creator_lib.curdir_file_win()

var_list = ["op1_stage1_pmos_w","op1_stage1_pmos_l","op1_stage1_pmos_m","op1_stage1_nmos_w","op1_stage1_nmos_l","op1_stage1_nmos_m","op1_stage2_nmos_w","op1_stage2_nmos_l","op1_stage2_nmos_m","op1_vref_nmos_w","op1_vref_nmos_l","op1_vref_nmos_m","op1_vref_pmos_w","op1_vref_pmos_l","op1_vref_pmos_m","op1_stage2_pmos_w","op1_stage2_pmos_l","op1_stage2_pmos_m","psrr_stabilizer_w","psrr_stabilizer_l","psrr_stabilizer_m","res_rr0_w","res_rr0_l","res_rr0_m","res_rr2_w","res_rr2_l","res_rr2_m"]
original = dict(map(lambda k, v: (k, v), var_list, [10,18,4,12,1,8,8,1,4,14,8,4,5,8,4,8,4,8,2.75,4,8,1,14,1,1,8.55,2]))
original_unit = dict([(key, (0.025, 1)[key.endswith("m")]) for key in var_list])
original_min = dict([(key, ((1, 5)[key.endswith("ratio")], 0.18 + (0, 0.04)[key.endswith("w")])[key.endswith("w") or key.endswith("l")]) for key in var_list])
original_max = dict([(key, (20, 25)[key.endswith("w")]  / (1, 4)[key.startswith("op1_stage2")]) for key in var_list])
assert(len(var_list) == len(original))
for x in [original_unit, original_min, original_max]:
	assert(len(original) == len(x))
random_spread = 0.4

#current_context = dict([(context_keys, locals()[context_keys]) for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
current_context = {}
current_context["original"] = original
current_context["original_unit"] = original_unit
current_context["original_min"] = original_min
current_context["original_max"] = original_max
current_context["random_spread"] = random_spread

def file_exists(path):
	try:
		with open(path, "rb") as source:
			interestlist = pickle.load(source)
			assert type(interestlist) == type([])
			# interestlist is in format of [(FoM, parameter_dict)] -- changed in 20220905
			# open(summary_path).readlines()[1:] contains the variables used for the circuit
		return (True, interestlist)
	except Exception:
		dict_creator_lib.log_write("error: "+ str(sys.exc_info()[1]) + " @ accessing selectedparents.pickle - " + path)
		return (False, None)

check = file_exists(basedir + "selectedparents.pickle")
boolval = len(sys.argv) > 2 and sys.argv[2] == "True"
#floatval = float((0.4 + 0.5 * random.random(), sys.argv[-1])[len(sys.argv) > 2])
floatval = 0
#if True and check[0] and len(([], check[1])[check[0]]) >= 4: # Liu et al. (ACM 2009)
try:
	if check[0]:
		if len(check[1]) == 0:
			dict_creator_lib.log_write("warning: outlist list is empty")
		length = int(sys.argv[1])
		cma_es_len = int(2 * (4 + int(3 * math.log(len(var_list)))))
		dict_creator_lib.log_write("dictionary creator - CMA-ES using recommended length (" + str(cma_es_len) + ")" + ("", ", accompanied by GA - " + str(max(length - cma_es_len, 0)))[1 - int(boolval)])
		best_param = None
		best_param_temp = {}
		#dict_creator_lib.log_write("dictionary creator - BO")
		try:
			list_order = None
			col_names, table = csv_data_sifter.create_fromcsv(basedir + "result/result.csv")
			fom1_index = csv_data_sifter.namesearch_index(col_names, ["FoM_1"])[0]
			fom2_index = csv_data_sifter.namesearch_index(col_names, ["FoM_2"])[0]
			time_index = csv_data_sifter.namesearch_index(col_names, ["time"])[0]
			idex_index = csv_data_sifter.namesearch_index(col_names, ["identifier"])[0]
			#for fom_list in csv_data_sifter.select([fom2_index, fom1_index], table): # this includes the initial circuit
			for fom_list in csv_data_sifter.select([fom2_index, fom1_index], table, {idex_index: lambda id: id != 10000000}): # this omits the initial circuit
				if list_order is None:
					list_order = sorted(range(len(fom_list)), key=lambda z: fom_list[z], reverse=True)
				else:
					list_order = sorted(list_order, key=lambda z: fom_list[z], reverse=True)
			best_time, best_idex = tuple([value[list_order[0]] for value in csv_data_sifter.select([time_index, idex_index], table, {idex_index: lambda id: id != 10000000})])
			fom1_best, fom2_best = tuple([value[list_order[0]] for value in csv_data_sifter.select([fom1_index, fom2_index], table, {idex_index: lambda id: id != 10000000})])
			dict_creator_lib.log_write("DEBUG: best current (time, id):\t" + str((best_time, best_idex)))
			dict_creator_lib.log_write("DEBUG: best FoM(s):\t\t" + str((measform3_short(fom1_best), measform3_short(fom2_best))).replace("\'", ""))
			col_names, table = csv_data_sifter.create_fromcsv(basedir + "result/library.csv")
			best_param_temp = dict([(var_list[z], keyvalues[0]) for z, keyvalues in enumerate(csv_data_sifter.select([csv_data_sifter.namesearch_index(col_names, [var])[0] for var in var_list], table, {time_index:lambda x: x == best_time, idex_index: lambda x: x == best_idex}))])
			best_param_temp["title"] = "sim_"+str(best_idex)+".sp" # remnants of old program config in the form of sim_#.sp
			# dict_creator_lib.log_write("DEBUG: best param constructed")
			# dict_creator_lib.log_write("DEBUG: all var_list data taken into account: " + str(not False in [var in best_param_temp for var in var_list]))
			best_param = best_param_temp
		except Exception as err:
			dict_creator_lib.log_write("error: " + str(err) + " @ CSV sifting - " + str(best_param_temp))
			best_param = None
		#dict_creator_lib.log_write("dictionary creator - GA")
		#dict_creator_lib.regenerate_cma_es(length, check[1], current_context, best_param, True, boolval)
		cma_es_dict = dict_creator_lib.regenerate_cma_es(cma_es_len, check[1], current_context, best_param, True, boolval, True, 20)
		if boolval or cma_es_len >= length: # True forces CMA-ES only
			dict_creator_lib.dict_pickle_dump(cma_es_dict)
		else:
			dict_creator_lib.dict_pickle_dump(cma_es_dict, dict_creator_lib.regenerate_ga(length - cma_es_len, check[1], current_context))
		#dict_creator_lib.regenerate_ga(sys.argv[1], check[1], current_context)
	else:
		if boolval:
			dict_creator_lib.log_write("dictionary creator - random - does not use base value - trigger: " + str(sys.argv[2:]))
		else:
			dict_creator_lib.log_write("dictionary creator - random - uses base value")
		dict_creator_lib.dict_pickle_dump(dict_creator_lib.create(int(sys.argv[1]), boolval, current_context))
except Exception as err:
	dict_creator_lib.log_write("Error during dict creation: " + str(err))
	raise
