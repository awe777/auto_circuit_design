import sys, os, pickle, dict_creator_lib, random
basedir = os.getcwd() + "/"
# last cycle was 1705524931406114
# [10.5, 20.0, 4.0, 15.1, 0.8, 9.0, 16.2, 4.9, 4.0, 19.3, 11.8, 3.0, 8.9, 15.1, 3.0, 11.6, 7.3, 5.0, 6.9, 8.7, 4.0, 1.4, 12.1, 1.0, 2.5, 18.8, 2.0]
var_list = ["vref_op_pch3_w", "vref_op_pch3_l", "vref_op_pch3_m", "vref_op_nch_w", "vref_op_nch_l", "vref_op_nch_m", "vref_op_stg2_w", "vref_op_stg2_m", "vref_op_stg2_l", "vref_op_join_l", "vref_op_outn_l", "vref_op_out3_l", "stage0_pch_w", "stage0_pch_l", "stage0_pch_m", "op1_stack0_w", "op1_stack0_l", "op1_stack0_m", "op1_stack1_w", "op1_stack1_l", "op1_stack1_m", "op1_join_w", "op1_join_l", "op1_join_m", "op2_input_w", "op2_input_l", "op2_input_m", "op2_base_w", "op2_base_l", "op2_base_m", "res_rr0_w", "res_rr0_l", "res_rr0_m", "res_rr1_w", "res_rr1_l", "res_rr1_m", "res_rr3_w", "res_rr3_l", "res_rr3_m", "res_rr5_w", "res_rr5_l", "res_rr5_m", "psrr_stabilizer_w", "psrr_stabilizer_l", "psrr_stabilizer_m"]
original = dict(map(lambda k, v: (k, v), var_list, [4, 8, 2, 1, 18.4, 2, 2.5, 4, 18.75, 15.1, 20, 15, 50, 20, 1, 25, 20, 8, 12.5, 5, 8, 0.6, 20, 2, 1, 8.5, 8, 1.5, 3, 4, 1, 20, 1, 1, 20, 1, 1, 20, 1, 1, 20, 1, 8.65, 13.4, 8]))
original_unit = dict([(key, (0.025, 1)[key.endswith("m")]) for key in var_list])
original_min = dict([(key, ((1, 5)[key.endswith("ratio")], 0.35 + (0, -0.13)[key.endswith("w")])[key.endswith("w") or key.endswith("l")]) for key in var_list])
original_max = dict([(key, (10, 500)[key.endswith("w")]) for key in var_list])

assert(len(var_list) == len(original))
for x in [original_unit, original_min, original_max]:
	assert(len(original) == len(x))
random_spread = 0.4

current_context = dict([(context_keys, locals()[context_keys]) for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
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
floatval = float((0.4 + 0.5 * random.random(), sys.argv[-1])[len(sys.argv) > 2])
#if True and check[0] and len(([], check[1])[check[0]]) >= 4: # Liu et al. (ACM 2009)
if check[0]:
	if len(check[1]) == 0:
		dict_creator_lib.log_write("warning: outlist list is empty")
	dict_creator_lib.log_write("dictionary creator - PSO")
	dict_creator_lib.regenerate_pso(sys.argv[1], check[1], floatval, 2, 2, current_context) 
	# initial idea of w = (1 - cycle/max_cycle) * (0.9 - 0.4) + 0.4, phi_p = phi_g = 2 comes from PSO-LDIW as written by X. L. Wang et al., "A High-Efficiency Design Method of TSV Array for Thermal Management of 3-D Integrated System", TCAD 2023
	# end up unused
else:
	if boolval:
		dict_creator_lib.log_write("dictionary creator - random - does not use base value - trigger: " + str(sys.argv[2:]))
	else:
		dict_creator_lib.log_write("dictionary creator - random - uses base value")
	dict_creator_lib.create(int(sys.argv[1]), boolval, current_context)

