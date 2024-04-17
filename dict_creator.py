import os, random, pickle, sys, time, math, cmath
basedir = os.getcwd() + "/"
def log_write(string):
	output_str= str(int(time.time())) + ",;\t,;" + string
	print(output_str[:output_str.index(",")] + ":\t" + string)
	logfile = open(basedir + "run.log", "at")
	logfile.write("\n" + output_str)
	logfile.close()
def avg(inputlist):
	return sum(inputlist) / float(len(inputlist))
def std(inputlist, avg_init=None):
	if avg_init is None:
		avg_ = avg(inputlist)
	else:
		avg_ = avg_init
	return math.pow(avg([pow(float(x - avg_), 2) for x in inputlist]), 0.5)
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
None # log_write("debug: dictionary length assertion success")
random_spread = 0.4
def normal_2():
	base = [1 - random.random(), random.random()]
	out = cmath.rect(math.sqrt(2 * math.log(1/base[0])), 2 * cmath.pi * base[1])
	return (out.real, out.imag) # returns 2 samples from normal distribution
def create(length, full_random = False):
	assert type(int(length)) == type(0)
	count = 0
	dictpickle = {'title': []}
	for key in list(original):
		dictpickle[key] = []
	for index in range(int(length)):
		title = "sim_" + str(int(time.time() * 1e6))
		randomized = full_random or (random.random() < 0.2 and index > 0)
		if randomized and not full_random:
			count = count + 1
		for key in list(original):
			assert(original[key] != None)
			if randomized:
				value = round(int(round(original_max[key] * random.random() / original_unit[key])) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
			else:
				value = round(int(original[key] * (1 - (0, random_spread)[index > 0] + 2 * (0, random_spread)[index > 0] * random.random()) / original_unit[key]) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
			if original_min[key] is not None:
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				value = (value, original_max[key])[value > original_max[key]]
			#if key.endswith("m"):
			#	value = 2 ** int(max(0, math.floor(math.log10(value) / math.log10(2))))
			dictpickle[key].append(value)
			# value = round(value, 2)
			# title = title + str(key) + "-" + (str(value) + "u", str(value).replace(".", "u"))["." in str(value)] + "-"
		dictpickle['title'].append(title[0:-1]+".sp")
		time.sleep(1e-6)
	if count > 0:
		log_write("".join([str(x) for x in ["DEBUG: ", count, " of ", length, " is fully randomized"]]))
	with open(basedir + "dict.pickle", "wb") as dest:
		pickle.dump(dictpickle, dest, 0)

def crossover(dict0, dict1, nu = 13, cr = 0.8, mr = 0.8, mm = 0.01, print_debug = False):
	None # log_write("debug: crossover starts")
	if print_debug:
		log_write("mutation rate: " + str(mr * 100) + "%")
		log_write("mutation multiplier: " + str(mm))
	rand = random.random()
	beta = math.pow((2 * rand, 0.5 / (1 - rand))[rand > 0.5], 1 / (nu + 1))
	dict2 = {}
	dict3 = {}
	for key in [x for x in list(dict0) if x != "title"]:
		normals = normal_2()
		# Liu et al. (ACM 2009) - Equation 11 mutation rate / multiplier (mr/mm)
		dict0_value = dict0[key] * (1, 1 + mm * normals[0])[random.random() < mr]
		dict1_value = dict1[key] * (1, 1 + mm * normals[1])[random.random() < mr]
		# Liu et al. (ACM 2009) - crossover
		dict2[key] = 0.5 * (dict0_value + dict1_value) + (beta / 2) * (dict0_value - dict1_value)
		dict3[key] = 0.5 * (dict0_value + dict1_value) - (beta / 2) * (dict0_value - dict1_value)
	None # log_write("debug: crossover complete")
	return ((dict0, dict2)[random.random() < cr], (dict1, dict3)[random.random() < cr])
def getsingle(fomlist):
	acclist = [max(fomlist[0], 0)]
	for z in range(len(fomlist) - 1):
		acclist.append(acclist[z] + max(fomlist[z + 1], 0))
	if True:
		None # log_write("debug: getsingle - length of acclist is equal to fomlist: " + str(len(acclist) == len(fomlist)))
	if sum(fomlist) != 0:
		acclist = [value / sum(fomlist) for value in acclist]
	else:
		acclist = [(z + 1) / len(fomlist) for z in range(len(fomlist))] # equal opportunity if all FoM = 0
	index = 0
	rand = random.random()
	for z in range(len(acclist)):
		if acclist[z] >= rand:
			index = z
			break
	return index
def getsample(fomlist, number):
	result = []
	for z in range(number):
		result.append(getsingle([(x, 0)[z in result] for z, x in enumerate(fomlist)]))
	#log_write("DEBUG: getsample input - " + str(number))
	#log_write("DEBUG: getsample output - " + str(len(result)))
	return result
def getcouple(fomlist):
	return getsample(fomlist, 2)
def file_exists(path):
	try:
		with open(path, "rb") as source:
			interestlist = pickle.load(source)
			assert type(interestlist) == type([])
			# interestlist is in format of [(FoM, parameter_dict)] -- changed in 20220905
			# open(summary_path).readlines()[1:] contains the variables used for the circuit
		return (True, interestlist)
	except BaseException:
		log_write("error: "+ str(sys.exc_info()[1]) + " @ accessing selectedparents.pickle - " + path)
		return (False, None)
def regenerate_pso(length, outlist, w, phi_p, phi_g, force_ignore_vb_dict=False):
	# w is weight
	# phi_p is cognitive coefficient - self-improvement factor
	# phi_g is social coefficient - global memetic factor, because ~~memes are the DNA of the soul~~ PSO is inherently memetic algorithm
	# u_p, u_g ~ U(0, 1) - random.random()
	# v = w * v + phi_p * u_p * (p - x) + phi_g * u_g * (g - x)
	# x = x + v
	assert type(int(length)) == type(0)
	num = int(length)
	try:
		if not force_ignore_vb_dict:
			with open(basedir + "pso_vb.pickle", "rb") as source:
				vb_dict = pickle.load(source)
		else:
			vb_dict = {}
	except OSError:
		# velocity-best dict pickle file is not found
		vb_dict = {}
		# vb_dict is structured as title: [FoM, velocity vector, best position vector]
	if len(outlist) == 0 and len(vb_dict) == 0:
		log_write("warning: outlist and vb_dict length are both 0")
	# 	return create(int(length))
	sorted_outlist = sorted(outlist, key=lambda x: x[0], reverse=True)
	particle_title_outlist = [particle[1]['title'] for particle in sorted_outlist]
	for particle in sorted_outlist:
		if not particle[1]['title'] in vb_dict:
			# generate initial velocity
			velocity = []
			for key in var_list:
				if original_max[key] is not None:
					minimum = (0, original_min[key])[bool(original_min[key])]
					velocity.append((2 * random.random() - 1) * (original_max[key] - minimum))
				else:
					a, b = (0,0)
					while b == 0:
						a, b = normal_2()
					velocity.append(a / b) # sample from Cauchy distribution
			vb_dict[particle[1]['title']] = [particle[0], velocity, [particle[1][key] for key in var_list]]
		elif particle[0] > vb_dict[particle[1]['title']][0]:
			vb_dict[particle[1]['title']][0] = particle[0]
			vb_dict[particle[1]['title']][2] = [particle[1][key] for key in var_list]
	for pt_negfom in [pt0 for pt0 in vb_dict if vb_dict[pt0][0] < 0]:
		vb_dict.pop(pt_negfom)
	if len(vb_dict) > 0:
		vb_best_key = sorted(list(vb_dict), key=lambda x: vb_dict[x][0], reverse=True)[0]
		f_evol = []
		for z, key in enumerate(var_list):
			values = [vb_dict[particle_title][2][z] for particle_title in vb_dict]
			average = avg(values)
			stdev = std(values, average)
			z_values = [int(stdev != 0) * (value - average)/(1, stdev)[stdev != 0] for value in values]
			best_index = list(vb_dict).index(vb_best_key)
			f_evol.append(2 * math.atan(max(abs(max(z_values) - z_values[best_index]), abs(z_values[best_index] - min(z_values)))) / math.pi)
			# f_evol = 2/pi * atan(max Z-distance from best), 0 <= f_evol < 1, inspired from Zhan's constant (d_g - d_min)/(d_max - d_min)
		w_adj = [math.pow(1 + 1.5 * math.pow(13.5, 0 - f_e), -1) for f_e in f_evol] # 0.4 <= w <= 0.9, directly from Zhan's paper
		# w_adj = [w + (math.pow(1 + 1.5 * math.pow(13.5, 0 - f_evol[z]), -1) - w) * random.random() for z, key in enumerate(var_list)] # 0.4 <= w <= 0.9
		s0 = [max(1, min(0, min(8 - 10 * f_e, 5 * f_e - 2))) for f_e in f_evol]	# exploration
		s1 = [max(1, min(0, min(3 - 5 * f_e, 10 * f_e - 2))) for f_e in f_evol]	# exploiation
		s2 = [max(1, min(0, 1.5 - 5 * f_e)) for f_e in f_evol]					# convergence
		s3 = [max(1, min(0, 5 * f_e - 3.5)) for f_e in f_evol]					# jumping-out
		# s0 to s3 states are directly from Zhan's paper
		delta = 0.05 + 0.05 * random.random() # delta directly comes from Zhan's paper
		phi_p_inc = [delta * (random.random() * s0[z] - random.random() * s3[z]) + delta * (random.random() * s2[z] + random.random() * s1[z]) / 2 for z in range(len(f_evol))]
		phi_g_inc = [delta * (random.random() * s3[z] - random.random() * s0[z]) + delta * (random.random() * s2[z] - random.random() * s1[z]) / 2 for z in range(len(f_evol))]
		phi_p_adj = [(phi_p + phi_p_inc[z]) * min(4 / (phi_p + phi_p_inc[z] + phi_g + phi_g_inc[z]), 1) for z in range(len(f_evol))]
		phi_g_adj = [(phi_g + phi_g_inc[z]) * min(4 / (phi_p + phi_p_inc[z] + phi_g + phi_g_inc[z]), 1) for z in range(len(f_evol))]
		# APSO inspired by Z. H. Zhan et al. "Adaptive Particle Swarm Optimization", IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS 2009
		# new stuff: per-dimension evolutionary factor f_evol (Zhan's f_evol is scalar), f_evol using Z-distance (Zhan's f_evol uses absolute distance, which may not scale well), states being additive instead of definitive
		for particle_title in vb_dict:
			u_p = random.random()
			u_g = random.random()
			for z, key in enumerate(var_list):
				vb_dict[particle_title][1][z] = w_adj[z] * vb_dict[particle_title][1][z]
				if particle_title in particle_title_outlist:
					particle_1_key = sorted_outlist[particle_title_outlist.index(particle_title)][1][key]
				else:
					particle_1_key = vb_dict[particle_title][2][z]
				vb_dict[particle_title][1][z] = vb_dict[particle_title][1][z] + u_p * phi_p_adj[z] * (vb_dict[particle_title][2][z] - particle_1_key)
				vb_dict[particle_title][1][z] = vb_dict[particle_title][1][z] + u_g * phi_g_adj[z] * (vb_dict[vb_best_key][2][z] - particle_1_key)
				if original_max[key] is not None:
					maximum = original_max[key]
					minimum = (0, original_min[key])[bool(original_min[key])]
					vb_dict[particle_title][1][z] = max(min(vb_dict[particle_title][1][z], maximum - minimum), minimum - maximum)
	else:
		s0 = [1 for key in var_list]
		s1 = s0
		s2 = s0
		s3 = s0
	first_pass = True
	while len(vb_dict) < num:
		# add particles based from one of the greatest particle up to 80% of the time, random otherwise
		if random.random() < 0.1 or first_pass:
			log_write(''.join([str(x) for x in ["DEBUG: ", len(vb_dict), " out of ", num]]))
			first_pass = False
		new_title = "sim_" + str(int(random.random() * 1e7 + 7e7)) + ".sp"
		while new_title in vb_dict:
			new_title = "sim_" + str(int(random.random() * 1e7 + 7e7)) + ".sp"
		velocity = []
		for key in var_list:
			if original_max[key] is not None:
				minimum = (0, original_min[key])[bool(original_min[key])]
				velocity.append((2 * random.random() - 1) * (original_max[key] - minimum))
			else:
				a, b = (0,0)
				while b == 0:
					a, b = normal_2()
				velocity.append(a / b)
		if random.random() < 0.04 * min(sum([int(vb_dict[particle_title][0] > 0) for particle_title in vb_dict]) - 1, 20):
			selected_particle_title = list(vb_dict)[getsingle([vb_dict[particle_title][0] for particle_title in vb_dict])]
			vb_dict[new_title] = [vb_dict[selected_particle_title][0], velocity, vb_dict[selected_particle_title][2]]
		else:
			position = []
			for key in var_list:
				if original_max[key] is not None:
					minimum = (0, original_min[key])[bool(original_min[key])]
					value = minimum + random.random() * (maximum - minimum)
				else:
					a, b = (0,0)
					check_pass = b != 0
					while not check_pass:
						a, b = normal_2()
						value = abs(a / b)
						if b != 0:
							check_pass = value > float((bool(original_min[key]), original_min[key])[bool(original_min[key])])
				position.append(value)
			vb_dict[new_title] = [-float("inf"), velocity, position]
	with open(basedir + "pso_vb.pickle", "wb") as dest:
		pickle.dump(vb_dict, dest, 0)
	nextbatch = dict([(key, []) for key in var_list + ["title"]])
	fom_max = max([vb_dict[particle_title][0] for particle_title in vb_dict])
	#log_write("DEBUG: " + str(math.lcm(*[decimal.Decimal((1, vb_dict[pt1][0])[vb_dict[pt1][0] > 0]).as_integer_ratio()[1] for pt1 in vb_dict])))
	#log_write("DEBUG: " + str(max([vb_dict[pt1][0] for pt1 in vb_dict])))
	#log_write("DEBUG: " + str((lambda a: [a for pt0 in vb_dict])(1)))
	#log_write("DEBUG: " + str((lambda x: [int(vb_dict[pt0][0] * (2**63 - 1) / x) for pt0 in vb_dict if vb_dict[pt0][0] > 0])(max(list(map(lambda pt0: vb_dict[pt0][0], vb_dict))))))
	#log_write("DEBUG: " + str((lambda x: list(map(lambda pt0: max(0, round(vb_dict[pt0][0] * (2**31 - 1) / x)), vb_dict)))(max([vb_dict[pt1][0] for pt1 in vb_dict]))))
	#for particle_title in random.sample(list(vb_dict), counts=(lambda x:[max(0, int(round(vb_dict[pt0][0] * (2**31 - 1) / x))) for pt0 in vb_dict])(max([vb_dict[pt1][0] for pt1 in vb_dict])), k=num):
	new_titles = [pt1 for pt1 in vb_dict if vb_dict[pt1][0] == -float("inf")]
	if len(new_titles) >= num:
		old_titles = []
	else:
		#old_titles = (lambda list0: random.sample(list0, counts=(lambda x: [int(vb_dict[pt0][0] * (2**63 - 1) / (x * len(list0))) for pt0 in list0])(max(list(map(lambda pt0: vb_dict[pt0][0], list0)))), k=num - len(new_titles)))([pt1 for pt1 in vb_dict if vb_dict[pt1][0] > 0])
		old_titles = (lambda list0: [list0[z] for z in getsample([vb_dict[pt1][0] for pt1 in list0], num - len(new_titles))])(list(vb_dict))
		#old_titles = (lambda list0: [pt0 for z, pt0 in enumerate(list0) if z in getsample([vb_dict[pt1][0] for pt1 in list0], 60)])(list(vb_dict))
	#log_write("DEBUG: " + str((len(new_titles), len(old_titles), len(vb_dict))))
	for particle_title in new_titles + old_titles:
		vb = vb_dict[particle_title]
		nextbatch["title"].append(particle_title)
		lucky_index = int(len(var_list) * random.random())
		clamping = False
		for z, key in enumerate(var_list):
			value = vb[2][z] + (0, vb[1][z])[vb[0] > -float("inf")]
			value = value + (0, s2[z] * normal_2()[0] * (original_max[key] - original_min[key]) / 3)[(vb[0] >= 0.9 * fom_max and z == lucky_index) or fom_max == -float("inf")] 
			# elitist learning directly from Zhan's paper, but we used threshold of 90% * best FoM instead of "the best"
			# forced spread for all if somehow initial vb_dict is empty AND outlist is empty
			if original_unit[key] is not None:
				value = int(round(value / original_unit[key])) * original_unit[key]
			if original_min[key] is not None:
				# clamping = clamping or value < original_min[key]
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				# clamping = clamping or value > original_max[key]
				value = (value, original_max[key])[value > original_max[key]]
			nextbatch[key].append(value)
		if clamping:
			log_write("warning: clamping detected for " + particle_title)
	with open(basedir + "dict.pickle", "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
def regenerate(length, outlist, use_prev = False):
	assert type(int(length)) == type(0)
	num = int(length)
	sorted_outlist = sorted(outlist, key=lambda x: x[0])
	try:
		nextbatch = {}
		if use_prev:
			nextbatch['title'] = [x[1]['title'] for x in sorted_outlist]
		else:
			nextbatch['title'] = []
		fill = len(nextbatch['title'])
		for paramkey in var_list:
			if use_prev:
				templist = []
				for x in sorted_outlist:
					if paramkey in x[1]:
						templist.append(x[1][paramkey])
					else:
						log_write("warning: " + x[1]["title"] + " missing parameter: " + paramkey)
						log_write("appending random multiplier spread of " + str(random_spread * 100) + "% for parameter " + paramkey)
						templist.append(int(original[paramkey] * (1 - random_spread + 2 * random_spread * random.random()) / original_unit[paramkey]) * original_unit[paramkey])
				nextbatch[paramkey] = templist
			else:
				nextbatch[paramkey] = []
	except BaseException:
		log_write("error: issue with parent extraction - "+str(sys.exc_info()[1]))
		log_write("disabling parent transplant")
		nextbatch = dict([(x, []) for x in list(original) + ["title"]])
	toggle_once = True
	while fill < num:
		indexes = getcouple([max(x[0], 0) for x in sorted_outlist])
		None # log_write("debug: fill - "+str(fill))
		None # log_write("debug: rolled - "+str(indexes)[1:-1])
		childs = crossover(sorted_outlist[indexes[0]][1], sorted_outlist[indexes[1]][1], print_debug = toggle_once)
		if toggle_once:
			toggle_once = False
			None # log_write("debug: toggle_once disabled")
		if num - fill == 1:
			None # log_write("debug: add 1")
			selected = childs[(0, 1)[random.random() >= 0.5]]
			nextbatch['title'].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
			for key in [sel_key for sel_key in selected if sel_key != "title"]:
				value = round(int(round(selected[key] / original_unit[key])) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
				if original_min[key] is not None:
					value = (original_min[key], value)[value > original_min[key]]
				if original_max[key] is not None:
					value = (value, original_max[key])[value > original_max[key]]
				#if key.endswith("m"):
				#	value = 2 ** int(max(0, math.floor(math.log10(value) / math.log10(2))))
				nextbatch[key].append(value)
			fill = fill + 1
		else:
			None # log_write("debug: add 2")
			for selected in childs:
				nextbatch['title'].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
				for key in [sel_key for sel_key in selected if sel_key != "title"]:
					value = round(int(round(selected[key] / original_unit[key])) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
					None # log_write("debug: current key - " + str(key))
					if original_min[key] is not None:
						value = (original_min[key], value)[value > original_min[key]]
					if original_max[key] is not None:
						value = (value, original_max[key])[value > original_max[key]]
					#if key.endswith("m"):
					#	value = 2 ** int(max(0, math.floor(math.log10(value) / math.log10(2))))
					nextbatch[key].append(value)
			fill = fill + 2
	with open(basedir + "dict.pickle", "wb") as dest:
		pickle.dump(nextbatch, dest, 0)

check = file_exists(basedir + "selectedparents.pickle")
boolval = len(sys.argv) > 2 and sys.argv[2] == "True"
floatval = float((0.4 + 0.5 * random.random(), sys.argv[-1])[len(sys.argv) > 2])
#if True and check[0] and len(([], check[1])[check[0]]) >= 4: # Liu et al. (ACM 2009)
if check[0]:
	if len(check[1]) == 0:
		log_write("warning: outlist list is empty")
	log_write("dictionary creator - PSO")
	regenerate_pso(sys.argv[1], check[1], floatval, 2, 2) 
	# initial idea of w = (1 - cycle/max_cycle) * (0.9 - 0.4) + 0.4, phi_p = phi_g = 2 comes from PSO-LDIW as written by X. L. Wang et al., "A High-Efficiency Design Method of TSV Array for Thermal Management of 3-D Integrated System", TCAD 2023
	# end up unused
else:
	if boolval:
		log_write("dictionary creator - random - does not use base value - trigger: " + str(sys.argv[2:]))
	else:
		log_write("dictionary creator - random - uses base value")
	create(int(sys.argv[1]), boolval)

