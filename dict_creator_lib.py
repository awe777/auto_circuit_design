import os, random, pickle, time, math, cmath
import numpy as np
import purecma as pcma
import cma
from numpy import linalg as npLA
def context_builder(var_list_dicttuple = {}, rspread = 0.4):
	tuple_list = [('random_spread', rspread)]
	tuple_list = tuple_list + [(context_keys, dict([(key, var_list_dicttuple[z]) for key in var_list_dicttuple])) for z, context_keys in enumerate(["original", "original_unit", "original_min", "original_max"])]
	return dict(tuple_list)
def curdir_file_win(filename="", force_dir=None):
    return (str(force_dir), os.getcwd().replace("\\", "/"))[force_dir is None] + "/" + str(filename)
def log_write(string):
	output_str= str(int(time.time())) + ",;\t,;" + string
	print(output_str[:output_str.index(",")] + ":\t" + string)
	logfile = open(curdir_file_win("run.log"), "at") # "run.log"
	logfile.write("\n" + output_str)
	logfile.close()
# mathematic functions
def avg(inputlist):
	return sum(inputlist) / float(len(inputlist))
def std(inputlist, avg_init=None):
	if avg_init is None:
		avg_ = avg(inputlist)
	else:
		avg_ = avg_init
	return math.pow(avg([float(x - avg_) * float(x - avg_) for x in inputlist]), 0.5)
def normal_2():
	base = [1 - random.random(), random.random()]
	out = cmath.rect(math.sqrt(2 * math.log(1/base[0])), 2 * math.pi * base[1])
	return (out.real, out.imag) 
    # returns 2 samples from normal distribution
	# can be deprecated given random.gauss(0, 1) exists, but can be useful in airgapped systems
def getsingle(fomlist):
	acclist = [max(fomlist[0], 0)]
	for z in range(len(fomlist) - 1):
		acclist.append(acclist[z] + max(fomlist[z + 1], 0))
	if sum(fomlist) != 0:
		acclist = [value / sum(fomlist) for value in acclist]
	else:
		acclist = [(z + 1) / len(fomlist) for z in range(len(fomlist))]
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
	return result
def getcouple(fomlist):
	return getsample(fomlist, 2)
def crossover(dict0, dict1, nu = 13, cr = 0.8, mr = 0.8, mm = 0.01, print_debug = False):
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
	return ((dict0, dict2)[random.random() < cr], (dict1, dict3)[random.random() < cr])
# kernel functions are mostly helpers for Bayesian optimization methods
# value_input for kernel functions is stretched distance sum([alpha[z] * pow(x[z] - y[z], 2.0) for z in range(len(x))]) ** 0.5 of points x and y
"""
def power_exponential_kernel_func(value_input, alpha, gamma):
	assert gamma > 0 and not gamma > 2
	return math.exp(0 - math.pow(abs(value_input / alpha), gamma))
def rational_quadratic_kernel_func(value_input, alpha, gamma):
	assert alpha > 0 and gamma > 0
	return math.pow(1 + 2.0 * math.pow(0.5 * value_input / alpha, 2) / gamma, 0 - gamma)
def factorial(value_input):
	assert type(value_input) == type(0) and not value_input < 0
	result = 1
	for z in range(value_input):
		result = result * (z + 1)
	return result
def gamma_approx(value_input):
	assert value_input != 0
	if type(value_input) == type(0) and value_input > 0:
		# terminal condition
		return factorial(value_input - 1)
	elif int(2.0 * value_input) == 2.0 * value_input and value_input > 0:
		count = int(value_input)
		result = math.sqrt(math.pi)
		for z in range(count):
			result = result * (z + 0.5)
		return result
	elif value_input < 50:
		# domain (-inf, 50), non-integer (integer hits 0 check above)
		return gamma_approx(value_input + 1) / value_input
	else:
		# Stirling approximation for (50, inf)
		stirling_value = value_input
		result = math.sqrt(2 * math.pi * stirling_value) * math.pow(stirling_value / math.e, stirling_value) * sum([math.pow(stirling_value, -float(z)) / float(x) for z, x in enumerate([1, 12, 288, -372.95, -4357.83])])
		return result
def digamma_approx(value_input):
	assert type(value_input) != type(0) or value_input > 0
	dx = 0.05
	lower_bound = math.log(value_input) - 1.0/value_input
	upper_bound = math.log(value_input) - 0.5/value_input
	derivative = (gamma_approx(value_input + dx) - gamma_approx(value_input))/dx
	lower_diff = abs(derivative/gamma_approx(value_input) - lower_bound)
	upper_diff = abs(derivative/gamma_approx(value_input) - upper_bound)
	if lower_diff < upper_diff:
		return lower_bound
	else:
		return upper_bound
def basset_fourier_dirac_comb(order, limit_sum):
	# https://dlmf.nist.gov/10.32#E11
	# for order v > -1/2, x > 0:
	# (x/2)^v * K_v(x) = gamma(v + 0.5) / gamma(0.5) * \int_0^\infty cos(xt)/(1+t^2)^(v+1/2) dt
	assert limit_sum > 0
	fourier_sampling = [gamma_approx(order + 0.5)*math.pow(1 + math.pow(z / float(limit_sum), 2), -0.5 - order)/math.sqrt(math.pi) for z in range(limit_sum + 1)]
	return [(fourier_sampling[z + 1] - fourier_sampling[z]) * 0.5 / limit_sum for z in range(limit_sum)]
# def bessel_k_old(value_input, order):
# 	limit_sum = 1000
# 	absolute_order = abs(order)
# 	e_i_theta = [math.cos(2 * math.pi * (z + 0.5) * value_input/ limit_sum) for z in range(limit_sum)]
# 	if int(2.0 * absolute_order) == 2.0 * absolute_order:
# 		if int(absolute_order) == absolute_order:
# 			# integer order - using recurrence relation starting from K0 and K1 (from Basset's integral)
# 			order_count = absolute_order
# 			order_index = 0
# 			current = [sum([x * y for x, y in zip(e_i_theta, basset_fourier_dirac_comb(z, limit_sum))]) for z in range(2)]
# 		else:
# 			# half-integer order - using recurrence relation starting from K(-0.5)=K0.5 and K0.5 (K0.5 from Basset's integral)
# 			order_count = int(round(absolute_order + 0.5))
# 			order_index = -0.25
# 			current = [sum([x * y for x, y in zip(e_i_theta, basset_fourier_dirac_comb(0.5, limit_sum))]) for z in range(2)]
# 		while order_count > 1:
# 			mult = [1.0, 2.0 * (2 * order_index + 1) / value_input, 2.0 * (2 * order_index + 2) / value_input]
# 			mult.append(1.0 + mult[1]*mult[2])
# 			current = (lambda x,y: [mult[0] * x + mult[1] * y, mult[2] * x + mult[3] * y])(*current)
# 			order_count = order_count - 2
# 			order_index = order_index + 1
# 		return current[order_count]
# 	else:
# 		# 2 * order is not integer, will try hard approximation directly from Basset's integral
# 		return sum([math.pow(value_input / 2.0, -order) * x * y for x, y in zip(e_i_theta, basset_fourier_dirac_comb(order, limit_sum))])
def bessel_k(raw_value_input, order):
	value_input = float(raw_value_input)
	limit_sum = 1000
	e_i_theta = [math.cos(2 * math.pi * (z + 0.5) * value_input / limit_sum) for z in range(limit_sum)]
	absolute_order = abs(float(order))
	order_count = int(absolute_order) + int(absolute_order > int(absolute_order)) # functionally, this is ceil(absolute_order)
	# if order_count is odd, the requested value will be outputted by current[order_count] when order_count has been reduced to 1
	order_of_interest = absolute_order - float(order_count)
	order_index = order_of_interest * 0.5
	current = [sum([x * y for x, y in zip(e_i_theta, basset_fourier_dirac_comb(abs(z + order_of_interest), limit_sum))]) for z in range(2)]
	# with this, basset_fourier_dirac_comb() only needs to be very accurate at order between 0 and 1 (inclusive)
	while order_count > 1:
		mult = [1.0, 2.0 * (2 * order_index + 1) / value_input, 2.0 * (2 * order_index + 2) / value_input]
		mult.append(1.0 + mult[1]*mult[2])
		current = (lambda x,y: [mult[0] * x + mult[1] * y, mult[2] * x + mult[3] * y])(*current)
		order_count = order_count - 2
		order_index = order_index + 1
	return float(current[order_count])
# kernel wrapper and wrapped kernels
def kernel_wrapper(kernel_function, point0, point1, **kwargs): # if there is scaling, "axis_scale_list" must be a named input
	# note: C. E. Rasmussen & C. K. I. Williams, "Gaussian Processes for Machine Learning", MIT Press, p.85, 2006
	# in the above reference, l (length) is a positive parameter constantly re-occuring in its list of covariance functions
	# due to initial work being based of works of Peter I. Frazier, "A Tutorial on Bayesian Optimization", arXiv, 2018:
	# alpha is used in lieu of l (length), gamma is used in lieu of alpha
	assert len(point0) == len(point1)
	if "axis_scale_list" in kwargs:
		as_list = kwargs["axis_scale_list"]
	else:
		as_list = [1 for x in point0]
	distance = math.pow(sum([abs(m * (x - y)) ** 2 for m, x, y in zip(as_list, point0, point1)]), 0.5)
	if distance == 0:
		return 0
	else:
		return kernel_function(distance, **kwargs)
def se_kernel(value_input, **kernel_setup):
	assert not False in [keyword in kernel_setup for keyword in ["alpha"]], [keyword for keyword in ["alpha"] if not keyword in kernel_setup]
	return power_exponential_kernel_func(value_input, math.sqrt(2.0) * float(kernel_setup["alpha"]), 2)
def ou_kernel(value_input, **kernel_setup):
	assert not False in [keyword in kernel_setup for keyword in ["alpha"]], [keyword for keyword in ["alpha"] if not keyword in kernel_setup]
	return power_exponential_kernel_func(value_input, float(kernel_setup["alpha"]), 1)
def gamma_exp_kernel(value_input, **kernel_setup):
	assert not False in [keyword in kernel_setup for keyword in ["gamma", "alpha"]], [keyword for keyword in ["gamma", "alpha"] if not keyword in kernel_setup]
	return power_exponential_kernel_func(value_input, float(kernel_setup["alpha"]), float(kernel_setup["gamma"]))
def rq_kernel(value_input, **kernel_setup):
	assert not False in [keyword in kernel_setup for keyword in ["gamma", "alpha"]], [keyword for keyword in ["gamma", "alpha"] if not keyword in kernel_setup]
	return rational_quadratic_kernel_func(value_input, float(kernel_setup["alpha"]), float(kernel_setup["gamma"]))
def matern_kernel(value_input, **kernel_setup):
	assert not False in [keyword in kernel_setup for keyword in ["order", "alpha"]], [keyword for keyword in ["order", "alpha"] if not keyword in kernel_setup]
	absolute_order = abs(float(kernel_setup["order"]))
	alpha = float(kernel_setup["alpha"])
	assert absolute_order > 0
	normalized_input = float(value_input) * math.sqrt(2.0) / float(alpha)
	if absolute_order > 4: # value of 4 = ceil(3.5) came from C. E. Rasmussen & C. K. I. Williams, "Gaussian Processes for Machine Learning", MIT Press, p.85, 2006
		return power_exponential_kernel_func(value_input, math.sqrt(2.0) * float(kernel_setup["alpha"], 2))
	else:
		return 2.0 * math.pow(normalized_input / 2.0, absolute_order) * bessel_k(normalized_input, absolute_order) / gamma_approx(absolute_order)
# the penultimate function that uses the kernel functions above
def data_torture_normal_approximation(x_list, y_list):
	# for x \in R^d in x_list and its associated y = f(x) \in R in y_list,  
	# generates (multiple?) multivariate, len(x_list) = len(y_list) = n-dimensional normal distribution (mu \in R^n, sigma \in R^(n x n)) where getting y_list is very likely 
	# returns a list of (mu, sigma, chance)
	list_of_tuples = []
	return list_of_tuples
"""
# actual machine-learning functions
def create(length, full_random = False, context=context_builder()):
	assert type(int(length)) == type(0)
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	count = 0
	dictpickle = {'title': []}
	for key in list(original):
		dictpickle[key] = []
	for index in range(int(length)):
		title = "sim_" + str(int(time.time() * 1e6))
		randomized = (index > 0) and (random.random() < 0.2 or full_random)
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
			dictpickle[key].append(value)
		dictpickle['title'].append(title[0:-1]+".sp")
		time.sleep(1e-6)
	if count > 0:
		log_write("".join([str(x) for x in ["DEBUG: ", count, " of ", length, " is fully randomized"]]))
	with open(curdir_file_win("dict.pickle") , "wb") as dest:
		pickle.dump(dictpickle, dest, 0)
'''
def create_cma_es(length, full_random = False, context=context_builder()):
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	var_list = sorted(list(original))
	es = cma.CMAEvolutionStrategy([0.5 * (original_max[key] + original_min[key]) for key in var_list], 0.3 * max([0.5 * (original_max[key] - original_min[key]) for key in var_list]), {'popsize': int(length)})
	es_ask = []
	es_ask = es_ask + list(es.ask())
	while len(es_ask) < int(length):
		es_ask = es_ask + list(es.ask())
	nextbatch = dict([(key, []) for key in var_list + ["title"]])
	for generated in es_ask:
		for z0, key in enumerate(var_list):
			value = float(generated[z0])
			if original_unit[key] is not None:
				value = int(round(value / original_unit[key])) * original_unit[key]
			if original_min[key] is not None:
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				value = (value, original_max[key])[value > original_max[key]]
			nextbatch[key].append(value)
		nextbatch["title"].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
	with open(curdir_file_win("cma_es_param.pickle"), "wb") as cma_obj_source:
		pickle.dump(es, cma_obj_source, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("cma_es_param.pickle"))
	with open(curdir_file_win("dict.pickle"), "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("dict.pickle"))
def cma_es_helper(es, outlist, var_list, tell_limit, sample_count):
	es_ask = []
	es_x = []
	es_y = []
	entry_added = False
	try:
		for output in sorted(outlist, key=lambda x: x[0], reverse=True):
			if output[0] > 1:
				try:
					es_y.append(abs(1.0/output[0]))
					#es_y.append(1.0/abs(math.log(output[0])))
					es_x.append([output[1][key] for key in var_list])
				except Exception:
					es_x, es_y = (lambda esx, esy: (list(esx), list(esy)))(zip(*zip(es_x, es_y)))
			if len(es_y) == tell_limit:
				es.tell(es_x, es_y, False) # tell MUST be equal to es.params.mu
				es_x = []
				es_y = []
				entry_added = True
		if len(es_y) > 0:
			if len(es_y) < tell_limit:
				for output in sorted(outlist, key=lambda x: x[0], reverse=True):
					if output[0] > 1 and len(es_y) < tell_limit:
						try:
							es_y.append(abs(1.0/output[0]))
							#es_y.append(1.0/abs(math.log(output[0])))
							es_x.append([output[1][key] for key in var_list])
						except Exception:
							es_x, es_y = (lambda esx, esy: (list(esx), list(esy)))(zip(*zip(es_x, es_y)))
			if len(es_y) == tell_limit:
				es.tell(es_x, es_y, False)
				es_x = []
				es_y = []
				entry_added = True
	except Exception as err:
		if not entry_added:
			log_write("error during telling sequence with no successful tells: " + str(tell_limit))
			raise
		else:
			log_write("error during telling sequence with at least 1 successful tell")
			log_write("error details: " + str(err))
	es_ask = es_ask + list(es.ask())
	while len(es_ask) < sample_count:
		es_ask = es_ask + list(es.ask())
	log_write("DEBUG: # of entries in es_ask:\t" + str(len(es_ask)))
	if len(es_ask) == 0:
		log_write("ERROR: failed to execute cma_es_helper, output list is empty")
		raise RuntimeError("failed to execute cma_es_helper, output list is empty")
	else:
		return es_ask
def cma_es_helper_2(es, outlist, var_list, tell_limit, sample_count):
	es_ask = []
	es_x = []
	es_y = []
	#entry_added = False
	try:
		for output in sorted(outlist, key=lambda x: x[0], reverse=True):
			if output[0] > 1:
				try:
					#es_y.append(abs(1.0/output[0]))
					es_y.append(math.log(1e58-output[0]))
					#es_y.append(1.0/abs(math.log(output[0])))
					es_x.append([output[1][key] for key in var_list])
				except Exception:
					es_x, es_y = (lambda esx, esy: (list(esx), list(esy)))(zip(*zip(es_x, es_y)))
			if len(es_y) == tell_limit:
				es.tell(es_x, es_y, False) # tell MUST be equal to es.params.mu
				es_x = []
				es_y = []
				es_ask = es_ask + list(es.ask())
		if len(es_y) > 0:
			if len(es_y) < tell_limit:
				for output in sorted(outlist, key=lambda x: x[0], reverse=True):
					if output[0] > 1 and len(es_y) < tell_limit:
						try:
							#es_y.append(abs(1.0/output[0]))
							es_y.append(math.log(1e58-output[0]))
							#es_y.append(1.0/abs(math.log(output[0])))
							es_x.append([output[1][key] for key in var_list])
						except Exception:
							es_x, es_y = (lambda esx, esy: (list(esx), list(esy)))(zip(*zip(es_x, es_y)))
			if len(es_y) == tell_limit:
				es.tell(es_x, es_y, False)
				es_x = []
				es_y = []
				es_ask = es_ask + list(es.ask())
	except Exception as err:
		if len(es_ask) == 0:
			log_write("error during telling sequence with no successful tells: " + str(tell_limit))
			raise
		else:
			log_write("error during telling sequence with at least 1 successful tell")
			log_write("error details: " + str(err))
	while len(es_ask) < sample_count:
		es_ask = es_ask + list(es.ask())
	log_write("DEBUG: # of entries in es_ask:\t" + str(len(es_ask)))
	if len(es_ask) == 0:
		log_write("ERROR: failed to execute cma_es_helper_2, output list is empty")
		raise RuntimeError("failed to execute cma_es_helper_2, output list is empty")
	else:
		return es_ask
'''
def regenerate_cma_es(length, outlist, context=context_builder(), maximize=True, force_length=False, raw_weight_multfunc=None, a_cov=2, c_m=1, force_reset=False):
	# WARNING: VERY DIFFICULT TO PORT TO AIR-GAPPED SYSTEMS
	# thinking of implementing this instead: https://github.com/CMA-ES/pycma/tree/development
	# note to self: pickling the ES object is possible: https://github.com/CMA-ES/pycma/issues/126
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	var_list = sorted(list(original))
	#num = int(length)
	# try:
	# 	if not force_reset:
	# 		with open(curdir_file_win("cma_es_param.pickle"), "rb") as cma_obj_source:
	# 			es = pickle.load(cma_obj_source)
	# 			assert type(es) == pcma.CMAES
	# 	else:
	# 		log_write("force reset enabled")
	# 		log_write("initializing cma_es prior hyperparameters by taking current results as initial batch")
	# 		es = pcma.CMAES([0.5 * (original_max[key] + original_min[key]) for key in var_list], 0.5 * max([0.5 * (original_max[key] - original_min[key]) for key in var_list]))
	# except OSError as err:
	# 	log_write("encountered error while trying to read cma_es_param.pickle file: " + str(err))
	# 	log_write("initializing cma_es prior hyperparameters by taking current results as initial batch")
	# 	#es = pcma.CMAES([original[key] for key in var_list], 0.5, 2 * len(var_list))
	# 	#es = pcma.CMAES([original[key] for key in var_list], 0.5)
	# 	es = pcma.CMAES([0.5 * (original_max[key] + original_min[key]) for key in var_list], 0.4 * max([0.5 * (original_max[key] - original_min[key]) for key in var_list]))

	# try:
	# 	es_ask = cma_es_helper_2(es, outlist, var_list, len(es.params.weights), num)
	# except Exception as err:
	# 	log_write("retrying with es.params.mu")
	# 	es_ask = cma_es_helper_2(es, outlist, var_list, es.params.mu, num)
	
	# outlist is a list of [FoM, dict([(key in var_list + ["title"], value)])]
	# a_cov is usually set to 2, setting less than 2 could be useful in noisy functions
	# c_m is usually set to 1, setting less than 1 could be useful in noisy functions
	#'''
	ndim = len(var_list)
	default_num = 4 + int(3 * math.log(ndim))
	num = min(int(default_num * 2), len(outlist))
	num_out = num
	if force_length:
		num_out = max(num_out, length)
	log_write("DEBUG: CMA-ES sample size is set to " + str(num_out))
	var_list_ordered = sorted(var_list)
	sorted_outlist = sorted(outlist, key=lambda x: x[0], reverse=maximize)
	make_new = force_reset
	try:
		if not force_reset:
			with open(curdir_file_win("cma_es_selfparam.pickle"), "rb") as source:
				mean, step_sigma, cov_mat, p_sigma, p_cov, gen_count = pickle.load(source)
	except OSError as err:
		log_write("encountered error while trying to read cma_es_param.pickle file: " + str(err))
		log_write("initializing cma_es prior hyperparameters by taking current results as initial batch")
		make_new = True
	if raw_weight_multfunc is None:
		rw_func = lambda z: 1
	else:
		try:
			rw_func_list = [float(raw_weight_multfunc[z]) for z in range(num)]
			rw_func = lambda z: rw_func_list[z]
		except Exception:
			rw_func = raw_weight_multfunc
	enoi = math.sqrt(2) * math.gamma((1 + ndim)/2) / math.gamma(ndim/2) 
	# math.gamma() is introduced in 3.2, approximation is sqrt(n) * (1 - (4n)^-1 + (21n^2)^-1)
	default_weight = [math.log((1 + num)/ 2) - math.log(1 + z) for z in range(num)]
	raw_weight = [rw_func(z) * x for z, x in enumerate(default_weight)]
	mark = len([w for w in raw_weight if w >= 0])
	sum_w_pos = sum(raw_weight[:mark])
	sum_w_neg = -sum(raw_weight[mark:])
	mu_eff_pos = sum_w_pos * sum_w_pos / sum([w * w for w in raw_weight[:mark]])
	mu_eff_neg = sum_w_neg * sum_w_neg / sum([w * w for w in raw_weight[mark:]])
	c_1 = a_cov / ((ndim + 1.3) ** 2 + mu_eff_pos)
	c_cov = (4 + mu_eff_pos/ndim) / (ndim + 4 + 2 * mu_eff_pos/ndim)
	c_mu = min(1 - c_1, a_cov * (mu_eff_pos + 1 / mu_eff_pos - 1.75) / ((ndim + 2) ** 2 + a_cov * mu_eff_pos / 2))
	c_sigma = (mu_eff_pos + 2) / (ndim + mu_eff_pos + 5)
	d_sigma = 1 + max(0, math.pow((mu_eff_pos - 1) / (ndim + 1), 0.5) - 1) + c_sigma
	alpha_mu = 1 + c_1/c_mu
	alpha_mueff = 1 + 2 * mu_eff_neg / (mu_eff_pos + 2)
	alpha_pdef = (1 - c_1 - c_mu) / (ndim * c_mu)
	act_weight = [(min(alpha_mu, alpha_mueff, alpha_pdef), 1)[z < mark] * w / (sum_w_neg, sum_w_pos)[z < mark] for z, w in enumerate(raw_weight)]
	
	if make_new:
		mean = [avg([entry[1][key] for entry in sorted_outlist]) for key in var_list_ordered]
		step_sigma = avg([original_max[key] - original_min[key] for key in var_list_ordered]) / 3.0
		cov_mat = [[math.pow((0, avg([entry[1][var_list_ordered[row]] for entry in sorted_outlist]) / step_sigma)[row == col], 2) for col in range(ndim)] for row in range(ndim)]
		p_sigma = [0 for z in range(ndim)]
		p_cov = [0 for z in range(ndim)]
		gen_count = 0
	
	# up to this point, the method doesn't use third-party libraries
	cov_eigen_lam, cov_eigen_vec = npLA.eig(np.array(cov_mat))
	# if cov_mat is a symmetric positive definite matrix, then cov_eigen_vec @ np.diag(cov_eigen_lam) @ cov_eigen_vec.transpose() == cov_mat holds true
	# Examples of symmetric positive definite matrixes are np.diag([1,2,3]) or [[4, 12, -16], [12, 37, -43], [-16, -43, 98]] or [[10, 5, 2], [5, 3, 2], [2, 2, 3]]
	# note that numpy may have slight calculation errors: e.g. with [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
	cov_invsqrt = cov_eigen_vec @ np.diag(np.power(cov_eigen_lam, -0.5)) @ cov_eigen_vec.transpose()
	y_iw = [[(entry[1][key] - mean[z0]) / step_sigma for z0, key in enumerate(var_list_ordered)] for entry in sorted_outlist]
	y_w = np.array([sum([w * y_iw[z1][z0] for z1, w in enumerate(act_weight[:mark])]) for z0, key in enumerate(var_list_ordered)])
	mean = np.array([(1 - c_m) * mean[z0] + c_m * sum([w * sorted_outlist[z1][1][key] for z1, w in enumerate(act_weight[:mark])]) for z0, key in enumerate(var_list_ordered)])
	p_sigma = (1 - c_sigma) * np.array(p_sigma) + np.sqrt(c_sigma * (2 - c_sigma) * mu_eff_pos) * cov_invsqrt @ y_w
	step_sigma = step_sigma * np.exp(c_sigma / d_sigma * (npLA.norm(p_sigma) / enoi - 1))
	
	h_sigma = int(npLA.norm(p_sigma) / np.sqrt(1 - np.power(1 - c_sigma, 2 * (gen_count + 1))) < (1.4 + 2 / (ndim + 1)) * enoi)
	p_cov = (1 - c_cov) * np.array(p_cov) + h_sigma * np.sqrt(c_sigma * (2 - c_sigma) * mu_eff_pos) * y_w
	new_cov_0 = (1 + c_1 * (1 - h_sigma) * c_cov * (2 - c_cov) - c_1 - c_mu * sum(act_weight)) * np.array(cov_mat)
	new_cov_1 = c_1 * (lambda x: x.transpose() @ x)(np.atleast_2d(p_cov))
	act_weight = [(ndim / np.power(npLA.norm(cov_invsqrt @ np.array(y_iw[z])), 2), 1)[z < mark] * w for z, w in enumerate(act_weight)]
	new_cov_2 = c_mu * sum([w * (lambda x: x.transpose() @ x)(np.atleast_2d(y_iw[z])) for z, w in enumerate(act_weight)])
	cov_mat = new_cov_0 + new_cov_1 + new_cov_2
	gen_count = gen_count + 1
	with open(curdir_file_win("cma_es_selfparam.pickle"), "wb") as dest:
		pickle.dump((mean.tolist(), float(step_sigma), cov_mat.tolist(), p_sigma.tolist(), p_cov.tolist(), int(gen_count)), dest, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("cma_es_selfparam.pickle"))
	cov_eigen_lam, cov_eigen_vec = npLA.eig(np.array(cov_mat))
	#'''
	nextbatch = dict([(key, []) for key in var_list + ["title"]])
	for z in range(num_out):
		generated = mean + step_sigma * (cov_eigen_vec @ np.diag(np.sqrt(cov_eigen_lam)) @ np.array([random.gauss(0, 1) for key in var_list]))
	# for generated in es_ask:
		for z0, key in enumerate(var_list):
			value = float(generated[z0])
			if original_unit[key] is not None:
				value = int(round(value / original_unit[key])) * original_unit[key]
			if original_min[key] is not None:
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				value = (value, original_max[key])[value > original_max[key]]
			nextbatch[key].append(value)
		nextbatch["title"].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
	# with open(curdir_file_win("cma_es_param.pickle"), "wb") as cma_obj_source:
	# 	pickle.dump(es, cma_obj_source, 0)
	# 	log_write("DEBUG: successful write to " + curdir_file_win("cma_es_param.pickle"))
	for z in range(int(num_out * 0.4)):
		for z0, key in enumerate(var_list):
			value = random.uniform(original_min[key], original_max[key])
			if original_unit[key] is not None:
				value = int(round(value / original_unit[key])) * original_unit[key]
			if original_min[key] is not None:
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				value = (value, original_max[key])[value > original_max[key]]
			nextbatch[key].append(value)
		nextbatch["title"].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
	log_write("DEBUG: total sample size is set to " + str(len(nextbatch["title"])))
	with open(curdir_file_win("dict.pickle"), "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("dict.pickle"))
'''
def regenerate_cma_es_lib(length, outlist, context=context_builder(), force_reset=False):
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	var_list = sorted(list(original))
	num = int(length)
	try:
		if not force_reset:
			with open(curdir_file_win("cma_es_param.pickle"), "rb") as cma_obj_source:
				es = pickle.load(cma_obj_source)
				assert type(es) == cma.CMAEvolutionStrategy
		else:
			log_write("force reset enabled")
			log_write("initializing cma_es prior hyperparameters by taking current results as initial batch")
			es = cma.CMAEvolutionStrategy([0.5 * (original_max[key] + original_min[key]) for key in var_list], 0.3 * max([0.5 * (original_max[key] - original_min[key]) for key in var_list]))
	except OSError as err:
		log_write("encountered error while trying to read cma_es_param.pickle file: " + str(err))
		log_write("initializing cma_es prior hyperparameters by taking current results as initial batch")
		#es = pcma.CMAES([original[key] for key in var_list], 0.5, 2 * len(var_list))
		#es = pcma.CMAES([original[key] for key in var_list], 0.5)
		es = cma.CMAEvolutionStrategy([0.5 * (original_max[key] + original_min[key]) for key in var_list], 0.4 * max([0.5 * (original_max[key] - original_min[key]) for key in var_list]))
	#log_write("DEBUG: current sigma shape: " + str(np.shape(es.D)))
	log_write("DEBUG: current variances: " + str(es.sm.variances))
	try:
		es_ask = cma_es_helper_2(es, outlist, var_list, es.N_pheno, num)
	except Exception as err:
		log_write("retrying with es.params.mu")
		es_ask = cma_es_helper_2(es, outlist, var_list, es.popsize, num)
	#log_write("DEBUG: current sigma shape: " + str(np.shape(es.D)))
	log_write("DEBUG: current variances: " + str(es.sm.variances))
	nextbatch = dict([(key, []) for key in var_list + ["title"]])
	for generated in es_ask:
		for z0, key in enumerate(var_list):
			value = float(generated[z0])
			if original_unit[key] is not None:
				value = int(round(value / original_unit[key])) * original_unit[key]
			if original_min[key] is not None:
				value = (original_min[key], value)[value > original_min[key]]
			if original_max[key] is not None:
				value = (value, original_max[key])[value > original_max[key]]
			nextbatch[key].append(value)
		nextbatch["title"].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
	with open(curdir_file_win("cma_es_param.pickle"), "wb") as cma_obj_source:
		pickle.dump(es, cma_obj_source, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("cma_es_param.pickle"))
	with open(curdir_file_win("dict.pickle"), "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
		log_write("DEBUG: successful write to " + curdir_file_win("dict.pickle"))
'''
def regenerate_pso(length, outlist, w, phi_p = 2, phi_g = 2, context=context_builder(), force_ignore_vb_dict=False):
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	var_list = sorted(list(original))
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
			with open(curdir_file_win("pso_vb.pickle"), "rb") as source:
				vb_dict = pickle.load(source)
		else:
			vb_dict = {}
	except OSError as err:
		# velocity-best dict pickle file is not found
		log_write("encountered error while trying to read pso_vb.pickle file: " + str(err))
		log_write("fallback to empty vb_dict")
		vb_dict = {}
		# vb_dict is structured as title: [FoM, velocity vector, best position vector]
	if len(outlist) == 0 and len(vb_dict) == 0:
		log_write("warning: outlist and vb_dict length are both 0")
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
	with open(curdir_file_win("pso_vb.pickle"), "wb") as dest:
		pickle.dump(vb_dict, dest, 0)
	nextbatch = dict([(key, []) for key in var_list + ["title"]])
	fom_max = max([vb_dict[particle_title][0] for particle_title in vb_dict])
	new_titles = [pt1 for pt1 in vb_dict if vb_dict[pt1][0] == -float("inf")]
	if len(new_titles) >= num:
		old_titles = []
	else:
		old_titles = (lambda list0: [list0[z] for z in getsample([vb_dict[pt1][0] for pt1 in list0], num - len(new_titles))])(list(vb_dict))
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
	with open(curdir_file_win("dict.pickle"), "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
def regenerate_ga(length, outlist, context=context_builder(), use_prev = False):
	assert type(int(length)) == type(0)
	original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
	var_list = sorted(list(original))
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
	except Exception as err:
		log_write("error: issue with parent extraction - " + str(err))
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
				nextbatch[key].append(value)
			fill = fill + 1
		else:
			None # log_write("debug: add 2")
			for selected in childs:
				nextbatch['title'].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
				for key in [sel_key for sel_key in selected if sel_key != "title"]:
					value = round(int(round(selected[key] / original_unit[key])) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
					if original_min[key] is not None:
						value = (original_min[key], value)[value > original_min[key]]
					if original_max[key] is not None:
						value = (value, original_max[key])[value > original_max[key]]
					nextbatch[key].append(value)
			fill = fill + 2
	with open(curdir_file_win("dict.pickle"), "wb") as dest:
		pickle.dump(nextbatch, dest, 0)
try:
	from bayes_opt import BayesianOptimization # IMPORTANT: REPLACE MODULE WITH UP-TO-DATE SOURCE: https://github.com/bayesian-optimization/BayesianOptimization
	#from bayes_opt import acquisition
	from bayes_opt.logger import JSONLogger
	from bayes_opt.event import Events
	from bayes_opt.util import load_logs, UtilityFunction
	def regenerate_bo(length, outlist, context=context_builder()):#, minimum_fom=0):
		assert type(int(length)) == type(0)
		original, original_unit, original_min, original_max, random_spread = tuple([context[context_keys] for context_keys in ["original", "original_unit", "original_min", "original_max", "random_spread"]])
		var_list = sorted(list(original))
		num = int(length)
		hyperparameter_boundaries = dict([(key, (original_min[key], original_max[key])) for key in var_list])
		
		acquisition_function = UtilityFunction(kind="ucb", kappa=2.5, xi=1)
		try:
			with open(curdir_file_win("bayes_opt_inst.pickle"), "rb") as source:
				bayes_opt_inst = pickle.load(source)
				assert type(bayes_opt_inst) == BayesianOptimization
			#load_logs(bayes_opt_inst, logs=[curdir_file_win("bo_log.json")])
			#log_write("".join([str(x) for x in ["Loaded ", bayes_opt_inst.space," prior points"]]))
			#bayes_opt_inst.subscribe(Events.OPTIMIZATION_STEP, JSONLogger(path=curdir_file_win("bo_log.json")))
		except OSError as err:
			log_write("encountered error while trying to read bayes_opt_inst.pickle file: " + str(err))
			log_write("fallback to a new run")
			bayes_opt_inst = BayesianOptimization(f=None, pbounds=hyperparameter_boundaries, verbose=2)
		sorted_outlist = sorted(outlist, key=lambda x: x[0])
		try:
			for x in sorted_outlist:
				if True: # x[0] >= minimum_fom:
					bayes_opt_inst.register(target=x[0], params=dict([(key, x[1][key]) for key in var_list]))
					# bayes_opt_inst.register(target=x[0], params=x[1]) # could've worked if x[1] doesn't have "title" key in it
		except Exception as err:
			log_write("encountered error while trying to read bo_log.json file: " + str(err))
			log_write("current number of tallied points is " + str(bayes_opt_inst.space))
		fill = 0
		nextbatch = dict([(key, []) for key in var_list + ["title"]])
		while fill < num:
			nextpoint = bayes_opt_inst.suggest(acquisition_function)
			nextbatch["title"].append("sim_" + str(int(time.time() * 1e6))[0:-1] + ".sp")
			for key in var_list:
				value = round(int(round(nextpoint[key] / original_unit[key])) * original_unit[key], max(0, 0 - int(math.log10(original_unit[key]))))
				if original_min[key] is not None:
					value = (original_min[key], value)[value > original_min[key]]
				if original_max[key] is not None:
					value = (value, original_max[key])[value > original_max[key]]
				nextbatch[key].append(value)
			fill = fill + 1
		with open(curdir_file_win("bayes_opt_inst.pickle"), "wb") as dest:
			pickle.dump(bayes_opt_inst, dest, 0)
		with open(curdir_file_win("dict.pickle"), "wb") as dest:
			pickle.dump(nextbatch, dest, 0)
except Exception:
	log_write("WARNING: Bayesian Optimization module is unavailable")
