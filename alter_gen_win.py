import pickle, os, sys, time, subprocess
basedir = os.getcwd().replace("\\", "/") + "/"
def log_write(string):
	output_str= str(int(time.time())) + ",;\t,;" + string
	print(output_str[:output_str.index(",")] + ":\t" + string)
	logfile = open(basedir + "run.log", "at")
	logfile.write("\n" + output_str)
	logfile.close()
def sub_call(string):
	log_write("calling \'" + str(string) + "\' on PowerShell 7")
	subprocess.call(str(string), shell=True)
def samelen_dict(inputdict):
	max_length = 0
	for key in inputdict:
		try:
			if type(inputdict[key]) == type(""):
				max_length = max(max_length, 1)
			else:
				max_length = max(max_length, len(inputdict[key]))
		except TypeError:
			max_length = max(max_length, 1)
	newdict = dict()
	for key in inputdict:
		try:
			if type(inputdict[key]) == type("") :
				newdict[key] = [inputdict[key] for index in range(max_length)]
			else:
				newdict[key] = [inputdict[key][index % len(inputdict[key])] for index in range(max_length)]
		except TypeError:
			newdict[key] = [inputdict[key] for index in range(max_length)]
	return (max_length, newdict)
def mainrun(destname):
	try:
		dictpickle = open(basedir + "dict.pickle", "rb")
		dictionary = pickle.load(dictpickle)
		dictpickle.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ pickle load: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise
	sub_call(''.join(["cd \"",basedir.replace("/", "\\"),"\"; rm dict.pickle"]))
	old_times = []
	try:
		row, dict_ = samelen_dict(dictionary)
		dest = open(basedir + destname, "wt")
		for z in range(row):
			dest.write(".alter\n")
			for key in dict_:
				if key == "title":
					trunc_time = dict_[key][z][max(dict_[key][z].index("_") + 1, dict_[key][z].index(".") - 7):dict_[key][z].index(".")]
					while int(trunc_time) in old_times:
						trunc_time = max(old_times) + 1
					old_times.append(int(trunc_time))
					dest.write(".param creation_time="+str(int((0, 0)["selectedparents.pickle" in os.listdir(basedir)]) + int(int(trunc_time) % int(1e7))).zfill(7)+"\n")
				else:
					dest.write(".param "+key+"="+str((round(dict_[key][z]*1000)/1000, int(dict_[key][z]))[dict_[key][z] == int(dict_[key][z])])+"\n")
		dest.close()
	except Exception:
		log_write(str(sys.exc_info()[0]) + " @ write: " + str(sys.exc_info()[1]) + " > " + str(sys.exc_info()[2]))
		raise	
mainrun(sys.argv[1])
