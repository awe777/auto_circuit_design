import os, subprocess, time, pickle
basedir = os.getcwd().replace("\\", "/") + "/"
with open(basedir + "hspice.env", "rb") as hspice_pickle:
	hspice_working_dir, hspice_exe_path = pickle.load(hspice_pickle)
cur_env = basedir[:-1].rpartition("/")[2]
generation_limit = 500
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
	return subprocess.Popen(str(string), shell=True, stdin=popen_in, stdout=subprocess.PIPE).communicate()[0].decode(encoding="utf-8")

def stage0(trigger=False):
	log_write("stage 0, "+str(trigger)+" - " + str(time.ctime()))
	sub_call("".join(["cd \'",basedir,"\'; python dict_creator.py 60"]) + ("", " True")[trigger])
def stage1():
	log_write("stage 1 - " + str(time.ctime()))
	sub_call("".join(["cd \'",basedir,"\'; python alter_gen_win.py BGR.alter"]))
def stage2():
	log_write("stage 2 - " + str(time.ctime()))
	sub_call("cd " + hspice_working_dir + "; mkdir " + cur_env)
	sub_call("echo \'\' > " + hspice_working_dir + cur_env + "/BGR.lis")
	sub_call("".join(["cd \'",basedir,"\'; " + hspice_exe_path + " -i BGRN_new.sp -o " + hspice_working_dir + cur_env + "/BGR.lis -mp 16 -mt 16 > NUL"]).replace("/", "\\"))
	sub_call("cd " + hspice_working_dir + cur_env + "/; rm *.dp*")
	with open(hspice_working_dir + cur_env + "/BGR.lis", "rt") as lis_file:
		for line in lis_file:
			if "**error**" in line:
				log_write("error in LIS file: " + line.partition("**error**")[2].lstrip().rstrip())
	sub_call("echo \'\' > " + hspice_working_dir + cur_env + "/BGR.lis")
	sub_call("".join(["cd \'",basedir,"\'; rm BGR.alter"]))
def stage3():
	log_write("stage 3 - " + str(time.ctime()))
	sub_call("".join(["cd \'",basedir,"\'; mkdir output"]))
	sub_call("".join(["cd \'",basedir,"\'; mkdir result"]))
	sub_call("".join(["cd \'",basedir,"output/\'; copy " + hspice_working_dir + cur_env + "/*.csv ./"]))
	#sub_call("".join(["cd \'",basedir,"output/\'; copy " + hspice_working_dir + cur_env + "/*.lis ./"]))
	sub_call("".join(["cd \'",basedir,"\'; python extract_analysis.py"]))
generation_count = 0
try:
	open(basedir + "dont_delete.txt", "rt").close()
except OSError:
	sub_call("".join(["cd \'",basedir,"\'; rm run.log"]))
	log_write("deleting old pickle files, logs, and folders")
	sub_call("".join(["cd \'",basedir,"\'; rmdir output -Force -Recurse"]))
	sub_call("".join(["cd \'",basedir,"\'; rmdir result -Force -Recurse"]))
	sub_call("".join(["cd \'",basedir,"\'; rm BGR.alter"]))
	sub_call("".join(["cd \'",basedir,"\'; rm *.pickle"]))
	sub_call("".join(["cd \'",basedir,"\'; echo \'\' > dont_delete.txt"]))
log_write("I am "+sub_popen("whoami").lstrip().rstrip()+", starting at " + str(time.ctime()))
while generation_count < generation_limit:
	log_write("generation count: " + str(generation_count) + " - " + str(time.ctime()))
	try:
		open(basedir + "BGR.alter", "rt").close()
		log_write("resuming interrupted session")
	except OSError:
		stage0(False)
		stage1()
	stage2()
	stage3()
	generation_count = generation_count + 1