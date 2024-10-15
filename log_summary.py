# import csv_data_sifter, os, math
import os, math
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
# print("\n")
for folder_root in sorted(os.scandir(os.getcwd().replace("\\","/")), key=lambda item: item.name):
	if folder_root.is_dir() and folder_root.name != "__pycache__":
		# print(folder_root.name + "\n")
		print(folder_root.name)
		line = ""
		# for folder_env in os.scandir(folder_root.path.replace("\\","/")):
		# 	if folder_env.is_dir():
		# 		print(folder_env.name)
		# 		title, table = csv_data_sifter.create_fromcsv(folder_env.path.replace("\\","/") + "/result/result.csv")
		# 		print(''.join(["(" + str(z) + ") " + str(x) + ", " for z, x in enumerate(title)])[:-2])
		# 		list_order = None
		# 		for fom_list in csv_data_sifter.select([4, 3], table, {1: lambda id: id != 10000000}):
		# 			if list_order is None:
		# 				list_order = sorted(range(len(fom_list)), key=lambda z: fom_list[z], reverse=True)
		# 			else:
		# 				list_order = sorted(list_order, key=lambda z: fom_list[z], reverse=True)
		# 		print("".join(["(" + str(z) + ") " + str((x[list_order[0]], measform3_short(x[list_order[0]]))[int(z == 3 or z == 4)]) + ", " for z, x in enumerate(csv_data_sifter.select("*", table, {1: lambda id: id != 10000000}))])[:-2] + "\n")
		with open(folder_root.path.replace("\\","/") + "/run.log", "rt") as logfile:
			z = 0
			for lines in logfile.readlines():
				if "generation count" in lines:
					z = z + 1
					line = lines.rstrip().lstrip()
		print(''.join(["(", str(z), ") ", line, "\n"]))
		