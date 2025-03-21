import math
def avg(inputlist):
	return sum(inputlist) / float(len(inputlist))
def std(inputlist, avg_init=None):
	if avg_init is None:
		avg_ = avg(inputlist)
	else:
		avg_ = avg_init
	return math.pow(avg([math.pow(float(x - avg_), 2) for x in inputlist]), 0.5)
def parse_floatint(text):
	try:
		result = (lambda x: (x, int(x))[x == int(x)])(float(text))
	except Exception:
		result = float('nan')
	return result
def select(sql_select, sql_from, sql_where={}):
	# sql_select is an int (column), a list of columns, or "*"
	# sql_from is a table in the form of table[row][col] (2D list)
	# sql_where is a dict of functions, with its key corresponding to rows to check conditions
	assert isinstance(sql_where, type({}))
	if isinstance(sql_select, int):
		select_index = [sql_select]
	elif sql_select == "*":
		select_index = range(min([len(row) for row in sql_from]))
	else:
		select_index = list(sql_select)
	line_pass = [True for row in sql_from]
	for key in sql_where:
		for z, row in enumerate(sql_from):
			if line_pass[z]:
				line_pass[z] = bool(sql_where[key](row[key]))
	return [[row[col] for z, row in enumerate(sql_from) if line_pass[z]] for col in select_index]
	# returns data from passed rows in table in a list per selected column
	# say, for a table like this:
	# 0, 3, 1, 4
	# 1, 1, 5, 9
	# 2, 2, 6, 5
	# 3, 3, 5, 8
	# 4, 9, 7, 9
	# 5, 3, 2, 3
	# select([1,3], table, {2: lambda x: x % 2 == 0}) yields [[2, 3], [5, 3]]
def create_fromcsv(filepath, force_string=False, first_line=None, delimiter=",", comm_sep="."):
	col_names = None
	table = []
	with open(filepath, "rt") as src:
		read_lines = src.readlines()
		if first_line is not None:
			start_index = [str(first_line) in line for line in read_lines].index(True)
			read_lines = read_lines[start_index:]
		for z, line in enumerate(read_lines):
			line_split = line.split(delimiter)
			splitlines = [phrase.lstrip().rstrip() for z, phrase in enumerate(line_split) if z < len(line_split) - 1 or len(phrase.lstrip().rstrip()) > 0]
			if z == 0:
				col_names = splitlines
			else:
				if force_string:
					table.append(splitlines)
				else:
					table.append([parse_floatint(data.replace(comm_sep, ".")) for data in splitlines])
	return col_names, table
def namesearch_index(col_names, criteria=[]):
	result = list(range(len(col_names)))
	for cond in criteria:
		if isinstance(cond, str):
			result = [z for z in result if cond in col_names[z]]
		else:
			result = [z for z in result if cond(col_names[z])]
		#print([col_names[z] for z in result])
	return sorted(result)

#namesearch_index(title, map(lambda y: lambda x: y in x, ["avg", "hi_time"]))[0]