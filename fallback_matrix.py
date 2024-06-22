import math, time
# https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
# matrix is 2D array of matrix[row][col]
# e.g.: [[1, 0], [0, 1]] is a 2 x 2 identity square matrix
def vec_op(vec0, vec1, operation):
	return [operation(x, y) for x, y in zip(vec0, vec1, strict=True)]
def vec_dot(vec0, vec1):
	return sum(vec_op(vec0, vec1, lambda x, y: x * y))
def vec_scale(sca, vec0):
	return [sca * item for item in vec0]
def matrix_almostequal(matrix0, matrix1, epsilon=1e-12):
	return not False in [not False in vec_op(vec0, vec1, lambda x,y: abs(x - y) <= epsilon) for vec0, vec1 in zip(matrix0, matrix1, strict=True)]
def matrix_equal(matrix0, matrix1):
	return matrix_almostequal(matrix0, matrix1, 0)
def matrix_copy(matrix):
	return [[rowcol for rowcol in row] for row in matrix]
def matrix_eye(matrix_or_len):
	if type(matrix_or_len) == type(0):
		rowlen = matrix_or_len
	else:
		rowlen = len(matrix_or_len)
	return [[int(row == col) for col in range(rowlen)] for row in range(rowlen)]
def matrix_trans(matrix):
	return [list(row) for row in zip(*matrix)]
def matrix_linop(matrix0, matrix1, operation):
	return [vec_op(row0, row1, operation) for row0, row1 in zip(matrix0, matrix1, strict=True)]
def matrix_mult(matrix0, matrix1):
	matrix1_trans = matrix_trans(matrix1)
	assert(sum([sum([int(len(row0) == len(row1)) for row0 in matrix0]) for row1 in matrix1_trans]) == len(matrix0) * len(matrix1_trans))
	return [[vec_dot(row0, row1) for row1 in matrix1_trans] for row0 in matrix0]
def matrix_trace(matrix):
	return sum([matrix[z][z] for z in range(min([len(matrix)] + [len(row) for row in matrix]))])
def gauss_jordan_triangle(input_matrix, lower, accumulate_matrix=None):
	assert(not False in [len(input_matrix) == len(row) for row in input_matrix])
	matrix0 = [[rowcol for rowcol in rowval] for rowval in input_matrix]
	rowlen = len(matrix0)
	if accumulate_matrix is not None:
		assert(not False in [rowlen == len(row) for row in accumulate_matrix])
		matrix1 = matrix_copy(accumulate_matrix)
	else:
		matrix1 = matrix_eye(input_matrix)
	try:
		for z0 in range(rowlen):
			index = (z0, rowlen - 1 - z0)[lower]
			if matrix0[index][index] != 0:
				for z1 in range(rowlen):
					if int(z1 != index) * (int(lower) + int(z1 > index)) == 1:
						scale = matrix0[z1][index] / matrix0[index][index]
						for matrix in [matrix0, matrix1]:
							matrix[z1] = vec_op(matrix[z1], vec_scale(scale, matrix[index]), lambda x,y: x - y)
							# equivalent to:
							# (if upper) matrix_mult(I - [[(0, scale)[z0 > index and z1 == index] for z1, col in enumerate(row)] for z0, row in enumerate(matrix)], matrix)
							# i.e. gauss_jordan to get upper triangle matrix is equivalent of prod(lower triangle matrixes) * input_matrix
							# (if lower) matrix_mult(I - [[(0, scale)[z0 < index and z1 == index] for z1, col in enumerate(row)] for z0, row in enumerate(matrix)], matrix)
							# i.e. gauss_jordan to get lower triangle matrix is equivalent of prod(upper triangle matrixes) * input_matrix
	except Exception:
		pass
	return (matrix0, matrix1)
def gauss_jordan_lowertriangle(input_matrix, accumulate_matrix=None):
	return gauss_jordan_triangle(input_matrix, True, accumulate_matrix)
def gauss_jordan_uppertriangle(input_matrix, accumulate_matrix=None):
	return gauss_jordan_triangle(input_matrix, False, accumulate_matrix)
def gauss_jordan_col_manip(matrix, force_assert=False):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	rowlen = len(matrix)
	col_dict = {}
	value_list = []
	ne_0_index = [[z for z, val in enumerate(row) if val != 0] for row in matrix_trans(matrix)]
	for z in range(len([len(x) for x in ne_0_index if len(x) > 0])):
		chosen_index = [len(x) for x in ne_0_index].index(min([len(x) for x in ne_0_index if len(x) > 0]))
		value = min(ne_0_index.pop(chosen_index))
		col_dict[chosen_index] = value
		value_list.append(value)
		for x in ne_0_index:
			if value in x:
				x.remove(value)
	for z in range(rowlen):
		if not z in col_dict:
			assert(not force_assert)
			col_dict[z] = min([x for x in range(rowlen) if not x in value_list])
	return [col_dict[z] for z in range(rowlen)]
def gauss_jordan_col_manip_det(col_order):
	count = 0
	col_copy = [x for x in col_order]
	for z in range(len(col_copy)):
		if col_copy[z] != z:
			index = col_copy[z]
			col_copy[z] = col_copy[index]
			col_copy[index] = index
			count = count + 1
	assert(not False in [x == y for x, y in zip(col_copy, range(len(col_copy)))])
	return (1, -1)[count % 2]
def gauss_jordan_det(input_matrix):
	assert(not False in [len(row) == len(input_matrix) for row in input_matrix])
	rowlen = len(input_matrix)
	col_order = gauss_jordan_col_manip(input_matrix, False)
	right_transform = [[int(col == col_order[row]) for col in range(rowlen)] for row in range(rowlen)]
	right_transform_det = gauss_jordan_col_manip_det(col_order) # work-in-progress, (root_unity(matrix), matrix_trace(matrix), np.linalg.det(matrix))
	# failed, on 8 x 8 matrix, exists two permutation matrices A and B where x^3 == I and trace(x) == 0 for both, and det(A) == -det(B)
	# what is still guaranteed is x^(2z) == I for z integer and det(x) == 1
	# perm_mat = lambda col_order: [[int(col == col_order[row]) for col in range(len(col_order))] for row in range(len(col_order))]
	# perms = lambda z: [perm_mat(col_order) for col_order in itertools.permutations(range(z))]
	# root_unity = lambda perm: [not False in [perm_power[z][z] == 1 for z in range(len(perm_power))] for perm_power in itertools.accumulate([perm for z in range(len(perm) * 2)], lambda x,y: matrix_mult(x,y))].index(True)
	# problem: accumulate only exists in 3.3+, permutations only exists in 2.6+, source code provides sample code for permutations
	# check = lambda z: sorted(sorted(set([(root_unity(perm), np.trace(perm), np.linalg.det(perm)) for perm in perms(z)]), key=lambda x:x[1]), key=lambda x:x[0])
	not_exactly_eigenvalues = [row[z] for z, row in enumerate(gauss_jordan_lowertriangle(matrix_mult(input_matrix, right_transform))[0])]
	if 0 in not_exactly_eigenvalues:
		return 0
	else:
		negatives = len([ersatz_eigen for ersatz_eigen in not_exactly_eigenvalues if ersatz_eigen < 0])
		abs_mult = math.exp(sum([math.log(abs(ersatz_eigen)) for ersatz_eigen in not_exactly_eigenvalues]))
		return right_transform_det * (1, -1)[(negatives) % 2] * abs_mult
def matrix_inverse_gaussjordanable(matrix):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	matrix0, matrix1 = gauss_jordan_uppertriangle(*gauss_jordan_lowertriangle(matrix))
	return [vec_scale(1/matrix0[z][z], row) for z, row in enumerate(matrix1)]
def matrix_inverse(matrix):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	rowlen = len(matrix)
	col_order = gauss_jordan_col_manip(matrix, True)
	right_transform = [[int(col == col_order[row]) for col in range(rowlen)] for row in range(rowlen)]
	# right_inv_tform = [[int(col == col_order[col_order[row]]) for col in range(rowlen)] for row in range(rowlen)]
	# turns out not all permutation matrices have itself as its inverse, only symmetric permutation matrices do
	gauss_jordanable_matrix = matrix_mult(matrix, right_transform)
	if matrix_equal(gauss_jordanable_matrix, matrix_eye(rowlen)):
		return [[int(col == col_order[col_order[row]]) for col in range(rowlen)] for row in range(rowlen)]
	else:
		return matrix_mult(right_transform, matrix_inverse_gaussjordanable(gauss_jordanable_matrix))
def givens_rotation(ndim, i, j, theta):
	return [[((math.sin(theta) * (2 * int(row == i) - 1), math.cos(theta))[col == row], int(col == row))[sum([sum([int(z0 == z1) for z1 in (col, row)]) for z0 in (i, j)]) != 2] for col in range(ndim)] for row in range(ndim)]
def matrix_spectral_ratio(matrix_symmetric):
	return sum([math.pow(matrix_symmetric[z][z], 2) for z in range(len(matrix_symmetric))]) / sum([sum([math.pow(rowcol, 2) for rowcol in row]) for row in matrix_symmetric])
def jacobi_similar(matrix_symmetric, soft_timeout=float("inf")):
	current = [[rowcol for rowcol in row] for row in matrix_symmetric]
	possible = current
	check = lambda a, b: (lambda trace_b: abs(int(trace_b != 0) - matrix_trace(a)/(1, trace_b)[trace_b != 0]) < 0.001)(matrix_trace(b))
	start_time = time.time()
	while check(current, possible) and (matrix_spectral_ratio(current) - matrix_spectral_ratio(possible)) / matrix_spectral_ratio(current) < 0.05 and time.time() - start_time < soft_timeout:
		# loose runtime bound of 180 seconds
		current = possible
		record = (-1, -1)
		cur_max = -float("inf")
		for i, row in enumerate(current):
			for j, rowcol in enumerate(row):
				if i != j and rowcol > cur_max:
					record = (i, j)
					cur_max = rowcol
		x, y = record
		rotation_matrix = givens_rotation(len(current), i, j, 0.5 * math.atan2(current[y][y] - current[x][x], 2 * current[x][y]))
		rotation_matinv = matrix_trans(rotation_matrix)
		possible = matrix_mult(rotation_matinv, matrix_mult(current, rotation_matrix))
	return current
def nullspace(matrix):
	rowlen = max([len(row) for row in matrix])
	square_matrix = [[0 for col in range(rowlen)] for row in range(rowlen)]
	for rowid, row in enumerate(matrix):
		for colid, rowcol in enumerate(row):
			square_matrix[rowid][colid] = rowcol
	col_order = gauss_jordan_col_manip(square_matrix, False)
	right_transform = [[int(col == col_order[row]) for col in range(rowlen)] for row in range(rowlen)]
	right_inv_tform = [[int(col == col_order[col_order[row]]) for col in range(rowlen)] for row in range(rowlen)]
	attempted_diagonal = matrix_mult(gauss_jordan_uppertriangle(gauss_jordan_lowertriangle(matrix_mult(square_matrix, right_transform))), right_inv_tform)
	zero_row_id = [z for z, row in attempted_diagonal if not False in [rowcol == 0 for rowcol in row]]
	null_basis = []
	return null_basis
