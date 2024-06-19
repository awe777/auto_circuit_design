import math
# https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
# matrix is 2D array of matrix[row][col]
# e.g.: [[1, 0], [0, 1]] is a 2 x 2 identity square matrix
def vec_op(vec0, vec1, operation):
	return [operation(x, y) for x, y in zip(vec0, vec1, strict=True)]
def vec_dot(vec0, vec1):
	return sum(vec_op(vec0, vec1, lambda x, y: x * y))
def vec_scale(sca, vec0):
	return [sca * item for item in vec0]
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
		matrix1 = [[rowcol for rowcol in row] for row in accumulate_matrix]
	else:
		matrix1 = [[int(row == col) for col in range(rowlen)] for row in range(rowlen)]
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
def gauss_jordan_det(input_matrix):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	rowlen = len(matrix)
	col_order = gauss_jordan_col_manip(matrix, True)
	right_transform = [[int(col == col_order[row]) for col in range(rowlen)] for row in range(rowlen)]
	not_exactly_eigenvalues = [row[z] for z, row in enumerate(gauss_jordan_triangle(input_matrix, True)[0])]
	if 0 in not_exactly_eigenvalues:
		return 0
	else:
		negatives = len([ersatz_eigen for ersatz_eigen in not_exactly_eigenvalues if ersatz_eigen < 0])
		abs_mult = math.exp(sum([math.log(abs(ersatz_eigen)) for ersatz_eigen in not_exactly_eigenvalues]))
		return (1, -1)[negatives % 2] * abs_mult
def matrix_inverse_gaussjordanable(matrix):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	matrix0, matrix1 = gauss_jordan_uppertriangle(*gauss_jordan_lowertriangle(matrix))
	return [vec_scale(1/matrix0[z][z], row) for z, row in enumerate(matrix1)]
def gauss_jordan_col_manip(matrix, force_assert=False):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	rowlen = len(matrix)
	col_dict = {}
	ne_0_index = [[z for z, val in enumerate(row) if val != 0] for row in matrix_trans(matrix)]
	for z in sorted([z for z, row in enumerate(ne_0_index) if len(row) > 0], key=lambda x: len(ne_0_index[z])):
		col_dict[z] = min([x for x in ne_0_index[z] if not x in list(zip(*col_dict.items()))[1]])
	for z in range(rowlen):
		if not z in col_dict:
			assert(not force_assert)
			col_dict[z] = min([x for x in range(rowlen) if not x in list(zip(*col_dict.items()))[1]])
	return [col_dict[z] for z in range(rowlen)]
def matrix_inverse(matrix):
	assert(not False in [len(row) == len(matrix) for row in matrix])
	rowlen = len(matrix)
	col_order = gauss_jordan_col_manip(matrix, True)
	right_transform = [[int(col == col_order[row]) for col in range(rowlen)] for row in range(rowlen)]
	right_inv_tform = [[int(col == col_order[col_order[row]]) for col in range(rowlen)] for row in range(rowlen)]
	return matrix_mult(right_inv_tform, matrix_inverse_gaussjordanable(matrix_mult(matrix, right_transform)))
def givens_rotation(ndim, i, j, theta):
	return [[((math.sin(theta) * (2 * int(row == i) - 1), math.cos(theta))[col == row], int(col == row))[sum([sum([int(z0 == z1) for z1 in (col, row)]) for z0 in (i, j)]) != 2] for col in range(ndim)] for row in range(ndim)]
def matrix_spectral_ratio(matrix_symmetric):
	return sum([math.pow(matrix_symmetric[z][z], 2) for z in range(len(matrix_symmetric))]) / sum([sum([math.pow(rowcol, 2) for rowcol in row]) for row in matrix_symmetric])
def jacobi_similar(matrix_symmetric):
	current = [[rowcol for rowcol in row] for row in matrix_symmetric]
	possible = current
	check = lambda a,b: abs(1 - matrix_trace(a)/matrix_trace(b)) < 0.001 if matrix_trace(b) != 0 else abs(matrix_trace(a) - matrix_trace(b)) < 0.001
	while check(current, possible) and 1 - matrix_spectral_ratio(possible)/matrix_spectral_ratio(current)) < 0.05:
		current = possible
		record = (-1, -1)
		cur_max = -float("inf")
		for i, row in enumerate(current):
			for j, rowcol in enumerate(row):
				if i != j and rowcol > cur_max:
					record = (i, j)
					cur_max = rowcol
		x, y = record
		if current[x][x] == current[y][y]:
			angle = math.pi/4
		else:
			angle = 0.5 * math.atan(2 * current[x][y] / (current[y][y] - current[x][x]))
		rotation_matrix = givens_rotation(len(current), i, j, angle)
		rotation_matinv = matrix_trans(rotation_matrix)
		possible = matrix_mult(rotation_matinv, matrix_mult(current, rotation_matrix))
	return current	
