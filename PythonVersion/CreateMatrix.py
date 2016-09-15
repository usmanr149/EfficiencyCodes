def createMatrix(dataX, rel_data_sets):
	matrix_X = []

	for i in range(len(dataX)):
		for j in range(4+rel_data_sets):
			matrix_X[i][j] = 0

#(x-x0)	1	x 	x0(2x-x0), x^2	rel_1	rel_2	rel_3

for i in range(len(dataX)):
	if (x < x_0):
		matrix_X[i][0] = (x[i] - x_0)
		matrix_X[i][3] = x0(2*x[i] - x0)
	else:
		matrix_X[i][3] = x[i]^2
	matrix_X[i][1] = 1
	matrix_X[i][2] = x[i]
	if x[i] in rel_1:
		matrix_X[i][4] = 1
	elif x[i] in rel_2:
		matrix_X[i][5] = 1
	elif x[i] in rel_3:
		matrix_X[i][6] = 1