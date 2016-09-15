from math import log
import numpy as np

def ReadData(filename):
	dataX = []
	dataY = []
	errorY = []

	file = open(filename, 'r')
	with file as f:
		for line in f:
			x, y, ey = line.split('\t')
			#x, y = line.split()
			dataX.append(float(x))
			dataY.append(float(y))
			errorY.append(float(ey))

	file.close()

	return dataX, dataY, errorY

def combineDataSets(data, abs_data_sets, rel_data_sets):
	mergedDataX = []
	mergedDataY = []
	mergedDataeY = []
	
	for i in range(abs_data_sets):
		mergedDataX += data['absX{0}'.format(i)]
		mergedDataY += data['absy{0}'.format(i)]
		mergedDataeY += data['absey{0}'.format(i)]

	for i in range(rel_data_sets):
		mergedDataX += data['relX{0}'.format(i)]
		mergedDataY += data['rely{0}'.format(i)]
		mergedDataeY += data['reley{0}'.format(i)]		

	return mergedDataX, mergedDataY, mergedDataeY

def createMatrixX(dataX, data, rel_data_sets):
	#Create a Matrix of zero
	x0 = 193.5
	matrix_X = [ [0 for x in range(4+rel_data_sets)] for y in range(len(dataX)) ]
	#This is format for matrix_X
	#(x-x0)^2	1	x 	x0(2x-x0), x^2	rel_1	rel_2	rel_3

	for i in range(len(dataX)):
		if (dataX[i] < x0):
			matrix_X[i][0] = (log(dataX[i]) - log(x0))**2
			matrix_X[i][3] = log(x0)*(2*log(dataX[i]) - log(x0))
		else:
			matrix_X[i][3] = log(dataX[i])**2
		matrix_X[i][1] = 1
		matrix_X[i][2] = log(dataX[i])
		for rel in range(rel_data_sets):
			if dataX[i] in data['relX{0}'.format(rel)]:
				matrix_X[i][rel + 4] = 1

	return matrix_X

def createMatrixY(dataY):
	matrix_Y = [0 for y in range(len(dataY))]

	for i in range(len(dataY)):
		matrix_Y[i] = log(dataY[i])
	return matrix_Y

#Define the chisquare function
def ChiSquared(dataX, dataY, x0, theta):
	#we need to define x_l and x_h
	chi_squared = 0
	for i in range(len(dataX)):
		if dataX[i] < x0:
			x_l = np.matrix([(log(dataX[i]) - log(x0))**2, 1, log(dataX[i]), log(x0)*(2*log(dataX[i]) - log(x0))])
			chi_squared += (x_l*theta[0:4,:] - log(dataY[i]))**2
		else:
			x_h = np.matrix([0, 1, log(dataX[i]), log(dataX[i])**2])
			chi_squared += (x_h*theta[0:4,:] - log(dataY[i]))**2
	return chi_squared

def derChiSquaredWRTx0(dataX, dataY, x0, theta):
	derChiSquared = 0
	for i in range(len(dataX)):
		if (dataX[i] < x0):
			derChiSquared += theta.item(0,0)*(-2*(log(dataX[i]) - log(x0))) + theta.item(3,0)*(2*log(dataX[i])-log(x0))
	return derChiSquared	

def evaluateTheta(X, Y):
	XT = X.transpose()
	return (np.linalg.inv(XT*X)*XT*(Y.transpose()))

def GradientDescent(dataX, dataY, init_x0, data, rel_data_sets):
	#theta has to evaluated for each new x0
	precision = 0.05
	gamma = 0.01
	x_new = init_x0 + 3*precision
	x_old = init_x0
	while abs(x_new - x_old) > precision:
		getX = createMatrixX(dataX, data, rel_data_sets)
		getY = createMatrixY(dataY)
		#print (getY)
		X = np.matrix(getX)
		Y = np.matrix(getY)
		theta = evaluateTheta(X, Y)
		x_old = x_new
		x_new = x_old - gamma*derChiSquaredWRTx0(dataX, dataY, x_old, theta)
	return x_new


def main():
	###Read the absolute and relative efficiency data sets
	data = {}

	abs_data_sets = int(input('How many absolute efficiency data sets: '))

	for i in range(abs_data_sets):
		filename = input('Input the name of the file: ')
		dataX, dataY, errorY = ReadData(filename)
		#print(dataX)
		data['absX{0}'.format(i)] = dataX
		data['absy{0}'.format(i)] = dataY
		data['absey{0}'.format(i)] = errorY

	rel_data_sets = int(input('How many relative efficiency data sets: '))

	for i in range(rel_data_sets):
		filename = input('Input the name of the file: ')
		dataX, dataY, errorY = ReadData(filename)
		#print(dataX)
		data['relX{0}'.format(i)] = dataX
		data['rely{0}'.format(i)] = dataY
		data['reley{0}'.format(i)] = errorY

	mergedDataX, mergedDataY, mergedDataeY = combineDataSets(data, abs_data_sets, rel_data_sets)
	getX = createMatrixX(mergedDataX, data, rel_data_sets)
	getY = createMatrixY(mergedDataY)
	#print (getY)
	X = np.matrix(getX)
	Y = np.matrix(getY)
	theta = evaluateTheta(X, Y)
	#theta = np.linalg.inv(XT*X)*XT*(Y.transpose())
	x0 = 150
	x_0 = GradientDescent(mergedDataX, mergedDataY, x0, data, rel_data_sets)
	#chi = ChiSquared(mergedDataX, mergedDataY, x0, theta)
	print(x_0)
	#print(np.linalg.inv(XT*X)*XT*(Y.transpose()))


main()