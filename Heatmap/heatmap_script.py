import os
import numpy as np
import matplotlib.pyplot as plt

def test_heatmap(n_of_external_threads, n_of_internal_threads, heat_matrix):
	for external in range(0, n_of_external_threads, 2):
		for internal in range(0, n_of_internal_threads, 2):
			cmd = './heat ' + str(external) + ' ' + str(internal) 
			heat_value = os.popen(cmd).read()
			print(heat_value)
			heat_matrix[external][internal] = heat_value


def draw(heat_matrix):
	plt.imshow(heat_matrix, cmap='hot', interpolation='nearest')
	plt.show()

if __name__ == "__main__":
	n_of_external_threads = 8
	n_of_internal_threads = 32
	heat_matrix = np.zeros((n_of_external_threads, n_of_internal_threads))
	
	test_heatmap(n_of_external_threads, n_of_internal_threads, heat_matrix)
	draw(heat_matrix)