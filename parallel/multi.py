import multiprocessing as mp
import time

start = time.time()

def cube(x):
	suma = 0
	for i in range(1,2):
		suma += i**3
	return suma


pool = mp.Pool(processes=2)
results = pool.map(cube, range(1,20000000))
print(results[len(results)-1])

end = time.time()
print (end-start)