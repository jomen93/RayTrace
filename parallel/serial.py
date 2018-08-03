import time

start = time.time()


def cube(x):
	suma = 0
	for i in range(1,2):
		suma += i**3
	return suma

results=[]
for i in range(1,20000000):
	results.append(cube(i))

print(results[len(results)-1])

end = time.time()
print (end-start)