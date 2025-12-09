import numpy as np

x=np.array([i for i in range(10)])
print(x)
print(x.reshape(2,5)[::-1,:])