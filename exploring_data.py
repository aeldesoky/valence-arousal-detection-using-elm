import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

TIME = 60

# write to file
f = open('data/s01.txt', 'rb')
dic = pickle.load(f)
f.close()

y = dic['data']
x = np.arange(0,TIME,TIME/len(y[1][1]))

for i in range(1,40):
   plt.plot(x,y[1][i])
   plt.title('First Video')
   plt.ylabel('EEG')

plt.show()