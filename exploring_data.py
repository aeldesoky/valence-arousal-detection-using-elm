import pickle
import os
import numpy as np
import matplotlib.pyplot as plt

TIME = 60

# write to file
f = open('data/s01.txt', 'rb')
dic = pickle.load(f)
f.close()

y = dic['data']
x = np.arange(0,TIME,TIME/len(y[1][1]))
print(np.shape(y))
print(len(y[1][1]))

plt.subplot(411)
plt.plot(x,y[1][1])
plt.title('First Video')
plt.ylabel('EEG')
plt.subplot(412)
plt.plot(x,y[1][2])
plt.subplot(413)
plt.plot(x,y[1][3])
plt.subplot(414)
plt.plot(x,y[1][4])
plt.subplots_adjust(hspace=0.5)
plt.show()