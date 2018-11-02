import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
TIME = 60
#data_file  = open(os.path.expanduser('s01.dat'), 'rb')
#s01 = pickle.load(data_file, encoding='latin1')
#pickle.dump(s01, open('test.txt', 'wb'))

# write to file
f = open('test.txt', 'rb')
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
plt.show()


#print((dic['labels']))
#print((dic['data']))