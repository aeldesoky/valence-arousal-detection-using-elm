import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from scipy.fftpack import fft, fftfreq
from loading_data import change_label_values_to_calss

CHANNELS_NUM = 40
VIDEOS_NUM   = 40
CHANNEL_DATA_POINTS = 8064
VIDEO_LENGTH_SECONDS = 60

for file in sorted(os.listdir('data')):

   # write to file
   data_file = open('data/' + os.fsdecode(file), 'rb')
   data_dic = pickle.load(data_file)
   data_file.close()

   all_channels_data  = data_dic['data']
   all_labels_values  = data_dic['labels']
   all_labels_classes = change_label_values_to_calss(all_labels_values)
   valence_labels     = all_labels_values[:, 0].reshape(VIDEOS_NUM, 1)
   arousal_labels     = all_labels_values[:, 1].reshape(VIDEOS_NUM, 1)
   dominance_labels   = all_labels_values[:, 2].reshape(VIDEOS_NUM, 1)
   liking_labels      = all_labels_values[:, 3].reshape(VIDEOS_NUM, 1)

   time = np.arange(0, VIDEO_LENGTH_SECONDS, VIDEO_LENGTH_SECONDS/CHANNEL_DATA_POINTS)
   time_step = 1/128
   samples_num = CHANNEL_DATA_POINTS
   frequency = np.linspace(0.0, 1.0/(2.0*time_step), samples_num//2)

   print(np.shape(frequency))

   for video_num in range(0, VIDEOS_NUM):
      for channel_num in range(0, CHANNELS_NUM):
         channel_data = all_channels_data[video_num][channel_num].reshape(CHANNEL_DATA_POINTS, 1)
         channel_data_fft = fft(channel_data)
         plt.plot(frequency, 2.0/samples_num * np.abs(channel_data_fft[0:samples_num//2]))
         plt.plot()
      
      plt.title('Video number ' + str(video_num))
      plt.ylabel('EEG Amplitude')
      plt.xlabel('Time in Seconds')   
      plt.show()