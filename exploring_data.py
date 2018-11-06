import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tqdm import tqdm

CHANNELS_NUM = 40
VIDEOS_NUM   = 40
VIDEO_LENGTH_SECONDS = 60

for file in sorted(os.listdir('data')):

   # write to file
   data_file = open('data/' + os.fsdecode(file), 'rb')
   data_dic = pickle.load(data_file)
   data_file.close()

   all_channels_data = data_dic['data']
   all_labels        = data_dic['labels']
   valence_labels    = all_labels[:, 0].reshape(VIDEOS_NUM, 1)
   arousal_labels    = all_labels[:, 1].reshape(VIDEOS_NUM, 1)
   dominance_labels  = all_labels[:, 2].reshape(VIDEOS_NUM, 1)
   liking_labels     = all_labels[:, 3].reshape(VIDEOS_NUM, 1)

   time = np.arange(0,VIDEO_LENGTH_SECONDS, VIDEO_LENGTH_SECONDS/len(all_channels_data[1][1]))

   for video_num in range(1,VIDEOS_NUM):
      for channel_num in range(1,CHANNELS_NUM):
         channel_data = all_channels_data[video_num][channel_num]
         plt.plot(time,all_channels_data[video_num][channel_num])
      
      plt.title('Video number ' + str(video_num) + ' (' + str(all_labels[video_num]) + ')')
      plt.ylabel('EEG')
      plt.xlabel('Time in Seconds')   
      plt.show()