import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import mne
import math
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
   print(all_labels_values)
   break
   ch_names = ['Fp1', 'AF3', 'F3', 'F7', 'FC5', 'FC1', 'C3', 'T7', 'CP5', 'CP1', 'P3', 'P7', 'PO3', 'O1', 'Oz', 'Pz', 
               'Fp2', 'AF4', 'Fz', 'F4', 'F8', 'FC6', 'FC2', 'Cz', 'C4', 'T8', 'CP6', 'CP2', 'P4', 'P8', 'PO4', 'O2']
   ch_types = ['eeg'] * 32
   montage  = mne.channels.read_montage('standard_1020', ch_names=ch_names)
   info =  mne.create_info(ch_names=ch_names, sfreq=128.0, ch_types=ch_types)


   time = np.arange(0, VIDEO_LENGTH_SECONDS, VIDEO_LENGTH_SECONDS/CHANNEL_DATA_POINTS)
   fs = 128.0
   time_step = 1/128.0
   samples_num = CHANNEL_DATA_POINTS
   frequency = np.linspace(0.0, 1.0/(2.0*time_step), samples_num//2)

   for video_num in range(0, VIDEOS_NUM):
      video_data = all_channels_data[video_num]
      alpha_data = mne.filter.filter_data(video_data, fs, 8, 13, method='iir')
      beta_data = mne.filter.filter_data(video_data, fs, 13, 30, method='iir')
      gamma_data = mne.filter.filter_data(video_data, fs, 30, 40, method='iir')
      theta_data = mne.filter.filter_data(video_data, fs, 4, 8, method='iir')

      raw = mne.io.RawArray(video_data[0:32][:], info)
      raw.set_montage(montage, set_dig=True)
      #raw.plot(n_channels=32, scalings='auto', title='Auto-scaled Data from arrays', show=True, block=True)
      # raw.plot_psd(show=True, proj=True, spatial_colors=True, dB=True)

      # raw = mne.io.RawArray(alpha_data[0:32][:], info)
      # raw.set_montage(montage, set_dig=True)
      # raw.plot_psd(show=True, proj=True, spatial_colors=True, dB=True)

      # raw = mne.io.RawArray(beta_data[0:32][:], info)
      # raw.set_montage(montage, set_dig=True)
      # raw.plot_psd(show=True, proj=True, spatial_colors=True, dB=True)

      # raw = mne.io.RawArray(gamma_data[0:32][:], info)
      # raw.set_montage(montage, set_dig=True)
      # raw.plot_psd(show=True, proj=True, spatial_colors=True, dB=True)

      # raw = mne.io.RawArray(theta_data[0:32][:], info)
      # raw.set_montage(montage, set_dig=True)
      # raw.plot_psd(show=True, proj=True, spatial_colors=True, dB=True)

      for channel_num in range(0, CHANNELS_NUM):
         channel_data = video_data[channel_num].reshape(CHANNEL_DATA_POINTS, 1)

         channel_data_fft = fft(channel_data)
         plt.plot(frequency, 2.0/samples_num * np.abs(channel_data_fft[0:samples_num//2]))
         plt.plot()
      
         plt.title('Video number ' + str(video_num))
         plt.ylabel('EEG Amplitude')
         plt.xlabel('Time in Seconds')   
         plt.show()