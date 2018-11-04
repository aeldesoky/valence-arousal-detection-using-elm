import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

directory_string = '~/Desktop/DEAP/data_preprocessed_python'
directory = os.path.expanduser(directory_string)
os.makedirs('data', exist_ok=True)

print("Importing data...")
for file in tqdm(sorted(os.listdir(directory))):
   filename = os.fsdecode(file)
   data_file_path = os.path.join(directory, filename)
   if filename.endswith(".dat"):
      data_file = open(data_file_path, 'rb')
      pickle_file = pickle.load(data_file, encoding='latin1')
      text_file = open(os.path.join('data/', os.path.splitext(filename)[0]) + ".txt", 'wb')
      pickle.dump(pickle_file, text_file)
      data_file.close()
      text_file.close()