#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 23:29:55 2019

@author: akash
"""
import numpy as np
import csv
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import time

#import matplotlib
#matplotlib.use("Qt5Agg")

start_time = time.time()
def index(i,j,t, Nx, Ny):
    return((Nx*Ny*(t)) + (i*Nx) + j)

Nx = 100
Ny = 100
num_ranks = 4
iterations = 590

filenames = ['0TOP.txt']
for i in np.arange(1,num_ranks,1):
    filenames.append(str(i)+'.txt') 
filenames.append('0BOTTOM.txt')


ny_sizes = []
base_ny = int(Ny / num_ranks)
ny_remainder = Ny % num_ranks

ny_sizes.append(int(base_ny/2))

for i in np.arange(1,num_ranks,1):
    if i <= ny_remainder:
        ny_sizes.append(base_ny+1)
    else:
        ny_sizes.append(base_ny)

ny_sizes.append(base_ny - int(base_ny/2))

        
        
        
        
filepath= r'/home/akash/4th_year_computing/FDTD_2D/data'
data = {}

for filename in filenames:
    holder = []
    path = filepath + '/' + filename
    
    with open(path,'r') as file:
        reader=csv.reader(file)
        for row in reader:
            holder.append(row)
    data[filename] = holder


image = np.zeros((Nx,Ny))

ims = []
fig = plt.figure()


for t in np.arange(0,iterations,1):
    global_i = 0
    for filename, local_ny in zip(filenames, ny_sizes):
        for i in range(local_ny):
            for j in range(Nx):
                E = data[filename][index(i,j,t,Nx,local_ny)][3]
                image[global_i, j] = E
            global_i += 1
    ims.append([plt.imshow(image**2, cmap='inferno', vmin=0, vmax=0.01)]) 
    plt.show()               

ani = animation.ArtistAnimation(fig,ims, interval=33)
ani.save('dynamic_images.mp4')
plt.show()

print(time.time() - start_time)
