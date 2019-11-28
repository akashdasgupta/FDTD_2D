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

import matplotlib
matplotlib.use("Qt5Agg")

Nx = 100
Ny = 100

filepath= r'/home/akash/4th_year_computing/FDTD_2D/data'
files = [filepath + '/' + str(i+1) for i in range(1000)]

ims = []
fig = plt.figure()

count=0

for filename in files:
    holderez=[]
    holderhx = []
    holderhy = []

    with open(filename + '.txt','r') as file:
        reader=csv.reader(file)
        for row in reader:
            holderez.append(float(row[0]))
            holderhx.append(float(row[2]))
            holderhy.append(float(row[3]))
    
    
    imageEz = np.zeros((Ny, Nx))
    imageHx = np.zeros((Ny, Nx))
    imageHy = np.zeros((Ny, Nx))

    for i in range(Ny):
        for j in range(Nx):
            index = (i * Nx) + j
            imageEz[i,j] = holderez[index]
            imageHx[i,j] = holderhx[index]
            imageHy[i,j] = holderhy[index]
#    if count ==10:
#        ims.append([plt.imshow(imageEz, cmap='inferno')])
#        count = 0
#    else:
#        count += 1
    ims.append([plt.imshow((imageEz)**2, cmap='inferno', vmin=0, vmax=1)])

ani = animation.ArtistAnimation(fig,ims, interval=33)
ani.save('dynamic_images.mp4')
plt.show()

