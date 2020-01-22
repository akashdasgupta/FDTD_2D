#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:21:11 2020
@author: akash
"""
import matplotlib
matplotlib.use('Qt4Agg')


import numpy as np
import csv
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import scipy.spatial as sp


def indexer(index, Nx, Ny):
    base = Nx * Ny
    t = int(index / base)
    i = int((index-(t*Nx*Ny)) / Nx)
    j = index - (t*Nx * Ny) - (Nx*i)
    
    return t,i,j


data_dir = './data/'
dump_dir = './temp/'

params = {}

with open(data_dir+'info.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        params[row[0]] = int(row[1])


ep_holder = []
with open (data_dir+"epmap.txt", 'r') as file:
    reader=csv.reader(file)
    for row in reader:
        ep_holder.append(float(row[0]))
        
image = np.zeros((params['Ny'],params['Nx']))

for k, element in enumerate(ep_holder):
    t,i,j = indexer(k,params['Ny'],params['Nx'])
    image[i,j] = element

points= []

for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        if image[i-1,j] != image[i,j]:
            points.append((j,i))
points = np.array(points)
hull = sp.ConvexHull(points)

# inner = [(params['PML_size'],params['PML_size']),
#          (params['Nx']-params['PML_size'],params['PML_size']), 
#          (params['Nx']-params['PML_size'],params['Ny']-params['PML_size']), 
#          (params['PML_size'],params['Ny']-params['PML_size'])]

# outer = [(0,0), 
#          (params['Nx'],0), 
#          (params['Nx'],params['Ny']), 
#          (0,params['Ny'])]

sizes = [params['RANK '+str(i)+ ' size'] for i in range(params['Ranks'])]

data = [np.zeros((params['Ny'],params['Nx'])) for i in range(params['Frames_saved'])]

global_i = 0
for data_file_number in range(params['Ranks']):
    with open(data_dir+str(data_file_number)+'.txt', 'r') as file:
        reader = csv.reader(file)
        for row_number, row in enumerate(reader):
            t, i, j = indexer(row_number,params['Nx'],sizes[data_file_number])
            data[t][global_i + i,j] = row[0]
    global_i += sizes[data_file_number]

fig = plt.figure()
im=plt.imshow(data[0], vmin=-0.07, vmax=0.07, cmap='jet', interpolation='none')
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')


def animate(t):
    im.set_array(data[t])
    return [im]       

anim = animation.FuncAnimation( fig, 
                                animate, 
                                frames = len(data),
                                interval = 1000 / 30, # in ms
                                )

anim.save('test_anim.mp4', fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)




