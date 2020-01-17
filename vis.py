#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 14:21:11 2020

@author: akash
"""

import numpy as np
import csv
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import time

def indexer(index, Nx, Ny, iterations):
    base = Nx * Ny
    t = int(index / base)
    i = int((index-(t*Nx*Ny)) / Nx)
    j = index - (t*Nx * Ny) - (Nx*i)
    
    return t,i,j


params = {}

with open('info.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        params[row[0]] = int(row[1])

sizes = [params['RANK '+str(i)+ ' size'] for i in range(params['Ranks'])]
data_dir = './data/'
dump_dir = './temp/'

data = [np.zeros((params['Ny'],params['Nx'])) for i in range(params['Frames_saved'])]

global_i = 0
for data_file_number in range(params['Ranks']):
    with open(data_dir+str(data_file_number)+'.txt', 'r') as file:
        reader = csv.reader(file)
        for row_number, row in enumerate(reader):
            t, i, j = indexer(row_number,params['Nx'],sizes[data_file_number], params['Frames_saved'])
            data[t][global_i + i,j] = row[0]
    global_i += sizes[data_file_number]

# for n, image in enumerate(data):
#     plt.imshow(image)
#     plt.savefig(dump_dir+str(n)+'.png', dpi=300)
#     plt.close('all')

fig = plt.figure()
im=plt.imshow(data[0], vmin=-0.1, vmax=0.1)
def animate(t):
    im.set_array(data[t])
    return [im]       

anim = animation.FuncAnimation(
                               fig, 
                               animate, 
                               frames = len(data),
                               interval = 1000 / 30, # in ms
                               )

anim.save('test_anim.mp4', fps=30, extra_args=['-vcodec', 'libx264'])





