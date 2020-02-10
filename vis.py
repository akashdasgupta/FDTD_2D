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
import scipy.spatial as sp

object_type = 'lens' # lens, tunnell or none (double slit wont need edge drawn)
scale =0.5           # min/max values of the plot (contrast tuning parameter)

def indexer(index, Nx, Ny):
    base = Nx * Ny
    t = int(index / base)
    i = int((index-(t*Nx*Ny)) / Nx)
    j = index - (t*Nx * Ny) - (Nx*i)
    
    return t,i,j

# Loads parameters for the run:
data_dir = './data/'
params = {}
with open(data_dir+'info.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        try:
            params[row[0]] = int(row[1])
        except:
            params[row[0]] = float(row[1])

# Loads ep map to see what the object's map is:
ep_holder = []
with open (data_dir+"epmap.txt", 'r') as file:
    reader=csv.reader(file)
    for row in reader:
        ep_holder.append(float(row[0]))
ep_image = np.zeros((params['Ny'],params['Nx']))

for k, element in enumerate(ep_holder):
    t,i,j = indexer(k,params['Ny'],params['Nx'])
    ep_image[i,j] = element

# creates set of points where permitivity changes:
points= []
for i in range(ep_image.shape[0]):
    for j in range(ep_image.shape[1]):
        if ep_image[i-1,j] != ep_image[i,j]:
            points.append((i,j))
points = np.array(points)

# creates hull if lens was chosen:
if object_type == 'lens':
    hull = sp.ConvexHull(points)

sizes = [params['RANK '+str(i)+ ' size'] for i in range(params['Ranks'])]
data = [np.zeros((params['Ny'],params['Nx'])) for i in range(params['Frames_saved'])]

# loads data into created list:
global_i = 0
for data_file_number in range(params['Ranks']):
    with open(data_dir+str(data_file_number)+'.txt', 'r') as file:
        reader = csv.reader(file)
        for row_number, row in enumerate(reader):
            t, i, j = indexer(row_number,params['Nx'],sizes[data_file_number])
            data[t][global_i + i,j] = row[0]
    global_i += sizes[data_file_number]

fig = plt.figure()
sizer = ((params['Nx']/2) * params['spacestep']) / 1e-6

im=plt.imshow(data[0], vmin=-scale, vmax=scale, cmap='jet', interpolation='none', extent=[-sizer,sizer,-sizer,sizer])


cbar=plt.colorbar()
cbar.set_label(r"E Field / $\sqrt{Js^{-1}}$")
plt.xlabel(r"x ($\mu$m)")
plt.ylabel(r"y ($\mu$m)")


plt.tight_layout()
stepscale = params['spacestep'] / 1e-6

if object_type == 'lens':
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0]*stepscale - sizer, -points[simplex, 1]*stepscale + sizer, 'k-')

if object_type == 'tunnelling':
    for i in range(points.shape[0]):
        if points[i,1]*stepscale - sizer > -3.5:
            plt.scatter(-points[i,0]*stepscale + sizer, points[i,1]*stepscale - sizer, color='k', s=2)
    lims = params['Ny'] * params['spacestep'] / 2
    plt.xlim((-lims,lims))
    plt.ylim((-lims,lims))

def animate(t):
    im.set_array(data[t])
    return [im]       

anim = animation.FuncAnimation( fig, 
                                animate, 
                                frames = len(data),
                                interval = 1000 / 30, # in ms
                                )

anim.save('test_anim.mp4', fps=15, extra_args=['-vcodec', 'libx264'], dpi=300)




