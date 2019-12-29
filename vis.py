
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

def indexer(index, Nx, Ny, iterations):
    base = Nx * Ny
    t = int(index / base)
    i = int((index-(t*Nx*Ny)) / Nx)
    j = index - (t*Nx * Ny) - (Nx*i)
    
    return t,i,j

Nx = 200
Ny = 200
num_ranks = 5
iterations = 300
pml=int(Nx/20)

filenames = ['0TOP.txt']
for i in np.arange(1,num_ranks,1):
    filenames.append(str(i)+'.txt') 
filenames.append('0BOTTOM.txt')

mpi_lines = []

ny_sizes = []
base_ny = int((Ny-2*pml) / (num_ranks-1))
ny_remainder = (Ny-2*pml) % (num_ranks-1)

ny_sizes.append(pml)


for i in np.arange(1,num_ranks,1):
    if i-1 < ny_remainder:
        ny_sizes.append(base_ny+1)
    else:
        ny_sizes.append(base_ny)

ny_sizes.append(pml)

line = 0
for i in ny_sizes:
    line += i
    mpi_lines.append(line)

mpi_linex = range(Nx)
     
#filepath= r'/newhome/ad16020/FDTD_2D/data'
filepath= r'/home/akash/4th_year_computing/FDTD_2D/data'
# data = {}

# for filename in filenames:
#     holder = []
#     path = filepath + '/' + filename
    
#     with open(path,'r') as file:
#         reader=csv.reader(file)
#         for row in reader:
#             holder.append(row)
#     data[filename] = holder


images = [np.zeros((Nx,Ny)) for i in range (int(iterations))]



global_i = 0
for filename, local_ny in zip(filenames, ny_sizes):
    path = filepath + '/' + filename
    with open(path,'r') as file:
        reader=csv.reader(file)
        counter = 0
        for row in reader:
            try:
                t,i,j = indexer(counter,Nx,local_ny, iterations)
                images[t][global_i+i,j] = float(row[0])
                counter += 1
                #print(j)
                #print(i)
            except:
                pass
        global_i += local_ny
        
ims = []
fig = plt.figure()

for i, image in enumerate(images):
    ims.append([plt.imshow(image, vmin=-0.001, vmax=0.001, cmap='jet')])
    #plt.imshow(image, vmin=-0.025, vmax=0.025, cmap='jet')
    #for size in mpi_lines:
        #plt.plot(mpi_linex, np.ones(len(mpi_linex)) * size, color='w')
    
    #plt.savefig('/newhome/ad16020/FDTD_2D/temp/'+str(i)+'.png', dpi=100)    
    #plt.savefig('/home/akash/4th_year_computing/FDTD_2D/temp/'+str(i)+'.png', dpi=100)

    #plt.close('all')
    
    
ani = animation.ArtistAnimation(fig,ims, interval=33)
ani.save('dynamic_images.mp4')
plt.show()
#plt.close('all')
print(time.time() - start_time)
