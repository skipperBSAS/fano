from ROOT import *

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks
# from scipy.signal import find_peaks_cw
import peakutils
import sys
import time
inroot = sys.argv[1]

# ===================================================
# Funcion para "suavizar"  
# ===================================================
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# ===================================================
# Importamos los datos  
# ===================================================

f = TFile(inroot)
skPixTree = f.Get('skPixTree')

#skPixTree.SetBranchStatus("pix",0)
#skPixTree.SetBranchStatus('pix',1)

# ===================================================
# Importamos los datos a un array 
# ===================================================

rows=50 # 50 Lineas 

for i in range(1):
    

    # ===================================================
    # Pasamos datos del Tree a un array  
    # ===================================================

    pixValue = []
    
    for index, event in enumerate(skPixTree):
        #
        if event.y==rows-1: break  # Acotamos el numero de lineas para medir
        if event.x>370 and event.x<7: continue # Quitamos el Overscan
            

        pixValue = np.append(pixValue,event.pix) # Importamos as este array los valores del Branch "pix"
        

    # ===================================================
    # Buscamos la media en la poissoneana del LED  
    # ===================================================

    plt.figure(1)
    
    meanP = np.mean(pixValue)
    range1 = meanP - 2*np.sqrt(meanP*200)   #min(pixValue)
    range2 = meanP + 2*np.sqrt(meanP*200)   #max(pixValue)
    
    #bin=int((range2-range1)/20)
    #n1,bins1,patches1 = plt.hist(pixValue,bins=bin,range=(range1,range2))
    #plt.show()

    bin=int((range2-range1)/50)
    #n1= plt.hist(pixValue,bins=bin,range=(range1,range2))[0]
    n1= plt.hist(pixValue,bins=bin,range=(range1,range1 + 10000))[0]
    plt.show()
    print n1
    indexes = peakutils.indexes(smooth(n1,5), thres=0.02/max(smooth(n1,5)), min_dist= 150)
    
    #peak1= peakutils.interpolate(np.linspace(range1,range2,bin), smooth(n1,5), ind=index1, width=1) #uses a gaussian function and the precedent indexes to enhance our peak finding
    
    print indexes
    #print peaks
    #peaks = [int(i) for i in peaks]
    #print peaks
         # plt.show()
    exit()

    # ===================================================
    # Buscamos la media y sigma de cada pico de carga 
    # ===================================================

    plt.figure(2)

    peak1 = [int(i) for i in peak1]
    range1=peak1[0]-int(np.sqrt(peak1[0]))*4
    range2=peak1[0]+int(np.sqrt(peak1[0]))*4
    bin=(range2-range1)/20
    n,bins,patches = plt.hist(pixValue,bins=bin,range=(range1,range2))
    indexes = peakutils.indexes(n, thres=0.02/max(smooth(n,20)), min_dist= 170*bin/(range2-range1))
    peaks= peakutils.interpolate(np.linspace(range1,range2,bin), n, ind=indexes, width=1) #uses a gaussian function and the precedent indexes to enhance our peak finding
    print indexes
    peaks = [int(i) for i in peaks]
    print peaks
    print peaks
    plt.plot(np.linspace(range1,range2,bin),smooth(n,5))
    # print n[indexes]
    plt.plot(peaks,n[indexes],'bo')
    #plt.show()
    x = np.linspace(1,len(peaks),len(peaks))
    slope = np.polyfit(x,peaks,1)
    print slope
    
    #i+=1
