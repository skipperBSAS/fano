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
# Importamos los datos
# ===================================================

rows=400*50 # 50 Lineas

for i in range(1):


    # ===================================================
    # Pasamos datos del Tree a un array
    # ===================================================


    pixValue = []
    #print 1
    for index, event in enumerate(skPixTree):
        #
        if index==rows: # Acotamos el numero de lineas para medir
            #print i, index
            break

        pixValue = np.append(pixValue,event.pix) # Importamos as este array los valores del Branch "pix"
        #
        # # print i
        #  # print f.skPixTree[i]
        #  # print event.pix, index
        #  # index = 1
        #  # print i, (i%40000)

        #
        #     # index = 40000*(i+1)
        #     # i+=1
        #     # print index, i





    # ===================================================
    # Buscamos la media en la poissoneana del LED
    # ===================================================

    plt.figure(1)

    range1=100000
    range2=550000
    bin=(range2-range1)/100
    n1,bins1,patches1 = plt.hist(pixValue,bins=bin,range=(range1,range2))
    index1 = peakutils.indexes(smooth(n1,5), thres=0.02/max(smooth(n1,5)), min_dist= bin)
    peak1= peakutils.interpolate(np.linspace(range1,range2,bin), smooth(n1,5), ind=index1, width=1) #uses a gaussian function and the precedent indexes to enhance our peak finding
    #print indexes
    #print peaks
    #peaks = [int(i) for i in peaks]
    #print peaks
         # plt.show()


    # ===================================================
    # Buscamos la media y sigma de cada pico de carga
    # ===================================================



    peak1 = [double(i) for i in peak1]
    range1=peak1[0]-int(np.sqrt(peak1[0]))*10
    range2=peak1[0]+int(np.sqrt(peak1[0]))*10
    BIN=int((range2-range1)/30) # enlarging the dividend, gives us thicker bins
    N=101
    bin_vector=np.linspace(BIN-N/2,BIN+N/2, num=N)
    m=0
    delta =[0]*N
    dist_peaks =((),)

    #chi=np.linspace(BIN-N/2,BIN+N/2, num=N)
    for h in bin_vector:
        #j=j+1
        print m
        #print "bines = "+str(bin_vector[m])
        plt.figure(m+2)
        n,bins,patches = plt.hist(pixValue,bins=bin_vector[m],range=(range1,range2))
        indexes = peakutils.indexes(n, thres=0.02/max(smooth(n,20)), min_dist= 150*bin_vector[m]/(range2-range1))

        peaks= peakutils.interpolate(np.linspace(range1,range2,bin_vector[m]), n, ind=indexes, width=1) #uses a gaussian function and the precedent indexes to enhance our peak finding

        print len(n)
        print len(bins)
        print len(patches)
        peaks = [double(i) for i in peaks]
        #M= np.linspace(1,len(peaks)-1,len(peaks)-1)
        
        #print len(peaks)
        dist_peaks=np.diff(peaks)
        #dist_peaks=[peaks[i+1]-peaks[i] for i in range(len(peaks)-1)]
        #print len(dist_peaks)
        #print np.mean(dist_peaks)
        x = np.linspace(1,len(peaks),len(peaks))
        range= np.linspace(range1,range2,len(peaks))
        slope = np.polyfit(x,peaks,1)
        x1=np.linspace(range1,range2)
        y1=slope[0]*x+slope[1]
        media=np.mean(dist_peaks)
        delta[m]=abs(slope[0]-media)
        #print "slope = "+str(slope[0])
        #print "media = "+str(media)
        if m==0:
	    factor=slope[0]
            factor2=np.mean(dist_peaks)
	    delta2=delta[m]

 
        else:
            if delta[m]<delta2:
                factor=slope[0]
                factor2=np.mean(dist_peaks)
		delta2=delta[m]
	    	M=m
	    	
	m=m+1
    print str(factor)+" delta = "+str(delta2)
#    plt.plot(np.linspace(range1,range2,bin_vector[M]),smooth(n,5))
    # print n[indexes]
    #plt.plot(peaks,n[indexes],'bo')
    #plt.show()
#    plt.figure(3)
#    plt.plot(x,peaks,'.')
#    plt.plot(x,y1)

#    plt.figure(4)
#    plt.hist(dist_peaks)
