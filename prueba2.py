from ROOT import *

import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks
# from scipy.signal import find_peaks_cw
import peakutils
import sys
import time
inroot = sys.argv[1]

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth



f = TFile(inroot)
skPixTree = f.Get('skPixTree')

skPixTree.SetBranchStatus("pix",0)
skPixTree.SetBranchStatus('pix',1)

for i in range(1):
    pixValue = []
    x = []
    # print x
    for j in range(19999):
        x = np.append(x,j)
        #print x
    #print 1
    for index, event in enumerate(skPixTree):
        #
        #  if index==10000:
        #print i, index
        #     break
        #
        #     # time.sleep(2)
        #
        #  # if  (40000*i > index and index >= 1000 + 40000*(i-1)):
        #     # print index
        #     # time.sleep(5)
        #
        #     # continue
        #  # print index
        pixValue = np.append(pixValue,event.pix)
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

    range1=250000
    range2=550000
    bin=(range2-range1)/100
    n,bins,patches = plt.hist(pixValue,bins=bin,range=(range1,range2))
    # plt.figure(1)
    #:plt.axis([-300,10000,0,1000
    indexes = peakutils.indexes(smooth(n,5), thres=0.02/max(smooth(n,5)), min_dist= bin)
    peaks= peakutils.interpolate(np.linspace(range1,range2,bin), smooth(n,5), ind=indexes, width=1) #uses a gaussian function and the precedent indexes to enhance our peak finding
    print indexes
    print peaks
    peaks = [int(i) for i in peaks]
    print peaks
         # plt.show()
    plt.figure(2)

    #plt.plot(n)
    plt.plot(np.linspace(range1,range2,bin),smooth(n,5))
    # print n[indexes]
    plt.plot(peaks,n[indexes],'bo')
    plt.show()
    i+=1
