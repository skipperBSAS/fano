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

    n,bins,patches = plt.hist(pixValue,bins=2000,range=(16000,18000))
    # plt.figure(1)
    #:plt.axis([-300,10000,0,1000])

    indexes = peakutils.indexes(smooth(n,20), thres=0.02/max(smooth(n,20)), min_dist=200)

    print indexes
         # plt.show()
    plt.figure(2)

    plt.plot(n)
    plt.plot(smooth(n,10))
    plt.plot(indexes,n[indexes],'bo')
    plt.show()
    i+=1
