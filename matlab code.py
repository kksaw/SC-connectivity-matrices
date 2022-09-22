# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 15:26:04 2018

@author: khaik
"""
'''
1. optimal order (bic)
'''

import numpy as np
from numpy import linalg as LA
from scipy.io import loadmat
from scipy import stats, signal
from matplotlib import pyplot as plt

def loadEEG(filename):
    mdata = loadmat(filename)['EEG']
    mdtype = mdata.dtype
    ndata = {n: mdata[n][0,0] for n in mdtype.names}
    return ndata['data']

def Gcausality():
    #only need E, find better function?
    #only 1 trial for each time point for Ca signalling...
    data = loadEEG('code\sampleEEGdata.mat')
    
    srate = 256
    timewin = 200
    order = 27
    ntrials = 99
    
    timewin_points = round(timewin/(1000.0/srate))
    order_points = round(order/(1000.0/srate))
    
    eegdata = np.zeros((2,640,ntrials))
    for i in range(640):    
        eegdata[0,i] = data[26][i] - np.mean(data[26][i])
        eegdata[1,i] = data[5][i] - np.mean(data[5][i])
    
    times2saveid = np.asarray([155,160,165,170,175,180,185,190,196,201,206,211,216,221,226,231,
                    237,242,247,252,257,262,267,272,277,283,288,293,298,303,308,313,
                    318,324,329,334,339,344,349,354,359,365,370,375,380,385,390,395,
                    400,405,411,416,421,426,431,436,441,446,452,457,462,467,472,477,
                    482,487,493,498,503,508,513])
    
    x2y = np.zeros(len(times2saveid))
    y2x = np.zeros(len(times2saveid))
    bic = np.zeros((len(times2saveid),15))
        
    for timei in range(len(times2saveid)):
        temparray = np.squeeze(np.copy(eegdata[:,times2saveid[timei]-int(np.ceil(timewin_points/2.0)):times2saveid[timei]+int(timewin_points/2.0)-int(np.remainder(timewin_points+1.0,2.0)),:]))
        
        for triali in range(len(temparray)):
            temparray[0,:,triali] = stats.zscore(signal.detrend(np.squeeze(np.copy(temparray[0,:,triali]))),ddof=1)
            temparray[1,:,triali] = stats.zscore(signal.detrend(np.squeeze(np.copy(temparray[1,:,triali]))),ddof=1)   
        
        '''check for stationarity'''
        #if not is_stationaryADF(temparray):
        #    break
            
        tempdata = np.zeros((2,np.size(temparray[0])))
        for i in range(2):
            tempdata[i] = np.concatenate(temparray[i].transpose())
        
        Ex = armorf(np.array([tempdata[0]]),ntrials,timewin_points,order_points)[1]
        Ey = armorf(np.array([tempdata[1]]),ntrials,timewin_points,order_points)[1]
        E = armorf(tempdata,ntrials,timewin_points,order_points)[1]
        
        y2x[timei] = np.log(Ex/E[0,0])
        x2y[timei] = np.log(Ey/E[1,1])
        
        for bici in range(np.shape(bic)[1]):
            E = armorf(tempdata,ntrials,timewin_points,bici+1)[1]    
            bic[timei,bici] = np.log(LA.det(E)) + (np.log(np.size(tempdata))*(bici+1)*2**2)/np.size(tempdata)
        
        return x2y, y2x, bic

def test_sigF():
    p-value = stats.levene()[1]
    return
def swindow():
    #O.o    
    return


def armorf(x, Nr, Nl, p):
    '''
    x - 2d array, each row is one variable's time series
    Nr - number of realisations
    Nl - length of every realisation
    p - order
    
    A = armorf(...) returns AR coefficients
    [A,E] = armorf(...) returns prediction error (covariance matrix)
    [A,E,K] = armorf(...) returns vector K of reflection coeff (parcor coeff)
    
    armorf([1x5049],99,51,7)
    armorf([2x5049],99,51,7)
    '''
    #if x in 1D
    #x = np.array([x])
    
    [L,N] = np.shape(x)
    pf = np.zeros((L,L))
    pb = np.zeros((L,L))
    pfb = np.zeros((L,L))
    ap = np.zeros((L,L,1))
    bp = np.zeros((L,L,1))
    En = np.zeros((L,L))
    
    for i in range(Nr):     #99
        En += np.dot(x[:,i*Nl:(i+1)*Nl], x[:,i*Nl:(i+1)*Nl].conj().transpose())
        ap[:,:,0] += np.dot(x[:,i*Nl+1:(i+1)*Nl], x[:,i*Nl+1:(i+1)*Nl].conj().transpose())
        bp[:,:,0] += np.dot(x[:,i*Nl:(i+1)*Nl-1], x[:,i*Nl:(i+1)*Nl-1].conj().transpose())
    
    
    ap[:,:,0] = LA.inv((LA.cholesky(ap[:,:,0]/Nr*(Nl-1))).conj())
    bp[:,:,0] = LA.inv((LA.cholesky(bp[:,:,0]/Nr*(Nl-1))).conj())
    
    for i in range(Nr):
        efp = np.dot(ap[:,:,0], x[:,i*Nl+1:(i+1)*Nl])
        ebp = np.dot(bp[:,:,0], x[:,i*Nl:(i+1)*Nl-1])
        pf += np.dot(efp, efp.conj().transpose())
        pb += np.dot(ebp, ebp.conj().transpose())
        pfb += np.dot(efp, ebp.conj().transpose())
    
    En = LA.cholesky(En/N).conj()
    
    coeff = []
    kr = []
    
    for m in range(p):
        ck = np.dot(LA.inv(LA.cholesky(pf).conj().transpose()), pfb.dot(LA.inv(LA.cholesky(pb))))
        kr = [kr,ck]
        ef = np.eye(L) - np.dot(ck, ck.conj().transpose())
        eb = np.eye(L) - np.dot(ck.conj().transpose(), ck)
        
        En = np.dot(En, LA.cholesky(ef).conj())
        E = (ef+eb)/2   #what is this
        
        ap = np.dstack((ap,np.zeros((L,L))))
        bp = np.dstack((ap,np.zeros((L,L))))
        
        pf = np.zeros((L,L))
        pb = np.zeros((L,L))
        pfb = np.zeros((L,L))
    
        a = np.zeros((L,L,m+2))
        b = np.zeros((L,L,m+2))        
        for i in range(m+2):
            a[:,:,i] = np.dot(LA.inv(LA.cholesky(ef).conj()), (ap[:,:,i]-ck.dot(bp[:,:,m+1-i])))
            b[:,:,i] = np.dot(LA.inv(LA.cholesky(eb).conj()), (bp[:,:,i]-ck.dot(ap[:,:,m+1-i])))
        
        for k in range(Nr):
            efp = np.zeros((L,Nl-m-2))
            ebp = np.zeros((L,Nl-m-2))
           
            for i in range(m+2):
                k1 = m-i+k*Nl+2
                k2 = Nl-i+k*Nl
                efp += np.dot(a[:,:,i], x[:,k1:k2])
                ebp += np.dot(b[:,:,m+1-i], x[:,k1-1:k2-1])
            
            pf += np.dot(efp, efp.conj().transpose())
            pb += np.dot(ebp, ebp.conj().transpose())
            pfb += np.dot(efp, ebp.conj().transpose())
        
        ap = a
        bp = b
    
    for j in range(p):
        coeff.append(np.dot(LA.inv(a[:,:,0]), a[:,:,j+1]))
    
    error = np.dot(En, En.conj().transpose())
    
    return coeff, error, kr    

def is_stationaryADF(test, sigf='5%'):
    from statsmodels.tsa.stattools import adfuller
    
    result = adfuller(test)
    adfstat = result[0]
    cvalue = result[4][sigf]
    
    return True if adfstat < cvalue else False



if __name__ == '__main__':
    main()