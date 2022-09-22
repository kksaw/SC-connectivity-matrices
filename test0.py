# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

#frame rate is 10Hz

def importDF(filename):
    with open(filename,'r') as f:
        lines = f.readlines()

    n =  len(lines[0].split('\t'))    
    mylist=[]
    for i in range(n):
        mylist.append([])
        
        for x in lines:
            mylist[i].append(x.split('\t')[i])
    
    mylist[-1] = [x.rstrip() for x in mylist[-1]]       	
    
    return np.asarray(mylist,float)

spikearray = importDF('df\df1.txt')

'''random scripts'''
rougharray = np.zeros((79,600),float)

for i in range(79):
    for j in range(600):
        rougharray[i,j]=float(spikearray[i,j*10])
''''''

'''spike-time correlation'''

'''power based connectivity aka paper?'''
#can use scipy.stats.pearsonr (returns r and pvalue)
from scipy import special, stats
def calc_correlation(sarray):
    nN = len(sarray)
    pmat = np.zeros((nN,nN),float)
    smat = np.zeros((nN,nN),float)

    for i in range(nN):
        for j in range(nN):
            if i < j:
                pmat[i,j]=pearson_corr(sarray[i],sarray[j])[0]
                smat[i,j]=stats.spearmanr(sarray[i],sarray[j])[0]
    return pmat,smat

def calc_corr(sarray, swindow):
    nN, nT = np.shape(sarray)
    nW = nT/swindow if nT%swindow==0 else nT/swindow+1

    pmat = np.zeros((nW,nN,nN),float)
    smat = np.zeros((nW,nN,nN),float)
    mmat = np.zeros((nW,nN,nN),float)
    gmat = np.zeros((nW,nN,nN),float)
    
    for i in range(nW):
        newdata = np.zeros((nN,swindow),float) if i!=nW-1 else np.zeros((nN,(nT-(nW-1)*swindow)),float)
        for n in range(nN):
            newdata[n] = sarray[n, i*swindow:(i+1)*swindow] if i!=nW-1 else sarray[n*swindow::]
    
        for a in range(nN):
            for b in range(nN):
                if a < b:
                    pmat[i,a,b]=stats.peaersonr(newdata[a],newdata[b])[0]
                    smat[i,a,b]=stats.spearmanr(newdata[a],newdata[b])[0]
                    mmat[i,a,b]=mutual_info_score(newdata[a],newdata[b])
                    gmat[i,a,b]=Gcausality(newdata[a],newdata[b])[0]
                elif a > b:
                    gmat[i,a,b]=Gcausality(newdata[a],newdata[b])[1]
                    
    return pmat, smat, mmat, gmat
        
'''                   
def pearson_corr(a1,a2):
    x=a1
    y=a2
    nt = len(x)
    xm, ym = x-x.mean(), y-y.mean()
    r_num = np.sum(xm*ym)
    r_den = np.sqrt(np.sum(xm*xm)*np.sum(ym*ym))
    r = r_num/r_den
    
    r = max(min(r,1.0),-1.0)
    df = nt-2
    if abs(r)==1.0:  #floating point arithmetic or same neuron
        prob=0.0
    else:
        t_squared = r*r * (df / ((1.0 - r) * (1.0 + r)))
        prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))   #assumes df / (df + t_squared) always < 1.0. if not check scipy.stats.pearson
    return r, prob

def spearman_corr(a1,a2):
    return
'''

'''phase-based - get phase through FFT?'''
#FFT then convolve - hmmm
from scipy import signal, fftpack
def calc_phase(sarray):
    nN = len(sarray)
    imat = np.zeros((nN,nN),float)
    
    for i in range(nN):
        for j in range(nN):
            if i < j:
                imat[i,j]=FFTconv(sarray[i],sarray[j])
    return imat
    
def FFTconv(a1,a2,wlength=10):
    nT = len(a1)
    win = signal.hann(wlength)
    signals=[a1,a2]
    phase=np.zeros((2,nT//2),float)
    
    for i in range(2):
        f = fftpack.fft(signals[i])
        p = np.angle(f[0:nT//2])
        pF = signal.convolve(p,win,mode='same')/sum(win)
        phase[i]=pF

    ispc = np.abs(np.mean(np.exp(1j*(phase[0]-phase[1]))))

    return ispc

def main():
    imat = calc_phase(spikearray)
       
'''
    N = 6000     #length of signal (number of timesteps)
    T = 1.0/800.0   #800 is sampling rate. 1 is time
    x = np.linspace(0.0,N*T,N)
    xf = np.linspace(0.0,1.0/(2/0*T),N//2)
    magnitude = 2.0/N*np.abs(yf[0:N//2])    #ignore -ve values, 2x +ve values, normalise w sample points
    phase = 2.0/N*np.angle(yf[0:N//2])
    
    - somehow phase angles always the same at 1200 and 2400... ???
     
    x = np.arange(nT//2)   #floor
    y = fftpack.fft(array1)
    power = 2.0/nT*np.abs(y[0:nT//2])
    phase = np.angle(y[0:nT//2])
    
    powerF = signal.convolve(power,win,mode='same')/sum(win)
    phaseF = signal.convolve(phase,win,mode='same')/sum(win)

'''

'''scipy correlate'''
#correlation at each time point btw 2 signals
corr = signal.signaltools.correlate(spike1,spike2,mode='same')


'''power-based connectivity'''


'''Granger prediction'''
#only 1 trial..?
from numpy import linalg as LA
from scipy import stats, signal

def main(sarray):
    nN = len(sarray)
    gmat = np.zeros((nN,nN),float)
    
    for i in range(nN):
        for j in range(nN):
            if i < j:
                gmat[i,j]=Gcausality(sarray[i],sarray[j])[0]
            elif i > j:
                gmat[i,j]=Gcausality(sarray[i],sarray[j])[1]
    
    return gmat
    
def Gcausality(a1,a2, ntrial=1, orderpt=7):
    #only need E
    #regularity??? set TH with sinusoidal + noise
    #set order???
    nT = len(a1)
    tempdata = np.zeros((2,nT))
    
    tempdata[0] = stats.zscore(a1, ddof=1)
    tempdata[1] = stats.zscore(a2, ddof=1)
    
        Ex = armorf(np.array([tempdata[0]]),ntrial,nT,orderpt)
        Ey = armorf(np.array([tempdata[1]]),ntrial,nT,orderpt)
        E  = armorf(tempdata,ntrial,nT,orderpt)
    
        y2x = np.log(Ex/E[0,0])
        x2y = np.log(Ey/E[1,1])
    
#        for bici in range(np.shape(bic)[1]):
#           E = armorf(tempdata,ntrial,timept,bici+1)[1]    
#          bic[i,bici] = np.log(LA.det(E)) + (np.log(np.size(tempdata))*(bici+1)*2**2)/np.size(tempdata)
        
    return y2x, x2y

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
    
    for m in range(p):
        ck = np.dot(LA.inv(LA.cholesky(pf).conj().transpose()), pfb.dot(LA.inv(LA.cholesky(pb))))
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
    
    error = np.dot(En, En.conj().transpose())
    
    return error

def is_stationaryADF(test, sigf='5%'):
    from statsmodels.tsa.stattools import adfuller
    
    result = adfuller(test)
    adfstat = result[0]
    cvalue = result[4][sigf]
    
    return True if adfstat < cvalue else False

            
'''Mutual information'''
#cheat method
#from sklearn.metrics import mutual_info_score

#long method (sklearn function - but ignores labelling hmm)
#how to choose bins, default=10
def mutual_info_score(x,y,bins=10):
    cxy= np.histogram2d(x,y,bins)[0]
    nzx, nzy = np.nonzero(cxy)
    nz_val = cxy[nzx,nzy]
    
    csum = cxy.sum()
    pi = np.ravel(cxy.sum(axis=1))
    pj = np.ravel(cxy.sum(axis=0))
    log_cnm = np.log(nz_val)
    cnm = nz_val/csum
    outer = pi.take(nzx).astype(np.int64)*pj.take(nzy).astype(np.int64)
    log_outer = -np.log(outer) + np.log(pi.sum()) + np.log(pj.sum())
    mi = cnm*(log_cnm - np.log(csum)) + cnm*log_outer
    
    return mi.sum()

def main(sarray):
    nN = len(sarray)
    mmat = np.zeros((nN,nN),float)
    
    for i in range(nN):
        for j in range(nN):
            if i < j:
                mmat[i,j]=mutual_info_score(sarray[i],sarray[j])
    
    return mmat

    

'''plotting script'''
#plot polar plot
fig = plt.figure()
ax = plt.subplot(111,projection='polar')
for x in a:     #a is array of complex numbers
    ax.polar([0,np.angle(x)],[0,np.abs(x)])
    
#rotate 3d plot
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(xf, magnitude, phase)
ax.set_xlabel('frequency')
ax.set_ylabel('power')
ax.set_zlabel('phase')

plt.show()

#imshow
from matplotlib import colors
class MptNorm(colors.Normalize):
	def __init__(self, vmin=None, vmax=None, vmid=None, clip=False):
		self.vmid = vmid
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.vmid, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

cval = 0.01
fig, ax = plt.subplots()
cax = ax.imshow(mmat, cmap='Reds', clim=(-cval,cval), 
                norm=MptNorm(vmin=-cval, vmax=cval, vmid=0))
cbar = fig.colorbar(cax)





#signal
corr = scipy.signal.correlate2d(spikearray, spikearray)   #hmm

'''min daily temp'''
from pandas import Series, DataFrame, concat
from pandas.tools.plotting import lag_plot, autocorrelation_plot
from statsmodels.graphics.tsaplots import plot_acf
from sklearn.metrics import mean_squared_error
from statsmodels.tsa.ar_model import AR

series = Series.from_csv('daily-minimum-temperatures.csv', header=0)
lag_plot(series)    #check for correlation

#create lagged dataset (1 day lag)
values = DataFrame(series.values)
dataframe = concat([values.shift(1), values], axis=1)
dataframe.columns = ['t-1','t']
result = dataframe.corr()

autocorrelation_plot(series)    #check autocorrelation across diff lags
plot_acf(series, lags=31)   #another way of checking correlation

#split into train and test sets (1 day lag)
X = dataframe.values
train, test = X[1:len(X)-7], X[len(X)-7:]
train_X, train_y = train[:,0], train[:,1]
test_X, test_y = test[:,0], test[:,1]

test_score = mean_squared_error(test_y,test_X)

#autoregression
X = series.values
train, test = X[1:len(X)-7], X[len(X)-7:]
model = AR(train)
model_fit = model.fit()

predictions = model_fit.predict(start=len(train), end=len(train)+len(test)-1, dynamic=False)
error = mean_squared_error(test, predictions)



