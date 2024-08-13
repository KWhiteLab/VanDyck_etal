import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from scipy import stats

# style stuff for plotting
font_names = [f.name for f in fm.fontManager.ttflist]
# mpl.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
params = {'mathtext.default': 'regular'}          
plt.rcParams.update(params)

## Functions ##
def plot_distance(ph,thresh):
    plt.figure(figsize = (8,5))
    p = 0
    for file in files:
        if float(phs[p]) == ph:
            l = 0
            curve = []
            frames = []
            with open(file,'r') as f:
                for line in f:
                    linestr = line.split()
                    
                    if l > 0: 
                        frames.append(float(linestr[0]))
                        curve.append(float(linestr[1]))
                            
                    l = l + 1
            
            frames = [2*ele/100 for ele in frames]
            
            ## truncate simulation based on threshold
            frames_trunc = [ele for ele in frames[int(thresh*len(frames))-1:]]
            curve_trunc = [ele for ele in curve[int(thresh*len(frames))-1:]]
            
            fd_bins = freedman_diaconis(data=curve, returnas="bins")
            
            plt.hist(curve_trunc,bins=10,density=True)
            plt.axvline(x=np.mean(curve_trunc),color='black',linewidth=3,ymax=0.9)
            t = plt.text(0.87, 0.9,'pH = ' + str(ph), horizontalalignment='center',
                verticalalignment='center',
                transform=plt.gca().transAxes)
            t.set_bbox(dict(facecolor='yellow', alpha=0.75, edgecolor='yellow'))
            plt.xticks([10,20,30,40,50,60])
            plt.ylabel('Density')
            plt.xlabel(r'NSH2:PTP Distance ($\AA$)')
            plt.tight_layout()
            plt.show()
            
            stat,pval = stats.shapiro(curve_trunc)
            print('pH ' + str(ph) + ' p-value = ',pval)
        
        p += 1

def freedman_diaconis(data, returnas="width"):
    """
    Use Freedman Diaconis rule to compute optimal histogram bin width. 
    ``returnas`` can be one of "width" or "bins", indicating whether
    the bin width or number of bins should be returned respectively. 

    Parameters
    ----------
    data: np.ndarray
        One-dimensional array.

    returnas: {"width", "bins"}
        If "width", return the estimated width for each histogram bin. 
        If "bins", return the number of bins suggested by rule.
    """
    data = np.asarray(data, dtype=np.float_)
    IQR  = stats.iqr(data, rng=(25, 75), scale="raw", nan_policy="omit")
    N    = data.size
    bw   = (2 * IQR) / np.power(N, 1/3)

    if returnas=="width":
        result = bw
    else:
        datmin, datmax = data.min(), data.max()
        datrng = datmax - datmin
        result = int((datrng / bw) + 1)
    return(result)

## set variables and call functions ##
phs = ['4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','8.5']

# list of cpptraj data files to analyze #
# Format:
## frame1  value1 ##
## frame2  value2 ##
files = ['NH2_PTP_dist_pH' + ph + '0.dat' for ph in phs]

sample_phs = [4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5]
threshold = 0.5

for ph in sample_phs:
    plot_distance(ph,threshold)
