import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# style stuff for plotting
font_names = [f.name for f in fm.fontManager.ttflist]
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['legend.title_fontsize'] = '20'

phs = ['4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','8.5']

# list of cpptraj data files to analyze #
# Format:
## frame1  value1 ##
## frame2  value2 ##
files = ['CSH2_NSH2_dist_pH' + ph + '0_v2.dat' for ph in phs]

viridis = mpl.colormaps['viridis'].resampled(8)
colors = [viridis(ele/9) for ele in range(9)]
fs = 27
def plot_distance(window):
    end_avg = []
    end_std = []
    plt.figure(figsize = (6,5))
    c = 0
    for file in files:
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

        end_avg.append(np.mean(curve[window[0]:window[1]]))
        end_std.append(np.std(curve[window[0]:window[1]]))
        plt.plot(frames,curve,alpha = 0.8,linewidth = 4, color = colors[c])

        c = c + 1

    plt.legend(phs,title = 'pH',bbox_to_anchor = (1.05,1),fontsize=fs)
    plt.xlabel('Time (ns)',fontsize = fs)
    plt.ylabel(r'Distance ($\AA$)',fontsize = fs)
    plt.axvspan(2*window[0]/100,2*window[1]/100,alpha = 0.5)
    plt.xticks([0,1,2,3,4,5,6,7,8],fontsize = fs)
    plt.yticks(fontsize = fs)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize = (6,5))

    plt.plot([float(ele) for ele in phs],end_avg,linewidth = 4,marker='o')
    plt.fill_between([float(ele) for ele in phs],
                     [end_avg[i] + end_std[i]/2 for i in range(len(end_avg))],
                     [end_avg[i] - end_std[i] for i in range(len(end_avg))],
                     alpha = 0.5)

    plt.xlabel('pH',fontsize = fs)
    plt.ylabel(r'Average Distance ($\AA$)',fontsize = fs)
    plt.tight_layout()
    plt.xticks([4,5,6,7,8,9],fontsize = fs)
    plt.yticks(fontsize = fs)
    plt.tight_layout()
    plt.show()

## set variables and call functions ##
window = [5,8]  # define a window in nanoseconds to take the average domain distance from 
window = [int(100*ele/2) for ele in window]

plot_distance(window)
