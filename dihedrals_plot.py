import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import numpy as np

# style stuff for plotting
font_names = [f.name for f in fm.fontManager.ttflist]
# mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

phs = [4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5]

# list of cpptraj data files to analyze #
# Format:
## frame1  value1 ##
## frame2  value2 ##
files1 = ['H201_chi1_pH' + str(ph) + '.dat' for ph in phs]

def plot_curves(files,mark,c):
    dihedrals = []
    frames = []
    for f in files:
        l = 0
        data_temp = []
        with open(f,'r') as exp:
            for line in exp:
                linestr = line.split()
                if l > 0:
                    data_temp.append(float(linestr[1]))
                    
                l = l + 1
                
        dihedrals.append(data_temp)
        
    dih_avgs = []
    dih_std = []

    for i in range(len(dihedrals)):
        dih_avgs.append(np.mean(dihedrals[i][350:]))
        dih_std.append(np.std(dihedrals[i][350:]))

    plt.plot(phs,dih_avgs,marker = mark,markersize=10,color = c)

    plt.fill_between(phs,
                        [dih_avgs[i] + dih_std[i]/2 for i in range(len(dih_avgs))],
                        [dih_avgs[i] - dih_std[i]/2 for i in range(len(dih_avgs))],
                        alpha = 0.2,label='_nolegend_',color=c)

plot_curves(files1,'s','#ff7f00')

plt.xlabel('pH')
plt.ylabel('Diherdral angle (degrees)')
plt.tight_layout()
plt.show()