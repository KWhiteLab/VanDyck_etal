import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import scipy.stats as st

# style stuff for plotting
font_names = [f.name for f in fm.fontManager.ttflist]
# mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['legend.title_fontsize'] = '20'

colors = plt.cm.viridis([0,0.175,0.32,0.45,0.6,0.75,0.83,0.9,1])
fs = 14

def plot_distance(window,files):
    end_avg = []
    end_ci = []

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

        m, c = mean_confidence_interval(curve[window[0]:window[1]])
        end_avg.append(m)
        end_ci.append(c)

    return end_avg, end_ci

def plot_triangle_bisector(a, b, c, col):
    """
    Plots a triangle with sides a, b, and c, with vertex A fixed at (0,0)
    and vertex B initially on the x-axis at (c, 0). Then it rotates the
    entire triangle so that the angle bisector at A lies on the x-axis.
    
    The sides correspond to: AC = a, BC = b, and AB = c.
    """
    # Validate the triangle inequality
    if a + b <= c or a + c <= b or b + c <= a:
        raise ValueError("The provided side lengths do not form a valid triangle.")
    
    # 1. Place vertices A, B, C in an initial position
    A = np.array([0.0, 0.0])
    B = np.array([c,   0.0])  # side AB = c
    # Law of cosines for angle at A
    cos_theta = (a**2 + c**2 - b**2) / (2 * a * c)
    cos_theta = np.clip(cos_theta, -1, 1)  # avoid floating precision issues
    sin_theta = np.sqrt(1 - cos_theta**2)
    # side AC = a
    C = np.array([a * cos_theta, a * sin_theta])

    # 2. Angle at A is theta = angle of AC from the x-axis
    #    B is on x-axis, so angle(AB) = 0
    #    angle(AC) = arctan2(y, x) = arctan2(a*sin_theta, a*cos_theta)
    theta = np.arctan2(C[1], C[0])  # same as sin_theta, cos_theta

    # 3. Rotate all points so the angle bisector is along the x-axis
    #    The angle bisector is at theta/2 from AB. 
    #    So we rotate by -theta/2.
    alpha = -theta / 2.0
    rot = np.array([
        [np.cos(alpha), -np.sin(alpha)],
        [np.sin(alpha),  np.cos(alpha)]
    ])
    A_rot = A  # A is the origin, unaffected by rotation
    B_rot = rot.dot(B)
    C_rot = rot.dot(C)

    # 4. Plot the triangle
    # Create array of rotated vertices for a single polyline
    vertices = np.array([A_rot, B_rot, C_rot, A_rot])
    
    plt.plot(vertices[:, 0], vertices[:, 1], '-', color=col,linewidth=2)
    plt.scatter(A_rot[0], A_rot[1], marker='o', s=220, color=col, label='_nolegend_')
    plt.scatter(B_rot[0], B_rot[1], marker='s', s=220, color=col, label='_nolegend_')
    plt.scatter(C_rot[0], C_rot[1], marker='^', s=220, color=col, label='_nolegend_')

    # Optional: you can adjust x/y limits or just let matplotlib auto-scale
    # For example, to give some margin:
    all_x = vertices[:, 0]
    all_y = vertices[:, 1]
    x_min, x_max = all_x.min(), all_x.max()
    y_min, y_max = all_y.min(), all_y.max()
    margin = 0.2 * max(x_max - x_min, y_max - y_min)

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m = np.mean(a)
    se = st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    return m, h

# -----------------------------

## set variables and call functions ##
window = [5,8]  # define a window in nanoseconds to take the average domain distance from 
window = [int(100*ele/2) for ele in window]
colors_plot = ['cornflowerblue','green','magenta']
path = '/path/'
mutants = ['SHP-2 WT','SHP-2 HAEA','SHP-2 H116A','SHP-2 E252A','SHP-2 E249A']
prefixes = [['2shp_monomer/analysis/'],['2shp_HAEA/analysis/'],['2shp_H116A/analysis/'],['2shp_E252A/analysis/'],['2shp_E249A/analysis/']]
for m, prefix in enumerate(prefixes):
    phs = ['4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','8.5']
    phs = phs[::-1] 
    files = []
    pairs = ['CSH2_NSH2','NSH2_PTP','CSH2_PTP']
    markers = ['o','s','^']
    distances = []
    distances_e = []
    for p in pairs:
        distance = []
        distance_e = []
        for i in range(len(prefix)):
            files_list = [path + prefix[i] + p + '_dist_pH' + ph + '0.dat' for ph in phs]
            mean, err = plot_distance(window,files_list)
            distance.append(mean)
            distance_e.append(err)
        distances.append(distance)
        distances_e.append(distance_e)

    plt.figure(figsize = (6,6))
    rotation_angle = 0
    for i in range(len(distance[0])):
        plot_triangle_bisector(distances[1][0][i],distances[0][0][i],distances[2][0][i],colors[::-1][i])
        
    plt.axis('equal')
    plt.title(mutants[m])
    plt.show()
    plt.clf()

    fig, axs = plt.subplots(1,3, figsize=(10,4))
    fig.suptitle(mutants[m])

    pairs = [[1,2],[1,0],[2,0]]
    xlabels = ['NSH2:PTP','NSH2:PTP','CSH2:PTP']
    ylabels = ['CSH2:PTP','CSH2:NSH2','CSH2:NSH2']
    
    dists_pdb_x = [31.3358,31.3358,41.3063]
    dists_pdb_y = [41.3063,36.5595,36.5595]
    xl = 0
    yl = 0
    xlims = [[28,66],[28,66],[45,54]]
    ylims = [[45,54],[35,44],[35,44]]
    xtick = [[30,40,50,60],[30,40,50,60],[46,48,50,52,54]]
    ytick = [[46,48,50,52,54],[36,38,40,42,44],[36,38,40,42,44]]
    for pair in pairs:
        for p in range(len(phs)):
            x = distances[pair[0]][0][p]
            xe = distances_e[pair[0]][0][p]
            y = distances[pair[1]][0][p]
            ye = distances_e[pair[1]][0][p]
            axs[xl].errorbar(x,y,xerr=xe,yerr=ye,capsize=2,fmt='o',markersize=2,color=colors[::-1][p])

        axs[xl].set_xlabel(xlabels[xl] + ' Distance ($\AA$)')
        axs[xl].set_ylabel(ylabels[yl] + ' Distance ($\AA$)')
        # axs[xl].axhline(dists_pdb_y[xl],color='grey',linestyle='--')
        # axs[xl].axvline(dists_pdb_x[xl],color='grey',linestyle='--')
        # axs[xl].scatter(dists_pdb_x[xl],dists_pdb_y[xl],facecolors='none',edgecolors='r')
        # axs[xl].set_xlim(xlims[xl])
        # axs[xl].set_ylim(ylims[xl])
        # axs[xl].set_xticks(xtick[xl])
        # axs[xl].set_yticks(ytick[xl])
        # axs[xl].legend(phs[::-1],frameon=False)
        xl = xl+1
        yl = yl+1
    plt.tight_layout()
    plt.show()
