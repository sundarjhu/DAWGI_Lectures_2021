from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('classic')

def list_ticks(x):
    x_tk=[]
    for i in x:
        if i%1.==0.:
            x_tk.append(str(int(i)))
        else:
            x_tk.append(str(i))
            
    return x_tk

Zw = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.020, 0.030, 0.040]
N_zw = len(Zw)
color_map = plt.cm.gist_rainbow(np.linspace(0., 1., N_zw))
colors={ Zw[i] : color_map[i] for i in range(N_zw)}


################################################################################
################################################################################

def sfh_plot_mode(name, Z0, age0, SFR0, typ, fig):

    Z = np.unique(Z0)
    Nz=len(Z)
    idx_Z = range(1,Nz+1)

    ages = np.unique(age0)
    age_list = list_ticks(np.round(ages,1))

    Nag=len(ages)
    ages_aux=np.arange(1,Nag+1)
    
    niso = Nz*Nag

    ### plot ###
    ax = fig.add_subplot(131, projection='3d',autoscale_on=True)

    nn=1
    for zn in Z:
        sfr = SFR0[np.where(Z0==zn)]
        cs = [colors[zn]] * len(ages)
        plt.bar(ages_aux-0.2, sfr, width=0.3, zs=nn, zdir='x', align='center', color=cs, alpha=0.8, linewidth=0)
        nn = nn + 1

    ax.view_init(30, -135)
    sz=10
    ax.set_xlabel('Z')
    ax.set_ylabel('Age (Gyr)')
    ax.set_zlabel('Stellar Fraction')


    idx_Z = np.arange(1,4)
    plt.xticks(idx_Z, ('0.014', '0.017', '0.020'))
    plt.yticks(ages_aux, ('0.03','0.06','0.10','0.18','0.32','0.56','1.00'))


    plt.xlim(min(idx_Z)-0.5,max(idx_Z)+0.5)
    
    
    plt.ylim(min(ages_aux)-0.5,max(ages_aux)+0.5)

    ax.set_zlim(0., 1.)

    plt.tight_layout()




