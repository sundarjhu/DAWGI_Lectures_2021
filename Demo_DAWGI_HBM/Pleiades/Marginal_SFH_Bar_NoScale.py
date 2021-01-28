import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import a_statistics_def_fun as st_def
plt.style.use('classic')


def list_ticks(x):
    x_tk=[]
    for i in x:
        if i%1.==0.:
            x_tk.append(str(int(i)))
        else:
            x_tk.append(str(i))
            
    return x_tk

##################

def marg_Z(ax, nag, n_z):
    a_mar = []
    x1, x2 = 0, nag
    for i in range(n_z):
        a_mar.append(np.sum(ax[x1:x2]))
        x1+=nag
        x2+=nag
    return np.array(a_mar)


def marg_AGE(ax, nag):
    a_mar = []
    x1, x2 = 0, nag
    for i in range(nag):
        a_mar.append(np.sum(ax[i::nag]))
    return np.array(a_mar)


################################################################################
################################################################################
################################################################################

def marg_sfh_bar_age(name,sfh,a_sp,fig):

    Z0, age0, mode = sfh[0], sfh[1], sfh[4]

    Z = np.unique(Z0)
    Nz=len(Z)
    idx_Z = range(1,Nz+1)

    age = np.unique(age0)
    age_list = list_ticks(np.round(age,1))
    Nag=len(age)
    age_aux= np.arange(1,Nag+1)

    SFR_mode_marg = marg_AGE(mode, Nag)
    
    ##
    a_age = []
    cont=0
    for ai in a_sp:
        a_aux = marg_AGE(ai, Nag)
        a_age.append(a_aux)

    a_age = np.array(a_age)

    perc = st_def.a_stat(a_age.T)
    ##

    sfh_mgl=[]
    for i in range(Nag):
        sfh_mgl.append([age_aux[i], SFR_mode_marg[i], perc[i][0], perc[i][1], perc[i][2]])
    sfh_mgl = np.array(sfh_mgl)
    
    sfh_mgl = sfh_mgl.T
    

    
    ###########################################

    ax = fig.add_subplot(132)

    violin_parts = ax.violinplot(a_age, positions=age_aux, showmedians=True)
    
    for partname in ('cbars','cmins','cmaxes','cmedians'):
        vp = violin_parts[partname]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    # Make the violin body blue with a red border:
    for vp in violin_parts['bodies']:
        vp.set_facecolor('y')
        vp.set_edgecolor('black')
        vp.set_linewidth(1)
        vp.set_alpha(0.3)

    
    labels = age_list
    ax.set_xticks(np.arange(1,len(labels) + 1))
    ax.set_xticklabels(labels)
    
    ax.set_xlim(age_aux[-1]+0.5,age_aux[0]-0.5)
    ax.set_ylim(0.,1.)

    ax.set_xlabel('Age(Gyr)')
    ax.set_ylabel('$a_{AGE}$', fontsize=15)
    

    

################################################################################
################################################################################
################################################################################


def marg_sfh_bar_Z(name,sfh,a_sp, niso,fig):

    Z0, age0, mode = sfh[0], sfh[1], sfh[2]

    Z = np.unique(Z0)
    Z_list = list_ticks(Z)
    Nz=len(Z)
    idx_Z = range(1,Nz+1)

    age = np.unique(age0)
    Nag=len(age)
    age_int = np.append(0.,age)
    
    SFR_mode_marg = marg_Z(mode, Nag, Nz)

    ##
    a_z = []
    for ai in a_sp:
        a_z.append(marg_Z(ai, Nag, Nz))
    a_z = np.array(a_z)
    perc = st_def.a_stat(a_z.T)

    ##
    sfh_mgl=[]
    for i in range(Nz):
        sfh_mgl.append([idx_Z[i], SFR_mode_marg[i], perc[i][0], perc[i][1], perc[i][2]])
    sfh_mgl = np.array(sfh_mgl)

    sfh_mgl = sfh_mgl.T
    
    ###########################################
    ###########################################
    
    ax = fig.add_subplot(133)

    violin_parts = ax.violinplot(a_z, positions=idx_Z, showmedians=True)

    for partname in ('cbars','cmins','cmaxes','cmedians'):
        vp = violin_parts[partname]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    # Make the violin body blue with a red border:
    for vp in violin_parts['bodies']:
        vp.set_facecolor('y')
        vp.set_edgecolor('black')
        vp.set_linewidth(1)
        vp.set_alpha(0.3)
    

    labels = ['0.014', '0.017', '0.020']

    tk=np.arange(1, len(labels) + 1)
    ax.set_xticks(tk)
    ax.set_xticklabels(labels)
    ax.set_xlim(tk[0]-0.5,tk[-1]+0.5)

    ax.set_ylim(0.,1.0)
    
    ax.set_xlabel('Z')
    ax.set_ylabel('$a_Z$', fontsize=15)

