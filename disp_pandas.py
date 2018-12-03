import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

mpl.use('Agg')


def analyze_halo(particle_file, halo_file, p_type=1):
    h = get_data(particle_file, 1,'Header','HubbleParam')
    coords = get_data(particle_file, 1, p_type, 'Coordinates')/h
    velocities  = get_data(particle_file, 1, p_type, 'Velocities')
    masses = get_data(particle_file, 1, p_type, 'Masses')*(10**10)/h

    f = h5py.File(halo_file + '.hdf5', 'r')
    halo_pos = np.array(f['position'])
    halo_mass = np.array(f['mass.vir'])

    idx = np.argmax(virial)
    host_cent = halo_pos[idx]

    radius = np.linalg.norm(coords - host_cent, None, 1)

    return dispersion_in_shells(velocities, radius)

def dispersion_in_shells(velocities, radius, rmax=300, bin_size=1):
    radial_bins = np.arange(0, rmax, bin_size)
    
    df = pd.DataFrame({"r" : radius,
                       "vx" : velocities[:,0],
                       "vy" : velocities[:,1], 
                       "vz" : velocities[:,2],
                       })
    shells = pd.cut(df.radius, radial_bins)
    groups = df.groupby(shells)

    return np.sqrt(groups.var().sum(axis=1))

def get_data(filename, num_of_file, key1, key2):
    if num_of_file == 1:
        f = h5py.File(filename+'.hdf5', 'r')
        if key1 == 'Header':
            return f[key1].attrs[key2]
        else:
            return f[key1][key2][:]
    else:
        for i in range(0,num_of_file):
            f = h5py.File(filename+'.'+str(i)+'.hdf5', 'r')
            if key1 == 'Header':
                return f[key1].attrs[key2]
            else:
                if ( len(f[key1][key2][:].shape)==1 ):
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.hstack( (result,f[key1][key2][:]) )
                else:
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.vstack( (result,f[key1][key2][:]) )
        return result

if __name__ == "__main__":
    sim = 1 #... Elvis (0) or LATTE (1)
    itype = 1 #...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
    isnap = 600 #...Snapshot
    irun = ['dark', 'm12i'] #... What run?
    my_bins = 150 #...Number of bins in histogram

    part_list = ['PartType0','PartType1','PartType2','PartType4']
    p_type = part_list[itype]

    type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
    type_tag = type_list[itype]

    if sim == 0:
        base = '/data11/home/dmckeown/output/m12i/ELVIS/'
        base1 = '/data11/home/dmckeown/output/m12i/ELVIS/'
    if sim == 1:
        base = '/data18/brunello/dmckeown/output'
        base1 = '/data18/brunello/dmckeown/plots/'


    fig,ax = plt.subplots(figsize=(12,12))
    mpl.rc('axes',linewidth=3)
    plt.yticks(fontsize = 10)
    plt.xticks(fontsize = 10)
    plt.tick_params(which='minor',width=2,length=10)
    plt.tick_params(which='major',width=2,length=15)

    for run_type in irun:
        particle_file = base + '/' + run_type + '/snapshot_'+ str(isnap)
        halo_file = base + '/' + run_type + '/halo_' + str(isnap)

        x,y = analyze_halo(particle_file, halo_file, p_type)
        label = "m12i"
        if "dark" in run_type.lower():
            label += " dark"
        ax.semilogx(x, y, linewidth=4.5, color='red', label=label)

    plt.legend(loc=0)
    plt.axis([0, 400, 0, 420])
    plt.title('Velocity Dispersion Comparison: M12i dark vs. M12i', fontsize = 25)
    plt.xlabel(r'$Radius\, [kpc]$', fontsize = 30)
    plt.ylabel(r'$\sigma_{disp}\, [ km s^{-1}]$', fontsize = 30)

    save_fig_file = base1+str(irun)+'/dm/'+str(isnap)+'.png'
    print "Saving : "
    print str(save_fig_file)
    fig.savefig(save_fig_file)


