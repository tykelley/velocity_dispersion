import h5py
import numpy as np
import pandas as pd


def analyze_halo(particle_file, halo_file, p_type='PartType1', num_of_file=1,
                 rmax=300, bin_size=1, use_pandas=False):
    r, v = get_radius_and_velocity(particle_file, halo_file, num_of_file,
                                    p_type)
    if use_pandas:
        return dispersion_in_shells_pd(r, v, rmax, bin_size)
    else:
        return dispersion_in_shells(r, v, rmax, bin_size)

def dispersion_in_shells(radius, velocities, rmax, bin_size):
    radial_bins = np.arange(0, rmax, bin_size)
    dispersion = np.zeros(radial_bins.size - 1)

    for i in range(1, len(radial_bins)): # Start @ 1 since we're making shells
        mask = (radius >  radial_bins[i-1]) & (radius < radial_bins[i])
        v_avg = velocities[mask].mean(axis=0)
        difference = velocities[mask] - v_avg  
        coord_var = np.var(difference, axis=0)
        dispersion[i-1] = np.sqrt(np.sum(coord_var))
    return radial_bins[1:], dispersion

def dispersion_in_shells_pd(radius, velocities, rmax, bin_size):
    radial_bins = np.arange(0, rmax, bin_size)
    df = pd.DataFrame({"r" : radius,
                       "vx" : velocities[:,0],
                       "vy" : velocities[:,1], 
                       "vz" : velocities[:,2],
                       })
    shells = pd.cut(df.radius, radial_bins)
    groups = df.groupby(shells)
    dispersion = np.sqrt(groups.var().sum(axis=1))
    return radial_bins[1:], dispersion.values 

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

def get_radius_and_velocity(particle_file, halo_file, num_of_file, p_type):
    h = get_data(particle_file, num_of_file,'Header','HubbleParam')
    coords = get_data(particle_file, num_of_file, p_type, 'Coordinates')/h
    velocities  = get_data(particle_file, num_of_file, p_type, 'Velocities')
    masses = get_data(particle_file, num_of_file, p_type, 'Masses')*(10**10)/h

    f = h5py.File(halo_file + '.hdf5', 'r')
    halo_pos = np.array(f['position'])
    halo_mass = np.array(f['mass.vir'])

    idx = np.argmax(halo_mass)
    host_cent = halo_pos[idx]
    radius = np.linalg.norm(coords - host_cent, None, 1)
    return radius, velocities
