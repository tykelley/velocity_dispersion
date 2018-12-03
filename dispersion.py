import numpy as np
import h5py
from astropy import constants as const
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pylab import *
import math
from operator import truediv

############################################### FUNCTIONS  ####################\
######################################                                         \
                                          
                                          
                                          
def get_data(filename,num_of_file,key1,key2):
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
        
#... Elvis (0) or LATTE (1)                                                                                                     \
                                                                                                                                 
sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                          \
                                                                                                                                 
itype = 1

#...Snapshot                                                                                                                    \
                                                                                                                                 
isnap = 600
#... What run?                                                                                                                  \
                                                                                                                                 
irun = 'dark'

#...Number of bins in histogram                                                                                                 \
                                                                                                                                 
my_bins = 150
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


fname = base + '/' + irun + '/snapshot_'+ str(isnap)
fname1 = base + '/' + irun + '/halo_' + str(isnap)

###############################################  LOAD DATA  #####################################################               \
                                                                                                                                 


f = h5py.File(fname1 + '.hdf5', 'r')


print f.keys()
print "these are the keys"

h = get_data(fname, 1,'Header','HubbleParam')
print 'h: '+ str(h)
dm_xyz = get_data(fname, 1, p_type, 'Coordinates')/h
velocities  = get_data(fname, 1, p_type, 'Velocities')
print "printing len of velocities"
print len(velocities)

radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])

dm_mass = get_data(fname, 1, p_type, 'Masses')*(10**10)/h
halo_pos = np.array(f['position'])#/h#/(10**3)                                                                                  \
                                                                                                                                 
# this picks most massive halo at index 0                                                                                       \
                                                                                                                                 
virial = np.array(f['mass.vir'])
host_velocity = np.array(f['host.velocity'])

import operator
index, value = max(enumerate(virial), key=operator.itemgetter(1))
virial_radius = np.array(f['radius'])
vir =  virial_radius[index]
host_v = host_velocity[index]
# Note, for some reason host v doesn't give the host velocity at this present time. I find this later manually.
host_cent = halo_pos[index]
# Gives coordinates of the central halo
my_rad = vir
xyz = dm_xyz - host_cent

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])
index1 = radius < 0.3
# Here I'm double checking the halo center given by rockstar to see how close it came!
xyz2 = xyz[index1]
dm_mass2 = dm_mass[index1]
mp = 0
for j in range (0, len(xyz2)):
# sum all the mass times coordinates                                                                                            \
                                                                                                                                 
    mp = dm_mass[j] * xyz2[j,:] + mp
DM_mass_tot1 = len(xyz2) * dm_mass[0]
DM_cm = mp/DM_mass_tot1

print DM_cm

xyz = xyz - DM_cm

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])


# First I define a master index that is all the particles within the virial radius of the given halo
master_index = radius <   my_rad                                                                                                \

xyz = xyz[master_index]
print "printing master xyz"

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])

print radius[0],radius[100], radius [10000]

master_velocities = velocities[master_index]

v_xav = np.sum((master_velocities[:,0])) / ( len(master_velocities))
v_yav = np.sum((master_velocities[:,1])) / ( len(master_velocities))
v_zav = np.sum((master_velocities[:,2])) / ( len(master_velocities))
v_av = [v_xav,v_yav,v_zav ]
velocities = velocities - v_av

# Now repeat to test and make sure its centered

velocities = velocities[master_index]

v_xav = np.sum((master_velocities[:,0])) / ( len(master_velocities))

v_yav = np.sum((master_velocities[:,1])) / ( len(master_velocities))
v_zav = np.sum((master_velocities[:,2])) / ( len(master_velocities))

print "printing average velocities again"
print v_xav, v_yav,v_zav
# These values should all be very nearly 0!!!

HOST_CENTER = np.array(host_cent)
#radius = np.linalg.norm(coords - HOST_CENTER, None, 1)                                                                          

RMAX = 300
radial_bins = np.arange(0, RMAX, 1)
dispersion = np.zeros(radial_bins.size - 1)

for i in range(1, len(radial_bins)): # Start at i = 1 since we're making shells               
    
    mask = (radius >  radial_bins[i-1]) & (radius < radial_bins[i])

    v_avg = velocities[mask].mean(axis=0)
    difference = velocities[mask] - v_avg  
    sig = np.sqrt(np.sum(np.square(difference), axis=1))
    dispersion[i-1] = np.mean(sig)

print radial_bins, dispersion
# we need to delete a radial bin so the len matches
radial_bins = radial_bins[:-1]




#########  NOW TIME FOR FULL BARYONIC RUN  ###################

sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                          \
                                                                                                                             

itype = 1

#...Snapshot                                                                                                                    \
                                                                                                                        

isnap = 600

#... What run?                                                                                                                  \                                                                                                                         

irun = 'm12i'

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


fname = base + '/' + irun + '/snapshot_'+ str(isnap)
fname1 = base + '/' + irun + '/halos_' + str(isnap)
f = h5py.File(fname1 + '.hdf5', 'r')

h = get_data(fname, 4,'Header','HubbleParam')
print 'h: '+ str(h)
dm_xyz = get_data(fname, 4, p_type, 'Coordinates')/h
velocities  = get_data(fname, 4, p_type, 'Velocities')

radius = np.sqrt(dm_xyz[:,0]*dm_xyz[:,0] + dm_xyz[:,1]*dm_xyz[:,1] + dm_xyz[:,2]*dm_xyz[:,2])

dm_mass = get_data(fname, 4, p_type, 'Masses')*(10**10)/h
halo_pos = np.array(f['position'])#/h#/(10**3)                                                                                  \
                                                                                                                                 
virial = np.array(f['mass.vir'])
host_velocity = np.array(f['host.velocity'])

import operator
index, value = max(enumerate(virial), key=operator.itemgetter(1))
print "printing index and value"
print index,value
virial_radius = np.array(f['radius'])
print "printing virial radius"
vir =  virial_radius[index]
print vir
host_v = host_velocity[index]
print "printing all velocities"
print host_velocity
host_cent = halo_pos[index]

my_rad = vir
xyz = dm_xyz - host_cent

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])
index1 = radius < 0.3
xyz2 = xyz[index1]

dm_mass2 = dm_mass[index1]
mp = 0
for j in range (0, len(xyz2)):
# sum all the mass times coordinates                                                                                            \
                                                                                                                                 
    mp = dm_mass[j] * xyz2[j,:] + mp
DM_mass_tot1 = len(xyz2) * dm_mass[0]
DM_cm = mp/DM_mass_tot1
print DM_cm


xyz = xyz - DM_cm

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])

print "printing radius stuff"
print radius[0],radius[100], radius [10000]


master_index = radius <   my_rad
xyz = xyz[master_index]
print "printing master xyz"

radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])
master_velocities = velocities[master_index]
print len(master_velocities)

v_xav = np.sum((master_velocities[:,0])) / ( len(master_velocities))
print "vxav"
print v_xav
v_yav = np.sum((master_velocities[:,1])) / ( len(master_velocities))
v_zav = np.sum((master_velocities[:,2])) / ( len(master_velocities))
print v_xav, v_yav,v_zav
v_av = [v_xav,v_yav,v_zav ]
velocities = velocities - v_av

velocities = velocities[master_index]

v_xav = np.sum((master_velocities[:,0])) / ( len(master_velocities))

v_yav = np.sum((master_velocities[:,1])) / ( len(master_velocities))
v_zav = np.sum((master_velocities[:,2])) / ( len(master_velocities))

print "printing average velocities again"
print v_xav, v_yav,v_zav
halo_mass = np.array(f['mass.vir'])

print 'DM data loaded...'
print "      Host halo center (DM): " +str(host_cent)
print " print length of dm particles "
print len(dm_mass)

HOST_CENTER = np.array(host_cent)
#radius = np.linalg.norm(coords - HOST_CENTER, None, 1)                                                                         \
                                                                                                                                 

RMAX = 300
radial_bins1 = np.arange(0, RMAX, 1)
dispersion1 = np.zeros(radial_bins1.size - 1)

for i in range(1, len(radial_bins1)): # Start at i = 1 since we're making shells  
    mask = (radius >  radial_bins1[i-1]) & (radius < radial_bins1[i])

    v_avg = velocities[mask].mean(axis=0)
    difference = velocities[mask] - v_avg
    sig = np.sqrt(np.sum(np.square(difference), axis=1))
    dispersion1[i-1] = np.mean(sig)
    
print radial_bins1, dispersion1
radial_bins1 = radial_bins1[:-1]

#...Set up figure box:                                                                                                          \
                                                                                                                                 
axes = plt.gca()
fig = plt.figure(figsize = (12,12))
rc('axes',linewidth=3)
plt.yticks(fontsize = 10)
plt.xticks(fontsize = 10)
plt.tick_params(which='minor',width=2,length=10)
plt.tick_params(which='major',width=2,length=15)

ax = fig.add_subplot(1,1,1)
plt.plot(radial_bins, dispersion, linewidth= 4.5, color = 'red',label = "m12i dark")
#plt.plot(rad2, dispersion2, linewidth= 4.5, color = 'red')                                                                      
plt.plot(radial_bins1, dispersion1, linewidth= 4.5, color = 'blue',label = "m12i" )
#plt.plot(rad0, dispersion0, linewidth= 4.5, color = 'blue' )                                                                    

ax.legend(loc=0)

#...Log scale axes                                                                                                              \
                                                                                                                                 
plt.xscale('log')

ax.set_xlim([0,400])
ax.set_ylim([0,420])

plt.title('Velocity Dispersion Comparison: M12i dark vs. M12i', fontsize = 25)
plt.xlabel(r'$Radius\, [kpc]$', fontsize = 30)

plt.ylabel(r'$\sigma_{disp}\, [ km s^{-1}]$', fontsize = 30)
save_fig_file = base1+str(irun)+'/dm/'+str(isnap)+'.png'

print "Saving : "
print str(save_fig_file)

#...Save Figure:                                                                                                                \
                                                                                                                                 
fig.savefig(save_fig_file)


