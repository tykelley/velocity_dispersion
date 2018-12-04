from argparse import ArgumentParser

import matplotlib as mpl
from matplotlib import pyplot as plt

from helpers import analyze_halo

parser = ArgumentParser()
parser.add_argument('inbase', help="Path to files.")
parser.add_argument('outbase',
                    help="Path to save figures to.")
parser.add_argument('snap_num', type=int,
                    help="Snapshot number. (%(default)s)",
                    default=1)
parser.add_argument('--rmax', dest="rmax", type=float,
                    help="Distance of furthest shell. (%(default)s kpc)",
                    default=300)
parser.add_argument('--binsize', dest="bin_size", type=float,
                    help="Size of radial bins. (%(default)s kpc)",
                    default=1)
parser.add_argument('--pandas', dest="pandas", action='store_true', 
                    help="Use pandas for shell calculations. (%(default)s)",
                    default=False)

args = parser.parse_args()

# Options 
sim = 1 #... Elvis (0) or LATTE (1)
itype = 1 #...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
isnap = 600 #...Snapshot


irun = ['dark', 'm12i'] #... What run?
part_list = ['PartType0','PartType1','PartType2','PartType4']
p_type = part_list[itype]

if sim == 0:
    base = '/data11/home/dmckeown/output/m12i/ELVIS'
    base1 = '/data11/home/dmckeown/output/m12i/ELVIS'
if sim == 1:
    base = '/data18/brunello/dmckeown/output'
    base1 = '/data18/brunello/dmckeown/plots'


fig,ax = plt.subplots(figsize=(12,12))
mpl.rc('axes',linewidth=3)
plt.yticks(fontsize = 10)
plt.xticks(fontsize = 10)
plt.tick_params(which='minor',width=2,length=10)
plt.tick_params(which='major',width=2,length=15)

for run_type in irun:
    particle_file = args.inbase + '/' + run_type + '/snapshot_'+ str(args.snap_num)
    halo_file = arg.inbase + '/' + run_type + '/halo_' + str(args.snap_num)

    label = "m12i"
    if "dark" in run_type.lower():
        label += " dark"
        num_of_file = 1
    else:
        num_of_file = 4

    x,y = analyze_halo(particle_file, halo_file, p_type, num_of_file, 
                        rmax=args.rmax, bin_size=args.binsize,
                        use_pandas=args.pandas)
    ax.semilogx(x, y, linewidth=4.5, color='red', label=label)

plt.legend(loc=0)
plt.axis([0, 400, 0, 420])
plt.title('Velocity Dispersion Comparison: M12i dark vs. M12i', fontsize = 25)
plt.xlabel(r'$Radius\, [kpc]$', fontsize = 30)
plt.ylabel(r'$\sigma_{disp}\, [ km s^{-1}]$', fontsize = 30)

save_fig_file = args.outbase+'dispersion_'+str(args.snap_num)+'.png'
fig.savefig(save_fig_file)
