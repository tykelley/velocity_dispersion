from argparse import ArgumentParser
from glob import glob
import os.path

import matplotlib as mpl
mpl.use("Agg")  # noqa
from matplotlib import pyplot as plt

from helpers import analyze_halo

parser = ArgumentParser(description="Create a velocity dispersion profile for\
                                     the most massive halo in a simulation.")
parser.add_argument('inbase', help="Path to files")
parser.add_argument('outbase', help="Path to save figures to")
parser.add_argument('snap_num', type=int,
                    help="Snapshot number. [%(default)s]",
                    default=1)
parser.add_argument('--subdirs', dest="subdirs",
                    help="Names of subdirs containing files. ['']",
                    default="")
parser.add_argument('--ptype', dest="ptype", type=int,
                    help="Particle type to use. [%(default)s]",
                    default=1)
parser.add_argument('--rmax', dest="rmax", type=float,
                    help="Distance of furthest shell. [%(default)s kpc]",
                    default=300)
parser.add_argument('--binsize', dest="bin_size", type=float,
                    help="Size of radial bins. [%(default)s kpc]",
                    default=1)
parser.add_argument('--pandas', dest="pandas", action='store_true',
                    help="Use pandas for shell calculations. [%(default)s]",
                    default=False)
args = parser.parse_args()

ptype_map = {0: 'PartType0',
             1: 'PartType1',
             2: 'PartType2',
             4: 'PartType4',
             }
p_type = ptype_map[args.ptype]

fig, ax = plt.subplots(figsize=(12, 12))
mpl.rc('axes', linewidth=3)
plt.yticks(fontsize=10)
plt.xticks(fontsize=10)
plt.tick_params(which='minor', width=2, length=10)
plt.tick_params(which='major', width=2, length=15)

for run_type in args.subdirs.split(','):
    root_path = os.path.join(args.inbase, run_type)
    particle_file = os.path.join(root_path, 'snapshot_' + str(args.snap_num))
    halo_file = os.path.join(root_path, 'halo_' + str(args.snap_num))
    num_of_file = len(glob(particle_file + '*'))
    x, y = analyze_halo(particle_file, halo_file, p_type, num_of_file,
                        rmax=args.rmax, bin_size=args.bin_size,
                        use_pandas=args.pandas)
    ax.semilogx(x, y, linewidth=4.5, color='red', label=run_type)

plt.legend(loc=0)
plt.axis([0, 400, 0, 420])
plt.title('Velocity Dispersion Comparison: M12i dark vs. M12i', fontsize=25)
plt.xlabel(r'$Radius\, [kpc]$', fontsize=30)
plt.ylabel(r'$\sigma_{disp}\, [ km s^{-1}]$', fontsize=30)

fig_name = 'dispersion_' + str(args.snap_num) + '.png'
save_fig_file = os.path.join(args.outbase, fig_name)
fig.savefig(save_fig_file)
