import matplotlib as mpl
from matplotlib import pyplot as plt

from helpers import analyze_halo

mpl.use('Agg')

if __name__ == "__main__":
    sim = 1 #... Elvis (0) or LATTE (1)
    itype = 1 #...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
    isnap = 600 #...Snapshot
    irun = ['dark', 'm12i'] #... What run?

    part_list = ['PartType0','PartType1','PartType2','PartType4']
    p_type = part_list[itype]

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
