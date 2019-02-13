coords = get_data(fname, 1, p_type, 'Coordinates')
velocities  = get_data(fname, 1, p_type, 'Velocities')

radius = ''# Distance from host center


# dispersion1 = []
# rad1 = []
# n = 0
# radii = 0.40
# while n < 29 :

# In general while loops are the slowest of the loop types in Python
# Most of the time you can replace them with a for loop. One of the more
# common times you can't is when you need a dynamic pointer to a list.
# Here, we can use a for loop with a little preparation. 

# Since you are dividing particles into bins, we can define the edges ahead 
# of time and then iterate through the list starting at the second value. 
# Additionally, we know the number of bins, and therefore, how large the
# dispersion list needs to be (number of bin edges - 1).

RMAX = 300
radial_bins = np.arange(0, RMAX, 1)
dispersion = np.zeros(radial_bins.size - 1)

for i in range(1, len(radial_bins)): # Start at i = 1 since we're making shells
    # index = radius > radii
    # master_velocity = velocities[index]

    # index1 = radius < radii + increment
    # master_velocity = velocities[index1]
    
    # I like your idea here of slicing the list but it can be done together.
    # We'll use two boolean masks to do this and combine them with an '&'.
    # Also, I think you have a typo, because you overwrite master_velocity
    # with the use of index1.

    mask = (radius >= radial_bins[i-1]) & (radius < radial_bins[i])

    # v_xav = np.sum((master_velocity[:,0])) / ( len(master_velocity))
    # v_yav = np.sum((master_velocity[:,1])) / ( len(master_velocity))
    # v_zav = np.sum((master_velocity[:,2])) / ( len(master_velocity))
    # v_av = [v_xav,v_yav,v_zav ]

    # difference = np.zeros(len(master_velocity))

    # average_v  = np.repeat(average_velocity,len(master_velocity) )
    # a = master_velocity
    # b = average_v

    # difference = [a_i - b_i for a_i, b_i in zip(a, b)]
    
    # We can take advantage of NumPy's matrix/vector math implementation
    # to do all of this in a few steps. Also, NumPy arrays have a set
    # of convenience methods allowing you to do common statistical 
    # calculations on them (.sum(), .mean(), .std(), .var() and more)

    v_avg = velocities[mask].mean()   
    difference = velocities[mask] - v_avg                                                                                                                           

    # sig = np.zeros(len(difference))
    # for j in range (0, len(difference)):
    #     sigma = np.square(difference[j])
    #     sig[j] = np.sum(sigma)

    # sigg = np.sqrt(sig)
    # disp = np.sum(sigg) / len(sig)
    
    # Again, with the power of NumPy, we can do elementwise calculations
    # Without ever needing to use a for loop. Occasionally, you will need
    # to specify an axis if you don't want to calculate something for the
    # whole array. Here, we need to specify 'axis=1' inside np.sum() to
    # make sure that it sums across (vx + vy + vz) versus summing down
    # each column (axis=0) or the whole thing (default behavior).

    sig = np.sqrt(np.sum(np.square(difference), axis=1))

    # dispersion1.insert(n,disp)
    # rad1.insert(n ,radii)
    # radii = radii + increment
    # n = n + 1
    
    # Finally, we can take our hard work and place it into the 
    # dispersion array we created earlier, remembering that we chose
    # i to start at 1 and not 0.

    dispersion[i-1] = np.mean(sig)
