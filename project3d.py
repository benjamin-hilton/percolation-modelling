import scipy as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt
from itertools import repeat

def touchleft(coordn, r):
    """
    Checks if a disk touches the left of the box.
    
    Arguments:
        coordn - an array of the coordinates of the disk
        r - the radius of the disks
    
    Returns:
        True if the disk touches the left of the box
        False if the disk does not touch the left of the box
    """
    if (coordn[0] < r):
        return True
    else:
        return False

def touchright(coordn, r):
    """
    Checks if a disk touches the right of the box.
    
    Arguments:
        coordn - an array of the coordinates of the disk
        r - the radius of the disks
    
    Returns:
        True if the disk touches the right of the box
        False if the disk does not touch the right of the box
    """
    if (coordn[0] >= 1-r):
        return True
    else:
        return False

def distances(N, coord):
    """
    Generates a list of the distances between disks.
    
    Arguments:
        N - the number of disks currently generated
        coord - an array of the coordinates of the disks
        
    Returns:
        dist - a list of arrays of the distances between disks
    """
    dist = []
    for i in xrange(N):
        displ = coord - coord[i]
        dist.append(displ[:,0]**2 + displ[:,1]**2 + displ[:,2]**2)
    return dist

def clusteradd(dist, r, N, clusters, disks):
    """
    Adds a disk N to its clusters or creates a new cluster.

    Arguments:
        dist - a list of arrays of the distances between disks
        r - the radius of the disks
        N - the number of disks currently generated
        clusters - a list of dictionaries describing each cluster containing the keys:
            'indices' - a list of the disk numbers for each disk in the cluster (starting at the first disk being 0)
            'left' - a boolean that is true iff the cluster contains a disk touching the left of the box
            'right' - a boolean that is true iff the cluster contains a disk touching the right of the box
        disks - a list of dictionaries describing each disk containing the keys:
            'coord' - an array of the coordinates of the disk
            'left' - a boolean that is true iff the disk touches the left of the box
            'right' - a boolean that is true iff the disk touches the right of the box

    Returns:
        clusters - an updated list of dictionaries describing each each cluster including the new disk
    """
    index = []
    for i in xrange(N):
        if dist[N-1][i] < 4 * r ** 2:
            k = next((j for j, cluster in enumerate([cluster.get('indices') for cluster in clusters]) if i in cluster), -1)
            if k != -1 and k not in index:
                index.append(k)
    return joincluster(index, clusters, N, disks)
    
def joincluster(index, clusters, N, disks):
    """
    Takes a list of indexes of clusters and joins them together.
    
    Arguments:
        index - a list of indexes of clusters
        clusters - a list of dictionaries describing each cluster containing the keys:
            'indices' - a list of the disk numbers for each disk in the cluster (starting at the first disk being 0)
            'left' - a boolean that is true iff the cluster contains a disk touching the left of the box
            'right' - a boolean that is true iff the cluster contains a disk touching the right of the box
        N - the number of disks currently generated
        disks - a list of dictionaries describing each disk containing the keys:
            'coord' - an array of the coordinates of the disk
            'left' - a boolean that is true iff the disk touches the left of the box
            'right' - a boolean that is true iff the disk touches the right of the box
            
    Returns:
        clusters - an updated list of dictionaries describing each cluster with the specified dictionaries combined
    """
    if index == []:
        clusters.append({'indices' : set([N-1]), 'left' : disks[N-1]['left'], 'right' : disks[N-1]['right']})
    elif index[1:] == []:
        clusters[index[0]] = {'indices' : clusters[index[0]]['indices'].union(set([N-1])), 'left' : clusters[index[0]]['left'] or disks[N-1]['left'], 'right' : clusters[index[0]]['right'] or disks[N-1]['right']}
    else:
        clusters[index[0]] = {'indices' : clusters[index[0]]['indices'].union(set([N-1])), 'left' : clusters[index[0]]['left'] or disks[N-1]['left'], 'right' : clusters[index[0]]['right'] or disks[N-1]['right']}
        for i in index[1:]:
            clusters[index[0]] = {'indices' : clusters[index[0]]['indices'].union(clusters[i]['indices']), 'left' : clusters[index[0]]['left'] or clusters[i]['left'], 'right': clusters[index[0]]['right'] or clusters[i]['right']}
        for i in sorted(index[1:])[::-1]:
            del clusters[i]
    return clusters

def NewDisk(r, N, clusters, coord, disks):
    """
    Creates a new disk and checks for spanning clusters.
    
    Arguments:
        r - the radius of the disks
        N - the number of disks currently generated
        clusters - a list of dictionaries describing each cluster containing the keys:
            'indices' - a list of the disk numbers for each disk in the cluster (starting at the first disk being 0)
            'left' - a boolean that is true iff the cluster contains a disk touching the left of the box
            'right' - a boolean that is true iff the cluster contains a disk touching the right of the box
        coord - a list of arrays of coordinates of each disk
        disks - a list of dictionaries describing each disk containing the keys:
            'coord' - an array of the coordinates of the disk
            'left' - a boolean that is true iff the disk touches the left of the box
            'right' - a boolean that is true iff the disk touches the right of the box

    Returns:
        clusters - an updated list of dictionaries describing each cluster
        coord - an updated list of arrays coordinates of each disk
    """
    coord = sp.append(coord, sp.array(sp.random.uniform(size = (1,3))), axis = 0)
    disks.append({'coord' : coord[-1], 'left' : touchleft(coord[-1], r), 'right' : touchright(coord[-1], r)})
    dist = distances(N, coord)
    clusters = clusteradd(dist, r, N, clusters, disks)
    return clusters,coord


def checkSpanning(clusters):
    """
    Returns the index of a spanning cluster if one exists.
    """
    for i in xrange(len(clusters)):
        if clusters[i]['left'] and clusters[i]['right']:
            return i

def findDensity(r):
    """
    Returns the percolation density for some radius r.
    """
    N=1
    coord = sp.random.uniform(size = (1,3))
    disks = [{'coord' : coord[0], 'left' : touchleft(coord[0], r), 'right' : touchright(coord[0], r)}]
    clusters = [{'indices' : set([0]), 'left' : disks[0]['left'], 'right' : disks[0]['right']}]
    while N > 0:
        clusters,coord = NewDisk(r, N, clusters, coord, disks)
        N += 1
        spanningcluster = checkSpanning(clusters)
        if spanningcluster != None:
            break
    density = (4./3) * (sp.pi * (N-1) * r**3)
    return density

def run (r, iterations):
    """
    Returns the mean percolation density and mean fractional volume covered by disks and their standard deviations.

    Arguments:
        r - the radius of the disks
        iterations - the number of times to find the density for each radius
    """
    print "r=",r
    densities = []
    for i in xrange(iterations):
        densities.append(findDensity(r))
    critvol = 1 - sp.exp(-1 * sp.array(densities))
    meandens = sp.mean(densities)
    stddevdens = sp.std(densities)
    meancrit = sp.mean(critvol)
    stddevcrit = sp.std(critvol)
    return [meandens, stddevdens, meancrit, stddevcrit],[densities,critvol]

def plot_r(x,iterations):
    """
    Plots the mean percolation density and mean fractional volume covered by disks and their standard deviations.
    
    Arguments:
        x - an array of values of r
        iterations - the number of times to find the density for each radius
    """
    y = sp.array([run(i,iterations) for i in x])
    y = sp.vstack(y[:,0][i] for i in xrange(y[:,0].size))
    y_dens = y[:,0]
    y_error_dens = y[:,1]/sp.sqrt(iterations)
    y_crit = y[:,2]
    y_error_crit = y[:,3]/sp.sqrt(iterations)

    plt.figure(figsize=(32,24))

    plt.subplot(2,1,1)
    plt.xlabel('r')
    plt.ylabel(r'$\Delta \bar{\eta_c}$')
    plt.title('Error of mean percolation threshold against radius for ' + str(iterations) + ' iterations')
    #plt.errorbar(x, y_dens, yerr=y_error_dens, color = 'c')
    plt.plot(x, y_error_dens,'b-')
    #plt.axis([0,1,0,5])

    plt.subplot(2,1,2)
    plt.xlabel('r')
    plt.ylabel(r'$\Delta \bar{\phi_c}$')
    plt.title('Error of mean critical volume against radius for ' + str(iterations) + ' iterations')
    #plt.errorbar(x, y_crit, yerr=y_error_crit, color = 'c')
    plt.plot(x, y_error_crit,'b-')
    #plt.axis([0,1,0,1])
    plt.show()

def save(name, filetypes):
    #saves a graph as "name" + ".filetype"
    for i in filetypes:
        plt.savefig(name + "." + i)


r = 0.5 # the radius of the disks for print_results
iterations = 1000 # the number of iterations to perform

#plot_r(sp.arange(0.05,1.00,0.001), iterations)
#save("2Dpercolation1000", ["svg","png","hires.png","pdf"])

x = sp.concatenate([sp.linspace(0.05,0.1,10), sp.linspace(0.1,1.00,1000)])
initial_guess = [3.,0.1,0.2, sp.pi-1.5]
y = [run(i,iterations) for i in x]
y = sp.array(y)
y1 = sp.vstack(y[:,0][i] for i in xrange(y[:,0].size))
y2 = sp.array(list(y[:,1]))

def fit_func(x,k,x0,c, L):
    return ((L) / (1+sp.exp(-k*(x-x0)))) + c

def plot_r_fit(x, y, initial_guess, iterations, fit_func, po, po_cov):
    """
    Fits a function to the data and plots it.
    
    Arguments:
        x - an array of values of r
        y - an array of [mean density, the standard deviation in density, critical volume, the standard deviation in the critical volume] for each value of r
        initial_guess - the initial guess for the fitting function
        iterations - the number of times to find the density for each radius
        fit_func - he function to fit to the data
        po - a list of the parameters to fit to the function
        po_cov - a covariance matrix for the above parameters
    """
    y_dens = y[:,0]
    y_error_dens = y[:,1]/sp.sqrt(iterations)
    y_crit = y[:,2]
    y_error_crit = y[:,3]/sp.sqrt(iterations)

    plt.figure(figsize=(32,24))
    
    plt.subplot(2,1,1)
    plt.errorbar(x, y_dens, yerr=y_error_dens, color = 'c')
    plt.plot(x, y_dens,'b-')
    xx=sp.linspace(0,1,100,0)
    plt.plot(xx,fit_func(xx,po[0],po[1],po[2],po[3]),'r-')
    plt.title('Mean percolation density against radius for ' + str(iterations) + ' iterations')
    plt.xlabel('r')
    plt.ylabel(r'$\bar{\eta}$')

    plt.subplot(2,1,2)
    plt.xlabel('r')
    plt.ylabel(r'$\bar{\phi}$')
    plt.title('Mean fractional volume against radius for ' + str(iterations) + ' iterations')
    plt.errorbar(x, y_crit, yerr=y_error_crit, color = 'c')
    plt.plot(xx,1-sp.exp(-1*(fit_func(xx,po[0],po[1],po[2],po[3]))),'r-')
    plt.plot(x, y_crit,'b-')
    plt.show()
    print "k =", po[0], "+/-", po_cov[0][0]
    print "x0 =", po[1], "+/-", po_cov[1][1]
    print "c =", po[2], "+/-", po_cov[2][2] 
    print "L =", po[3], "+/-", sp.sqrt(po_cov[3][3])
    print "y-intercept = ",fit_func(0,po[0],po[1],po[2],po[3]), "+/-", sp.sqrt(((1/(1-sp.exp(po[0]*po[1])))**2) * po_cov[3][3] + (((po[3]*po[1]*sp.exp(po[0]*po[1]))/(1-sp.exp(po[0]*po[1]))**2)**2) * po_cov[0][0] + ((po[0]*po[3]*sp.exp(po[0]*po[1]))/((1-sp.exp(po[0]*po[1]))**2)**2) * po_cov[1][1] + po_cov[2][2])
    return po,po_cov
    
def plot_scatter_fit(x,y,initial_guess, iterations, fit_func):
    """
    Fits a function to the data as a scatter graph and plots it.
    
    Arguments:
        x - an array of values of r
        initial_guess - the initial guess for the fitting function
        iterations - the number of times to find the density for each radius
        fit_func - he function to fit to the data
    """
    x = [i for item in x for i in repeat(item, iterations)]
    y_dens = [item for sublist in y[:,0] for item in sublist]
    y_crit = [item for sublist in y[:,1] for item in sublist]
    
    po,po_cov=spo.curve_fit(fit_func,x,y_dens,initial_guess)

    plt.figure(figsize=(32,24))

    plt.subplot(2,1,1)
    plt.plot(x, y_dens,'bx')
    xx=sp.linspace(0,1,1000)
    plt.plot(xx,fit_func(xx,po[0],po[1],po[2],po[3]),'r-')
    plt.title('Percolation density against radius for ' + str(iterations) + ' iterations')
    plt.xlabel('r')
    plt.ylabel(r'$\eta$')

    plt.subplot(2,1,2)
    plt.xlabel('r')
    plt.ylabel(r'$\phi$')
    plt.title('Fractional volume against radius for ' + str(iterations) + ' iterations')
    plt.plot(xx,1-sp.exp(-1*(fit_func(xx,po[0],po[1],po[2],po[3]))),'r-')
    plt.plot(x, y_crit,'bx')
    plt.show()
    print "k =", po[0], "+/-", sp.sqrt(po_cov[0][0])
    print "x0 =", po[1], "+/-", sp.sqrt(po_cov[1][1])
    print "c =", po[2], "+/-", sp.sqrt(po_cov[2][2])
    print "L =", po[3], "+/-", sp.sqrt(po_cov[3][3])
    print "y-intercept = ",fit_func(0,po[0],po[1],po[2],po[3]), "+/-", sp.sqrt(((1/(1-sp.exp(po[0]*po[1])))**2) * po_cov[3][3] + (((po[3]*po[1]*sp.exp(po[0]*po[1]))/(1-sp.exp(po[0]*po[1]))**2)**2) * po_cov[0][0] + ((po[0]*po[3]*sp.exp(po[0]*po[1]))/((1-sp.exp(po[0]*po[1]))**2)**2) * po_cov[1][1] + po_cov[2][2])
    return po,po_cov    

po,po_cov=plot_scatter_fit(x,y2,initial_guess,iterations,fit_func)
save("3Dpercolation1000scatterfit", ["svg","png","pdf"])                                                                        
plt.clf()
po,po_cov=plot_r_fit(x, y1, initial_guess,iterations,fit_func, po, po_cov)
save("3Dpercolation1000fitfunc", ["svg","png","pdf"])
