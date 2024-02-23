**A small computational physics project in numpy looking at continuum percolation modelling.**

Percolation theory is the study of long-distance connections in randomly generated systems.

It turns out that in an infinite system, there is some parameter pc such that for all p > pc (where p is the number of
objects in the system), there is an infinite connection spanning the system from one side to the other, and for all
p < pc, there is no such infinite connection. This parameter is called the ‘percolation threshold’ for that system.
In finite systems, this threshold is not so sharply defined: there may be a connection even if p < pc. As such, in order
to find the percolation threshold for a system, one must extrapolate from finite models to models for an infinite system.

By generating large numbers of disks or spheres, estimates for the percolation threshold and critical
fractional area occupied by these disks or spheres are obtained for three different 2 and 3 dimensional
continuous systems: disks of radius r (ηc = 1.126 ± 0.124 and φc = 0.675 ± 0.004), spheres of radius r
(ηc = 0.300±0.096 and φc = 0.256±0.007) and disks of radii r and 2r (ηc = 1.04±0.06 and φc = 0.65±0.02),
where ηc is the critical density and  φc is the critical fractional area (or volume) covered.

The basic method is the same for the modelling of all three systems. Starting with a square (or cube) of side
length 1 simplifies equations 1 to 4 by setting A = 1 and thus n = N. Disks (or spheres) are generated whose
centres are placed at a random point within the square (or cube). The coordinates of each object are recorded,
and each object is given an index (an integer from 1 to N). A list of clusters of disks is kept detailing the
indexes of the disks in that cluster and whether the cluster touches either the right or left side of the box.

For each new object generated, the list of clusters must be updated. First, we check if the new object
touches either the left or right edges of the box, as if a cluster contains an object touching an edge, that cluster
also touches an edge. Then, the distance between the new object and each other object is calculated. These
distances are compared with the sum of the radii of each object and the new object in order to find whether
the new object touches any other.

If no touching objects are found, the new object is placed into a new cluster of its own. If touching objects
are found, the clusters containing these objects are joined together and the new object is added to this enlarged
cluster. Either way, there is one new cluster.

Finally, we check if the new addition to the cluster has caused that cluster to touch an edge, either because
one of the clusters used to form it touched an edge or because the new disk touches an edge.
