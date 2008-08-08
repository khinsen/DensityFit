# Find a good initial orientation of a protein wrt a density map
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
#

from MMTK import *
from MMTK.Proteins import Protein
import Numeric; N = Numeric

def principalPoints(universe):
    conf = universe.configuration().array
    masses = universe.masses().array
    regions = N.sum(N.array([[1, 2, 4]])*N.greater(conf, 0.), -1)
    principal_points = []
    for i in range(8):
        weight = N.sum(N.equal(regions, i)*masses)
        cm = N.sum((N.equal(regions, i)*masses)[:, N.NewAxis]*conf)/weight
        principal_points.append(Vector(cm))
    return principal_points

def applyPermutation(list, permutation, times=1):
    for n in range(times):
        new = copy(permutation)
        for i in range(len(new)):
            new[i] = list[new[i]]
        list = new
    return list

zrot90 = [1, 3, 0, 2, 5, 7, 4, 6]
yrot180 = [5, 4, 7, 6, 1, 0, 3, 2]
ztoz = [0, 1, 2, 3, 4, 5, 6, 7]
ztox = [4, 0, 6, 2, 5, 1, 7, 3]
ztoy = [4, 5, 0, 1, 6, 7, 2, 3]

def generatePermutations(data):
    permutations = []
    for i in range(4):
        data_zrot = applyPermutation(data, zrot90, i)
        for j in range(2):
            data_yrot = applyPermutation(data_zrot, yrot180, j)
            for p in [ztoz, ztox, ztoy]:
                permutations.append(applyPermutation(data_yrot, p))
    return permutations

def tryOrientations(protein_points, map_points):
    u = InfiniteUniverse()
    for i in range(len(protein_points)):
        u.addObject(Atom('C'))
    c1 = copy(u.configuration())
    c2 = copy(u.configuration())
    for i in range(len(protein_points)):
        c1[i] = protein_points[i]
    fits = []
    for p in generatePermutations(map_points):
        for i in range(len(protein_points)):
            c2[i] = p[i]
        fits.append(u.findTransformation(c1, c2)[0])
    return fits


def initialOrientation(pdb_file, map_file, label):
    
    dmap = load(map_file)
    universe = InfiniteUniverse()
    universe.protein = Protein(pdb_file, model = 'calpha')
    universe.protein.normalizeConfiguration()

    pu = principalPoints(universe)
    pd = dmap.principalPoints()
    fits = tryOrientations(pu, pd)

    initial = copy(universe.configuration())
    for i in range(len(fits)):
        universe.applyTransformation(fits[i])
        fits[i] = (dmap.fit(universe), fits[i])
        universe.setConfiguration(initial)
    fits.sort(lambda a, b: cmp(a[0], b[0]))

    best = fits[0][0]
    i = 0
    while 1:
        fit, tr = fits[i]
        if fit > 1.5*best:
            break
        universe.setConfiguration(initial)
        universe.applyTransformation(tr)
        pdb_out = label + '%d.pdb' % (i+1)
        print "Writing %s, fit error: %10.3g" % (pdb_out, fit)
        universe.writeToFile(pdb_out)
        i = i + 1
