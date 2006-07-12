# Fit protein structure flexibly into density map.
#
# Written by Konrad Hinsen <hinsen@llb.saclay.cea.fr>
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.NormalModes import NormalModes, SparseMatrixSubspaceNormalModes
from MMTK.FourierBasis import FourierBasis, countBasisVectors
from MMTK.Random import randomParticleVector
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
import Numeric; N = Numeric


# Percentage of random force contribution
random = 20.

# Calculate a given number of normal modes
def calc_modes(universe, nmodes):
    natoms = universe.numberOfAtoms()
    nmodes = min(nmodes, 3*natoms)
    if nmodes > natoms:
        nmodes = 3*natoms
    if nmodes == 3*natoms:
        modes = NormalModes(universe, None)
        modes.cutoff = None
    else:
        p1, p2 = universe.boundingBox()
        cutoff_max = (p2-p1).length()
        cutoff = 0.5*cutoff_max
        nmodes_opt = nmodes
        nmodes = countBasisVectors(universe, cutoff)
        while nmodes > nmodes_opt:
            cutoff = cutoff + 0.1
            if cutoff > cutoff_max:
                cutoff = cutoff_max
                break
            nmodes = countBasisVectors(universe, cutoff)
        while nmodes < nmodes_opt:
            cutoff = cutoff - 0.1
            if cutoff < 0.1:
                cutoff = 0.1
                break
            nmodes = countBasisVectors(universe, cutoff)
        basis = FourierBasis(universe, cutoff)
        basis.may_modify = 1
        modes = SparseMatrixSubspaceNormalModes(universe, basis)
        modes.cutoff = cutoff
    return modes

# Remember the initial nearest-neighbour C-alpha distances
def set_distances(protein):
    for chain in protein:
        for i in range(1, len(chain)):
            d = (chain[i].position()-chain[i-1].position()).length()
            chain[i].neighbour_distance = d

# Re-set the initial nearest-neighbour C-alpha distances
def fix_structure(protein):
    for chain in protein:
        for i in range(1, len(chain)):
            d = chain[i].neighbour_distance
            r1 = chain[i-1].position()
            r2 = chain[i].position()
            chain[i].peptide.C_alpha.setPosition(r1+d*(r2-r1).normal())

# Weighting for the current mode window
def weight_function(mode, center):
    return 1./(1.-(mode.frequency-center)**2)

# The main fit function
def fit(pdb_file, map_file, energy_cutoff, label, trajectory_freq):
    dmap = load(map_file)
    universe = InfiniteUniverse(CalphaForceField(2.5))
    universe.protein = Protein(pdb_file, model='calpha')
    set_distances(universe.protein)

    modes = calc_modes(universe, 3*universe.numberOfAtoms())
    for nmodes in range(6, len(modes)):
        if modes[nmodes].frequency > energy_cutoff*modes[6].frequency:
            break

    if trajectory_freq > 0:
        trajectory_filename = 'fit_%s_m%d.nc' % (label, nmodes)
        print "Fit trajectory:", trajectory_filename
        trajectory = Trajectory(universe, trajectory_filename, 'w',
                                'Density fit %s using %d modes'
                                % (label, nmodes))
        actions = [TrajectoryOutput(trajectory,
                                    ["configuration", "fit_error"],
                                    0, None, 1)]
        snapshot = SnapshotGenerator(universe, actions=actions)
        snapshot.available_data.append('fit_error')

    amplitude = 0.1
    i = 0
    fit = []
    windex = 0
    protocol_filename = 'fit_%s_m%d_fiterror.txt' % (label, nmodes)
    print "Fit error protocol file:", protocol_filename
    protocol = file(protocol_filename, 'w')

    nmc = 2*nmodes
    modes = calc_modes(universe, nmc)
    windows = N.arange(10.)*modes[nmodes].frequency/10.

    while 1:
        f, gradient = dmap.fitWithGradient(universe)
        fit.append(f)
        print >> protocol, i, fit[-1]
        protocol.flush()
        if len(fit) > 1 and fit[-1] > fit[-2]:
            amplitude /= 2.
        random_width = gradient.norm()*(random/100.) \
                       / N.sqrt(universe.numberOfAtoms())
        random_vector = randomParticleVector(universe, random_width)
        gradient = gradient + random_vector
        displacement = ParticleVector(universe)
        for n in range(nmodes):
            m = modes[n]
            if n < 6:
                weight = 5.*weight_function(modes[6], windows[windex])
            else:
                weight = weight_function(m, windows[windex])
            displacement = displacement + \
                           weight*m.massWeightedDotProduct(gradient)*m
        displacement = displacement*amplitude/displacement.norm()
        universe.setConfiguration(universe.configuration() + displacement)
        fix_structure(universe.protein)
        if trajectory_freq > 0 and i % trajectory_freq == 0:
            snapshot(data = {'fit_error': fit[-1]})
        i = i + 1
        if i % 20 == 0:
            modes = calc_modes(universe, nmc)
        if i % 10 == 0 and i > 50:
            convergence = (fit[-1]-fit[-11])/(fit[30]-fit[20])
            if convergence < 0.05:
                break
            if convergence < 0.2:
                windex = min(windex+1, len(windows)-1)
            if convergence < 0.1:
                amplitude = 0.5*amplitude
    protocol.close()

    if trajectory_freq > 0:
        trajectory.close()

    pdb_out = 'fit_%s_m%d.pdb' % (label, nmodes)
    print "Writing final structure to", pdb_out
    universe.writeToFile(pdb_out)
