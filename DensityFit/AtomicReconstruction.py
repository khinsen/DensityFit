from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
from MMTK.Minimization import SteepestDescentMinimizer

def atomicReconstruction(ca_pdb, atomic_pdb, out_pdb):
    universe = InfiniteUniverse(Amber94ForceField())
    universe.protein = Protein(atomic_pdb)
    initial = copy(universe.configuration())
    final = Configuration(universe)

    calpha_protein = Protein(ca_pdb, model='calpha')
    for ichain in range(len(calpha_protein)):
        chain = universe.protein[ichain]
        ca_chain = calpha_protein[ichain]
        for iresidue in range(len(ca_chain)):
            ca = chain[iresidue].peptide.C_alpha
            final[ca] = ca_chain[iresidue].position()
        for iresidue in range(len(ca_chain)):
            if iresidue == 0:
                refs = [0, 1, 2]
            elif iresidue == len(ca_chain)-1:
                refs = [iresidue-2, iresidue-1, iresidue]
            else:
                refs = [iresidue-1, iresidue, iresidue+1]
            ref = Collection()
            for i in refs:
                ref.addObject(chain[i].peptide.C_alpha)
            tr, rms = ref.findTransformation(initial, final)
            residue = chain[iresidue]
            residue.applyTransformation(tr)
            ca = residue.peptide.C_alpha
            residue.translateBy(final[ca]-ca.position())
            for a in residue.atomList():
                final[a] = a.position()

    minimizer = SteepestDescentMinimizer(universe)
    minimizer(step = 1000)

    print "Writing configuration to file ", out_pdb
    universe.writeToFile(out_pdb)
