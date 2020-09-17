import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

def getForceByClass(system, klass):
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if isinstance(f, klass):
            return f

pdb = PDBFile('input.pdb')
forcefield = ForceField('tip3p.xml','opls-aam.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
for i in range(system.getNumForces()):
    force = system.getForce(i)
    force.setForceGroup(i)


print('How many harmonic bonds?  ', getForceByClass(system, HarmonicBondForce).getNumBonds())
print('How many harmonic angles? ', getForceByClass(system, HarmonicAngleForce).getNumAngles())
print('How many torsions?        ', getForceByClass(system, PeriodicTorsionForce).getNumTorsions())


integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
context = Context(system, integrator)
context.setPositions(pdb.positions)

for i in range(system.getNumForces()):
    mean_magnitude = np.mean(np.linalg.norm(context.getState(getForces=True, groups=1<<i).getForces(asNumpy=True)._value, axis=1))
    print('Force %d (%20s): mean magniude: %.3g' % (i, system.getForce(i).__class__.__name__, mean_magnitude))