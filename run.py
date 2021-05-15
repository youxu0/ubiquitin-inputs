#!/usr/bin/env python3
# usage: ./run.py start stop pressure temperature name forcefield

from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout, exit, stderr
from subprocess import call

start=int(sys.argv[1]) # the number to run
stop=int(sys.argv[2]) # the number to stop
pressure=float(sys.argv[3]) # pressure in bar
temp=float(sys.argv[4]) # temperature in Kelvin
jobname=sys.argv[5]+'-'+sys.argv[3]+'_'+sys.argv[4]+'-'+sys.argv[6]

#temp=300.0
#pressure=1.0
#jobname='ubi.c36m'

wkPath='../'+sys.argv[3]+'_'+sys.argv[4]+'/'
trjPath='../'+sys.argv[3]+'_'+sys.argv[4]+'/traj/' # directory to write dcd

# Read the PSF
psf = NewCharmmPsfFile('ubi_cubi.psf')

boxsize=5.85 # boxsize in nanometer
psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)

# Get the coordinates from the PDB
crd = app.CharmmCrdFile('ubi_cubi.crd') # in case a crd file instead of pdb is used

# Load the parameter set.
lib = '/home/xuyou/toppar/' # FULL path of force field topology and parameter files
params = app.CharmmParameterSet(lib+'top_all36_prot.rtf', lib+'par_all36m_prot.prm', lib+'toppar_water_ions.str')

platform = mm.Platform.getPlatformByName('CUDA') # make sure the GPU is used. Others: Reference, CPU & OpenCL

system = psf.createSystem(params, nonbondedMethod=app.LJPME, nonbondedCutoff=0.9*nanometer, ewaldErrorTolerance = 0.0001, constraints=app.HBonds) 
# The error tolerance is roughly equal to the fractional error in the forces due to truncating the Ewald summation. Default 0.0005. 
# The default vdwCutoff is as nonbondedCutof, so fswitch or fshift?

system.addForce(mm.AndersenThermostat(temp*kelvin, 1/picosecond)) # thermostat: temperature and collision frequency
system.addForce(mm.MonteCarloBarostat(pressure*bar, temp*kelvin)) # The barostat does not itself do anything to regulate the temperature
integrator = mm.VerletIntegrator(0.002*picoseconds)        # 2 fs time step

simulation = app.Simulation(psf.topology, system, integrator, platform)

# simulation.context.setPositions(pdb.getPositions()) # speicify the initial positions
simulation.context.setPositions(crd.positions) # in case the crd file is used 

simulation.minimizeEnergy(maxIterations = 100) # default tolerance=10*kJ/mol, maxIterations: until converged

if start > 0: # restart from the positions and velocities
    restart=str(start-1)
    with open(trjPath+jobname+'.'+restart+'.rst', 'r') as f:  # the .rst file should only be used if .chk restart file cannot be loaded
        simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))

nsavcrd=50000 # save coordinates every 100 ps
nstep=10000000 # write dcd files every 20 ns
nprint=50000 # report every 1 ns

for ii in range(start,stop+1):
    dcd=app.DCDReporter(trjPath+jobname+'.'+str(ii)+'.dcd', nsavcrd)
    firstdcdstep = ii*nstep + nsavcrd
    while (firstdcdstep > 2000000000): # reset frame number to avoid charmm overfloat
        firstdcdstep -= 2000000000
    dcd._dcd = app.DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0
    simulation.reporters.append(dcd) # "reporters" aka. write out
    simulation.reporters.append(app.StateDataReporter(wkPath+jobname+'.'+str(ii)+'.out', nprint, step=True, time=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))
    """StateDataReporter options: step, time, progress, remainingTime, potentialEnergy, kineticEnergy, totalEnergy, temperature, volume, density, speed, totalSteps. Default separator=','"""
    simulation.step(nstep) # run the simulation
    simulation.reporters.pop() # temporally keep to continue the loop
    simulation.reporters.pop()
    dcd._out.close() # close the dcd file to make sure the buffers are written.
 
    # write restart file
    state = simulation.context.getState( getPositions=True, getVelocities=True )
    with open(trjPath+jobname+'.'+str(ii)+'.rst', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))

    # move the dcd back to home
    #cmdstr='mv '+trjPath+jobname+'.'+str(ii)+'.dcd .'
    #call(cmdstr, shell=True)


