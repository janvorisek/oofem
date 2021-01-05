nonstat_nonlin_01.out.tm
Additive manufacturing - heat transfer
TransientTransport nsteps 2000 deltat 0.05 alpha 0.5 lumped exportfields 1 5 nmodules 1 smtype 5 lstype 1 lsprecond 1
vtkxml tstep_step 100 domain_all primvars 1 6
domain heattransfer
OutputManager
ndofman 0 nelem 0 ncrosssect 1 nmat 1 nbc 0 nic 0 nltf 1 nset 0
# Do not set any set for transport section!!!! done by voxel activator
SimpleTransportCS 1 thickness 1.0 mat 1
IsoHeat 1 d 1020. k 0.21 c 1600.0
ConstantFunction 1 f(t) 1.0