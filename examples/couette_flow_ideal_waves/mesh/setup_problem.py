#!/usr/bin/env python3
"""Setup the channel flow problem
"""

# ========================================================================
#
# Imports
#
# ========================================================================
import os
import yaml
import numpy as np
import subprocess as sp

# ========================================================================
#
# Setup
#
# ========================================================================
delta = 1
L = [6.,5.,1.]
N = [144,96,96]

polynomial_order = 1.
Re_tau = 130
Re=8000. 
viscosity = 1./Re
U0=1.

dbname = 'couette_ideal_waves'
msh_dbname = dbname + ".exo"
ic_dbname = dbname + "_ic.exo"

mpibin = '/Users/gdeskos/spack/opt/spack/darwin-mojave-x86_64/clang-10.0.1-apple/openmpi-3.1.4-rq4lu5lmrlvrisakvdqxg7ar73i36n7t/bin/mpiexec'
abl_mesh = '/Users/gdeskos/nalu-wind/wind-utils/build/src/mesh/abl_mesh'
nalu_preprocess = '/Users/gdeskos/nalu-wind/wind-utils/build/src/preprocessing/nalu_preprocess'

# ========================================================================
#
# Generate new mesh and IC files
#
# ========================================================================

# Load the skeleton data
msh_iname = "mesh.yaml"
msh_inp = yaml.load(open(msh_iname, 'r'))
msh_data = msh_inp['nalu_abl_mesh']

ic_iname = "ic.yaml"
ic_inp = yaml.load(open(ic_iname, 'r'))
ic_data = ic_inp['nalu_preprocess']

# New yaml mesh file
msh_oname = "msh_tmp.yaml"
msh_data['output_db'] = msh_dbname
vertices = msh_data['vertices']
vertices[1][0] = L[0]
vertices[1][1] = L[1]
vertices[1][2] = L[2]
msh_data['mesh_dimensions'] = N
z_spacing = msh_data['z_spacing']
z_spacing['stretching_factor'] = 1.2**(polynomial_order / 2.0)

with open(msh_oname, 'w') as of:
    yaml.dump(msh_inp, of, default_flow_style=False)

# New yaml IC file
ic_oname = "ic_tmp.yaml"
ic_data['input_db'] = msh_dbname
ic_data['output_db'] = ic_dbname
task = ic_data['init_couette_fields']
velocity = task['velocity']
velocity['Re_tau'] = Re_tau
velocity['U0'] = U0 
velocity['viscosity'] = viscosity

with open(ic_oname, 'w') as of:
    yaml.dump(ic_inp, of, default_flow_style=False)

# ========================================================================
#
# Run the utilities
#
# ========================================================================

proc = sp.Popen(mpibin + ' -np 1 ' + abl_mesh + ' -i ' + msh_oname,
                shell=True,
                stderr=sp.PIPE)
err = proc.communicate()
errcode = proc.returncode

proc = sp.Popen(mpibin + ' -np 1 ' + nalu_preprocess + ' -i ' + ic_oname,
                shell=True,
                stderr=sp.PIPE)
err = proc.communicate()
errcode = proc.returncode

# ========================================================================
#
# Clean up
#
# ========================================================================
os.remove(msh_oname)
os.remove(ic_oname)
