#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
hydroxylate.py

Introduce hydroxylation in molecules

@author: Fabián Suárez-Lestón
fabian.suarez.leston@usc.es
CiQUS - Center for Research in Biological Chemistry and Molecular Materials &
Physical Chemistry Department
Universidade de Santiago de Compostela (USC), Spain
'''

import os
import subprocess 

import MDAnalysis
import numpy as np
import argparse

# %%%%%%%%%%     Parse options     %%%%%%%%%%
parser = argparse.ArgumentParser(description='Introduce hydroxylation in molecules')

parser.add_argument('-i', '--name_initial', type=str, help='System name')
parser.add_argument('-f', '--name_final', type=str, help='Target system name')

parser.add_argument('-mi', '--molecule_initial', type=str, help='Name of the original molecule')
parser.add_argument('-mf', '--molecule_final', type=str, help='Name of the target molecule')

parser.add_argument('-ri', '--resname_initial', type=str, help='Resname of the initial molecule')
parser.add_argument('-rf', '--resname_final', type=str, help='Resname of the target molecule')

parser.add_argument('-ns', '--name_subs', type=str, help='Atomname to substitute')
parser.add_argument('-nr', '--name_ref', type=str, help='Atomname to use as reference')

parser.add_argument('-ff', '--forcefield', type=str, help='Forcefield')
parser.add_argument('-wm', '--water_model', type=str, help='Water model')

parser.add_argument('-l', '--lipids', nargs="+", help='List of lipids in the membrane')

parser.add_argument('-ps', '--path_system', type=str, default='./Systems',
                    help='Path to systems. Default: %(default)s')
parser.add_argument('-pt', '--path_topol', type=str, default='./Files/toppar',
                    help='Path to topologies. Default: %(default)s')
parser.add_argument('-pl', '--path_lipids', type=str, default='./Files/lipids',
                    help='Path to lipids. Default: %(default)s')

args = parser.parse_args()


# %%%%%%%%%%     Definitions     %%%%%%%%%%

# names for the original and target folder
SYSTEM = args.name_initial 
NEWSYSTEM = args.name_final 

# names for the original and target molecule
MOLECULE = args.molecule_initial
NEWMOLECULE = args.molecule_final 

# location of the reference atoms to introduce the hydroxylation
REF_RES = args.resname_initial
ATOM_SUBS = args.name_subs
ATOM_REF = args.name_ref

# target residue
NEW_REF_RES = args.resname_final

# force field and water model
FF = args.forcefield 
WM = args.water_model 

# A list of lipid molecules
# Must include the molecule to hydroxylate and any other potentially present in the membrane
LIPIDS = args.lipids 

# Path to relevant files
PATH_TO_SYS = args.path_system 
PATH_TO_ITP = args.path_topol
PATH_TO_LIPIDS = args.path_lipids 


# %%%%%%%%%%     Code     %%%%%%%%%%

if __name__ == "__main__":
    PATH_TO_SYSTEM = f"{PATH_TO_SYS}/{FF}_{WM}"
    print( "hydroxylate.py" )
    print( " Introduce hydroxylation in molecules" )

    print("~ Preparing the system")

    # Target hydroxylated molecule
    a = MDAnalysis.Universe( f"{PATH_TO_LIPIDS}/{NEWMOLECULE}_{FF}.gro" )

    # List of possible lipid residues
    RES_LIPIDS = []
    for lipid in LIPIDS:
        l = MDAnalysis.Universe(fr"{PATH_TO_LIPIDS}/{lipid}_{FF}.gro")
        RES_LIPIDS.extend( list(np.unique(l.residues.resnames)) )

    # Copy the reference system
    if not os.path.isdir( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}" ):
        os.mkdir( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}" )

    subprocess.run( f"cp -r '{PATH_TO_SYSTEM}'/{SYSTEM}/* '{PATH_TO_SYSTEM}/{NEWSYSTEM}'",  shell=True)

    # Change the ITP files
    subprocess.run( f"rm '{PATH_TO_SYSTEM}/{NEWSYSTEM}/toppar/{MOLECULE}.itp'",  shell=True)
    subprocess.run( f"cp '{PATH_TO_ITP}/{NEWMOLECULE}_{FF}.itp' '{PATH_TO_SYSTEM}/{NEWSYSTEM}/toppar/{NEWMOLECULE}.itp'",  shell=True)


    # Reference system
    u = MDAnalysis.Universe( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}/step5_input.gro" )

    new_names = []
    new_coords = []
    new_resindices = []
    new_res = []

    print("~ Starting the hydroxylation")

    for res in u.residues:

        # Modify the molecule to introduce the hydroxylation
        # TODO. Allow partial hydroxylation, both defining certain resids and a fraction of lipids randomly selected
        if res.resname == REF_RES:

            # Use the names of the target molecule
            names = a.atoms.names

            # Location of the atom to substitute
            subs_id = res.atoms.select_atoms( f"resname {REF_RES} and name {ATOM_SUBS}").indices[0] - res.atoms.indices[0]

            # The vector from the reference to the substituted atom 
            vec = res.atoms.select_atoms( f"resname {REF_RES} and name {ATOM_SUBS}").positions - res.atoms.select_atoms( f"resname {REF_RES} and name {ATOM_REF}").positions
            
            # A location for the added H. 
            # The OH vector is aligned with the reference, which is most likely wrong, but this is already fixed in the first minimization step
            # TODO. Introduce the H atom in a better position to prevent divergence in the minimzation for complex, crowded systems
            coord = res.atoms.select_atoms( f"resname {REF_RES} and name {ATOM_SUBS}").positions + vec/np.linalg.norm(vec)

            # Final coordinates for the molecule once hydroxylated
            coords = np.vstack( ( res.atoms.positions[:subs_id+1], coord, res.atoms.positions[subs_id+1:] ) )

            # Store relevant information
            new_names.append( names )
            new_coords.append( coords )
            new_resindices.append( [res.resindex]*len(names) )
            new_res.append( NEW_REF_RES )

        # Save all other molecules with no modification
        else:
            new_names.append( res.atoms.names )
            new_coords.append( res.atoms.positions )
            new_resindices.append( [res.resindex]*len(res.atoms) )
            new_res.append( res.resname )

    print("~ Assembling the final system")

    # Build a new universe with the modified molecules and save it
    v = MDAnalysis.Universe.empty( n_atoms = len( np.concat( new_names ) ), 
                                n_residues = len( new_res ),
                                atom_resindex = sum( new_resindices, [] ),
                                trajectory = True )

    v.add_TopologyAttr( "names", np.concat( new_names ) )
    v.add_TopologyAttr( "resnames", new_res )
    v.add_TopologyAttr( "resids", np.arange(1,len(new_res)+1) )
    v.load_new( np.vstack(new_coords) )

    v.dimensions = u.dimensions

    v.atoms.write( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}/step5_input.gro" )

    print("~ Modifying the topology")

    # Read the topology
    TOPlines = [] 
    with open( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}/topol.top", "r" ) as top:

        NameFound = False

        for line in top.readlines():
            TOPlines.append( line )

            # If the line is commented, continue to the next one
            if line.startswith( ";" ):
                continue

            # Change the hydroxylated molecule resname
            if NameFound:
                try:
                    values = line.split()
                    if REF_RES == values[0]:
                        values[0] = NEW_REF_RES

                    TOPlines[-1] = (8*" ").join( values ) + "\n"
                    NameFound = False
                except:
                    pass

            # Change the hydroxylated molecule ITP 
            if len( line.split(f"{MOLECULE}.itp") ) > 1:
                TOPlines[-1] = f"{NEWMOLECULE}".join( line.split(MOLECULE) )

            # Find the [ molecules ] record
            if line.startswith( "[ molecules ]" ):
                NameFound = True

    # Write the new topology
    with open( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}/topol.top", "w+" ) as top:
        for line in TOPlines:
            top.write( line )

    print("~ Building an index file")

    # Build a new index file
    with MDAnalysis.selections.gromacs.SelectionWriter( f"{PATH_TO_SYSTEM}/{NEWSYSTEM}/index.ndx", "w" ) as ndx:
        ndx.write( v.atoms, name="SYSTEM" )
        
        ndx.write( v.select_atoms( f"resname {' '.join(RES_LIPIDS)}" ), name="MEMB" )

        ndx.write( v.select_atoms( f"not resname {' '.join(RES_LIPIDS)}" ), name="SOLV" )

    print("DONE!")