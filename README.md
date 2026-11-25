# Bilayer system preparation
A collection of scripts used for preparing MD simulations of lipid bilayers.

Each script uses as a reference and modifies systems prepared with, or with the same folder and name structure, as generated with CHARMM-GUI.

If you use it, make sure to check the code works as intended, and that the results are correct for your systems. 


## File system structure
Scripts were designed to work and tested in a file system like the current repository

```
bilayer-system-preparation/
├── Files
│   ├── lipids
│   │   ├── <MOLECULE>_<FORCEFIELD>.pdb
│   │   ├── <MOLECULE>_<FORCEFIELD>.yaml
│   │   └── <MOLECULE>_<FORCEFIELD>.gro
│   ├── mdp
│   │   └── <FORCEFIELD>
│   │       ├── README
│   │       ├── step6.0_minimization.mdp
│   │       ├── step6.1_equilibration.mdp
│   │       ├── step6.2_equilibration.mdp
│   │       ├── step6.3_equilibration.mdp
│   │       ├── step6.4_equilibration.mdp
│   │       ├── step6.5_equilibration.mdp
│   │       ├── step6.6_equilibration.mdp
│   │       └── step7_production.mdp
│   └── toppar
│       └── <MOLECULE>_<FORCEFIELD>.itp
├── hydroxylate.py
├── README.md
└── Systems
    └── <FORCEFIELD>_<WATERMODEL>
        └── <SYSTEM_NAME>
```

`<MOLECULE>` defines and identifies a certain molecule, that can have the same resname, other, or even a certain number of residues.
`<FORCEFIELD>` can either refer to a defined forcefield or to a combination of them with a specific name, if that were the case (e.g. if both AMBER and GAFF2 parameters were to be used in the same simulation). If a label for a combination of forcefields is defined, files for each molecule labeled with the same label must also be present in the `lipids` and `toppar` folders.
`<SYSTEM_NAME>` is defined by the user; however, it is encouraged to maintain a consistent naming criterium that uniquely identifies the system composition.


## Features

#### hydroxylate.py
Hydroxylate molecules in a lipid membrane system.

```
python hydroxylate.py -i <SYSTEM_NAME_NON_OH> -f <SYSTEM_NAME_OH> \
  -mi <MOLECULE_NON_OH> -mf <MOLECULE_OH> \
  -ri <RESNAME_NON_OH> -rf <RESNAME_OH> -ns <NAME_TO_O> -nr <NAME_REF> \
  -ff <FORCEFIELD> -wm <WATER_MODEL> \
  -l <LIPID_MOLECULE_LIST> [ -ps <PATH_TO_SYSTEMS> -pt <PATH_TO_TOPOL> -pl <PATH_TO_LIPIDS> ]
```

`<SYSTEM_NAME_NON_OH>` is the name of the folder containing the system to hydroxylate. The results would be saved in `<SYSTEM_NAME_OH>`. The molecules `<MOLECULE_NON_OH>` and `<MOLECULE_OH>` differ only in the latter having an extra OH group replacing a certain atom in its structure. `<RESNAME_NON_OH>` indicates the residue in which the hydroxylation is applied, which is then substituted by `<RESNAME_OH>`. `<NAME_TO_O>` is the atom replaced by the OH in the original structure. `<NAME_REF>` is used as a reference to build the OH group; we recommend to indicate the atom to which the replace atom is bonded to. 

`<FORCEFIELD>` and `<WATER_MODEL>` define where the system is stored (see the file system structure). `<LIPID_MOLECULE_LIST>` is a list of `<MOLECULE>` which are present in the final membrane, so make sure **the hydroxylated species is included**.

**Note**: All residues with resname equal to `<RESNAME_NON_OH>` would be hydroxylated. There is no current plan modify the code to allow partial hydroxylation---however, it would not be difficult to do so--. A workaround to do it using the current code could be to manually modify the resname of certain molecules so only those are hydroxylated.
