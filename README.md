# MD_dihedral_tracker
Plot any dihedral angle along a .xyz trajectory file

Usage:

./MD_dihedral_tracker i j k l molecule.xyz

where:
- molecule.xyz is the trajectory in standard xyz format
- i, j, k and l are indexes of the atoms defining the torsion (indices are 1-based)

The easiest way to install MD_dihedral_tracker is with gfortran:

```
   git clone https://github.com/thomasjamespope/MD_dihedral_tracker.git
   cd MD_dihedral_tracker
   gfortran -ffree-form -O3 -o MD_dihedral_tracker.x MD_dihedral_tracker.f
```
