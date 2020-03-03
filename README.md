# traj2oniom
Python script to generate ONIOM (Gaussian) inputs from MD trajectory and topology

The scripts uses [MDAnalysis](https://www.mdanalysis.org/) module to select QM, MM and Point Charges layers. 
It can use the MDAnalysis [syntax](https://www.mdanalysis.org/docs/documentation_pages/selections.html) to set 
layers or use Gromacs index files. Using the former allows dynamic selections along a trajectory.

Example
-------

To generate a QM/MM input with residue DYE in the QM layer and all around 4 Angstrong in the MM layer, plus
the remaining environment within 10 Angs as point charges, one would use:

`traj2oniom.py -f traj.trr -s topol.tpr -selQM "resname DYE" -selMM "byres (around 4 (resname DYE))" -selPC "byres (around 10 (resname DYE))"`

The tested trajectory (traj.trr) and topology (topol.tpr) files correspond to Gromacs files, but in principle, 
any [topology](https://www.mdanalysis.org/docs/documentation_pages/topology/init.html#supported-topology-formats) and [coordinate](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id2) format supported by MDAnalysis should work. The trajectory should be processed before it is used in order
to place the QM molecules in the center and all MM/PC molecules whole and inside the box.
