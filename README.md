# traj2oniom
Python script to generate ONIOM (Gaussian) inputs from MD trajectory and topology

The script uses [MDAnalysis](https://www.mdanalysis.org/) module to select QM, MM and Point Charges layers. 
It can use the MDAnalysis [syntax](https://www.mdanalysis.org/docs/documentation_pages/selections.html) to set 
layers or, instead, Gromacs index files. Using the former allows dynamic selections along a trajectory.

Quick guide
-----------

Required input is:

| Option                | Description      |
|-----------------------|------------------|
| `-f` FILE(trr,gro,pdb...)   |Coordinates       |
| `-s` FILE(tpr...)           |Topology          |
| `-QMsel` SELECT-COMMAND |Sets the QM layer |

* Note: alternatively, the selection can be given through a gromacs index file, using the -indQM flag. The first group in the index file is used 

Most relevant optional input is:

| Option                | Description      |
|-----------------------|------------------|
| `-MMsel` SELECT-COMMAND |Sets the MM layer |
| `-PCsel` SELECT-COMMAND |Set region to be treated as point charges |


See all available options with:

`traj2oniom.py -h`


Example
-------

To generate a QM/MM input with residue DYE in the QM layer and all around 4 Angstrong in the MM layer, plus
the remaining environment within a sphere of 10 Angs as point charges, one would use:

`traj2oniom.py -f traj.trr -s topol.tpr -selQM "resname DYE" -selMM "byres (around 4 (resname DYE))" -selPC "byres (sphlayer 0 10 (resname DYE))"`


Additional notes
----------------

The tested trajectory (traj.trr) and topology (topol.tpr) files correspond to Gromacs files, but in principle, 
any [topology](https://www.mdanalysis.org/docs/documentation_pages/topology/init.html#supported-topology-formats) and [coordinate](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id2) format supported by MDAnalysis should work. 

It is preferable to process the trajectory before it is used in order 
to place the QM molecules in the center and all MM/PC molecules whole 
and inside the box. This can be done internally in the python script
with the `-compact` option, at the cost of extra ~50% time (when loading the system is not the bottleneck step). Moreover 
the option is not extensively tested.

If selections for QM, MM or PC overlap, duplications are removed from
the layer with lower priority. The order or priority is QM>MM>PC

By default, Gaussian inputs are written placing the QM atoms first, followed by MM ones. Optionally, the original order from the trajectory snapshot can be kept with the `-keep` keyword.

Bonds crossing QM/MM boundaries are treated with the linking atom approach. Bonds crossing MM/PC or QM/PC boundaries are not treated in any specific way; the same applies to bonds with atoms in no layer (WARNINGS are issued).

Atom names are kept from the trajectory file. Note that this may lead to errors in Gaussian not understanding some atom names (if they are not standard). A not-well tested option, -fixnames, can be used to change atom names into elements, but this can lead to wrong element assignemts in some case, so use it with caution. A better approach would be regenerating the topology files (tpr) enrusing proper naming of the elements, before running the script.
