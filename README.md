# traj2oniom
Python script to generate ONIOM (Gaussian) inputs from MD trajectory and topology

The scripts uses [MDAnalysis](https://www.mdanalysis.org/) module to select QM, MM and Point Charges layers. 
It can use the MDAnalysis [syntax](https://www.mdanalysis.org/docs/documentation_pages/selections.html) to set 
layers or use Gromacs index files. Using the former allows dynamic selections along a trajectory.

Quick guide
-----------

Required input is:

-f (trr,gro,pdb...)  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |Coordinates

-s (tpr...) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |Topology

-QMsel SELECT-COMMAND &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |Sets the QM layer

* Note: alternatively, the selection can be given through a gromacs index file, using the -indQM flag. The first group in the index file is used 

Most relevant optional input is:

-MMsel SELECT-COMMAND &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |Sets the MM layer

-PCsel SELECT-COMMAND &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |Set region to be treated as point charges


See all available options with:

`traj2oniom.py -h`


Example
-------

To generate a QM/MM input with residue DYE in the QM layer and all around 4 Angstrong in the MM layer, plus
the remaining environment within 10 Angs as point charges, one would use:

`traj2oniom.py -f traj.trr -s topol.tpr -selQM "resname DYE" -selMM "byres (around 4 (resname DYE))" -selPC "byres (around 10 (resname DYE))"`


Additional notes
----------------

The tested trajectory (traj.trr) and topology (topol.tpr) files correspond to Gromacs files, but in principle, 
any [topology](https://www.mdanalysis.org/docs/documentation_pages/topology/init.html#supported-topology-formats) and [coordinate](https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#id2) format supported by MDAnalysis should work. 

It is preferable to process the trajectory before it is used in order 
to place the QM molecules in the center and all MM/PC molecules whole 
and inside the box. This can be done internally in the python script
with the -compact option, at the cost of extra ~50% time. Moreover 
the option is not extensively tested.


Limitations
-----------

The script does not (yet) handle bonds between layers properly (i.e. linking atoms).

