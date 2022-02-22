#!/usr/bin/env python3

import MDAnalysis
import numpy as np
import sys
import re


def atname2element(name, mass):
    """Turn atomname from MD topology into element name, using the mass to guide to assigment
    
    Input:
     name, string: input name (from MD topology)
     mass, float : atomic mass

    Output:
     name, string: output name (converted if possible, or same as input otherwise)
     
    Note:
     It uses an array with atomic masses per mayor isotopes, instead of average atomic masses. 
     Since it finds the closest mass to the one given, it should no be an issue.
     Requires qcelemental
     """

    import qcelemental as qcel

    # Now uses masses to confirm the name assigment

    # Name from mass
    masses = np.array([qcel.periodictable.to_mass(i) for i in range(118)])
    # Find element with closer mass and compare with input name
    name_found = False
    while not name_found:
        Z = (np.abs(masses - mass)).argmin()
        if abs(masses[Z] - mass) > 10:
            print("Name not found. Returning input")
            return name
        name_qcel = qcel.periodictable.to_E(Z)
        if len(name_qcel) > len(name):
            # This name is not possible: skip picking from masses
            masses[Z] = -1.0
        elif name_qcel.upper() == name[: len(name_qcel)].upper():
            name_found = True
        else:
            # This name is not possible: skip picking from masses
            masses[Z] = -1.0

    return name_qcel


def read_ndx(fname, section=None):
    with open(fname) as f:
        if section:
            for line in f:
                if "[ %s ]" % section in line:
                    break
        else:
            # Reading first section
            line = f.readline()

        # Read all section (till next is encountered
        nums = []
        for line in f:
            if "[" in line:
                break
            nums += [int(c) for c in line.split()]

    return nums


def selections2universelayer(sels, labels):
    """Create a universe merging the atom groups in sels (list) that belong to the same universe. Each one has the resname giving in labels (list)

    Input:
     sels, list of AtomGroups: atom groups (selections) to be merged
     labels, list of str     : resname to be given to the whole selection

    Output:
     univ: Universe with merging the selections. Keeps the original resids
    """

    # Discard None elements
    sels_filt = [i for i in sels if i is not None]
    labels_filt = [j for i, j in zip(sels, labels) if i is not None]

    sel = sels_filt[0]
    for selection in sels_filt[1:]:
        sel += selection
    # Set labels in terms of resids
    resnames = []
    resids = []
    for label, selection in zip(labels_filt, sels_filt):
        resids += list(selection.residues.resids)
        resnames += [label] * selection.residues.n_residues

    # First recreate the resindex map from whole universe to reduced univ
    # NOTE: this will only work for whole residues (if broken it will fail)
    # RESINDICES
    ind = 0
    i_prev = sel.atoms.resindices[0]
    resindices = [ind]
    for i in sel.atoms.resindices[1:]:
        if i != i_prev:
            ind += 1
        resindices.append(ind)
        i_prev = i
    # SEGINDICES (mapping for residues)
    ind = 0
    i_prev = sel.residues.segids[0]
    segindices = [ind]
    for i in sel.residues.resindices[1:]:
        if i != i_prev:
            ind += 1
        segindices.append(ind)
        i_prev = i

    # Create new universe
    univ = MDAnalysis.Universe.empty(
        sel.n_atoms,
        n_residues=sel.n_residues,
        n_segments=sel.n_segments,
        trajectory=True,
        residue_segindex=segindices,
        atom_resindex=resindices,
    )
    n_residues = sel.residues.n_residues

    # Load topology attributes
    univ.add_TopologyAttr("name", sel.atoms.names)
    univ.add_TopologyAttr("type", sel.atoms.types)
    univ.add_TopologyAttr("resname", resnames)
    univ.add_TopologyAttr("resid", resids)
    univ.add_TopologyAttr("segid", sel.segments.segids)
    univ.add_TopologyAttr("segid", sel.segments.segids)
    # And coordinates
    univ.atoms.positions = sel.atoms.positions
    univ.dimensions = sel.dimensions

    return univ


def write_oniom(
    atomsQM,
    atomsMM,
    atomsPC,
    unit=sys.stdout,
    mem="2gb",
    nproc=16,
    chk="snap.chk",
    method="# hf/sto-3g",
    chargemult="0 1",
    title="Gaussian input from snapshot",
    FFfile=None,
    linking_atom="H-H-0.1",
    keep_traj_order=False,
    use_computed_charges=False,
):

    # Preliminary checks
    # Only ONIOM when atomsMM are specified
    if atomsMM:
        hl_flag = "H"
    else:
        hl_flag = ""

    if use_computed_charges:
        # Checking charges
        qreal = np.round(atomsQM.charges.sum() + atomsMM.charges.sum(), 5)
        qmodel = np.round(atomsQM.charges.sum(), 5)
        # Fractional?
        if np.round(np.modf(qreal)[0], 5) != 0.0:
            print("WARNING: fractional charge in real system %s" % qreal)
        if np.round(np.modf(qmodel)[0], 5) != 0.0:
            print("WARNING: fractional charge in model system %s" % qmodel)
        # Consistent with chargemult?
        if chargemult:
            # Take mult from input
            q = float(chargemult.split()[0])
            mreal = chargemult.split()[1]
            if qreal != q:
                pass
                # print('WARNING: sum of charges (real) inconsisten with charge on input: %s vs. %s',qreal,q)
            if len(chargemult.split()) >= 4:
                q = float(chargemult.split()[2])
                mmodel = chargemult.split()[3]
            else:  # same as real
                mmodel = mreal
            if qmodel != q:
                pass
                # print('WARNING: sum of charges (model) inconsisten with charge on input: %s vs. %s',qmodel,q)
        else:
            # Get mult from counting electrons (lowest mult is taken) - TODO (for now assume singlet)
            mreal = str(1)
            mmodel = str(1)

        # Update chargemult string
        qreal = str(int(qreal))
        qmodel = str(int(qmodel))
        chargemult = " ".join([qreal, mreal, qmodel, mmodel])

    # Write header
    header = """%%mem=%s
%%nproc=%s
%%chk=%s

%s

%s

%s""" % (
        mem,
        nproc,
        chk,
        method,
        title,
        chargemult,
    )
    print(header, file=unit)

    # Build QMMM layer
    atomsQMMM = atomsQM + atomsMM
    if keep_traj_order:
        atomsQMMM.atoms.ix_array.sort()

    # Initial preps for QMMM systems
    if atomsMM:
        # We first need a map between traj index and com index
        map_index = {}
        i = 0
        for atom in atomsQMMM.atoms:
            i += 1
            map_index[atom.index] = i
        # Detect qmmm bonds
        qmmm_bonds = {}
        for atom in atomsMM.atoms:
            i1 = map_index[atom.index]
            for bonded_atom in atom.bonded_atoms:
                if bonded_atom in atomsQM.atoms:
                    # QM/MM boundary
                    i2 = map_index[bonded_atom.index]
                    qmmm_bonds[i1] = i2

    # Print layers on Gaussian file
    for atom in atomsQMMM:
        # High layer
        if atom in atomsQM:
            if atomsMM:
                atomlabel = "%s" % (
                    atom.name + "-" + atom.type + "-" + str(round(atom.charge, 5))
                )
            else:
                atomlabel = "%s" % (atom.name)
            print(
                "%-20s %10.5f %10.5f %10.5f %s" % (atomlabel, *atom.position, hl_flag),
                file=unit,
            )
        # Low layer
        else:
            i1 = map_index[atom.index]
            if i1 in qmmm_bonds:
                i2 = qmmm_bonds[i1]
                ll_flag = "L " + linking_atom + " " + str(i2)
            else:
                ll_flag = "L"
            print(
                "%-20s %10.5f %10.5f %10.5f %s"
                % (
                    atom.name + "-" + atom.type + "-" + str(round(atom.charge, 5)),
                    *atom.position,
                    ll_flag,
                ),
                file=unit,
            )

    if atomsMM:
        # Generate and print connectivity
        print("", file=unit)
        for atom in atomsQMMM.atoms:
            i1 = map_index[atom.index]
            cnx_entry = str(i1) + " "
            for bonded_atom in atom.bonded_atoms:
                if bonded_atom.index in map_index:
                    i2 = map_index[bonded_atom.index]
                    cnx_entry += str(i2) + " 1.0 "
                elif bonded_atom in atomsPC.atoms:
                    if atom in atomsQM.atoms:
                        print(
                            "WARNING: atom %s in QM layer bonded to PCsol layer"
                            % atom.name
                        )
                    # else:
                    # print("NOTE: atom %s in MM layer bonded to PCsol layer"%atom.name)
                else:
                    if atom in atomsQM.atoms:
                        print(
                            "WARNING: atom %s in QM layer has an external bond"
                            % atom.name
                        )
                    else:
                        print(
                            "WARNING: atom %s in MM layer has an external bond"
                            % atom.name
                        )
            print("%s" % cnx_entry, file=unit)

        # Include FF
        if FFfile:
            print("", file=unit)
            with open(FFfile) as f:
                for line in f:
                    if len(line.lstrip()) != 0:
                        print(line, file=unit, end="")
    print("", file=unit)

    if atomsPC:
        # Print Point Charges
        for atom in atomsPC:
            print(
                "%10.5f %10.5f %10.5f %10.5f" % (*atom.position, round(atom.charge, 5)),
                file=unit,
            )
        print("", file=unit)

    return None


def write_oniom_gro(atomsQM, atomsMM, atomsPC, grofile, vmd_visualization=True):

    AllLayers = selections2universelayer([atomsQM, atomsMM, atomsPC], 
                                         ["QML", "MML", "PCL"])

    AllLayers.atoms.write(grofile)

    if vmd_visualization:
        with open("viewLayers.vmd", "w") as f:
            print(
                """# Use this as:
#  vmd GROFILE.gro -e viewLayers.vmd
mol delrep 0 top
# QM layer
mol representation CPK 1.200000 0.600000 12.000000 12.000000
mol color Name
mol selection {not resname MML PCL}
mol material Opaque
mol addrep top
# MM layer
mol representation Licorice 0.050000 12.000000 12.000000
mol color Name
mol selection {resname MML}
mol material Opaque
mol addrep top
# Point charges
mol representation Points 3.000000
mol color Name
mol selection {resname PCL}
mol material Opaque
mol addrep top""",
                file=f,
            )

    return None


def compact_atoms(box, atomsQM, atomsMM, atomsPC):

    # Use wrap function
    atomsQM.wrap()
    atomsMM.wrap()
    atomsPC.wrap()

    return None


def exclusive_selection(u, selection, ExclusionGroup):

    selection_exl = "(" + selection + ") and not group Exclusion"

    return u.select_atoms(
        selection_exl, updating=True, Exclusion=ExclusionGroup, periodic=True
    )


if __name__ == "__main__":
    import argparse
    import warnings

    # Input parser. Set flags
    parser = argparse.ArgumentParser(
        description="Generate ONIOM inputs from Gromacs trajectory."
    )
    parser.add_argument("-f", metavar="file.trr", help="Trajectory file", required=True)
    parser.add_argument(
        "-f_fmt", metavar="trr", help="Format of trajectory file", required=False
    )
    parser.add_argument(
        "-s", metavar="file.tpr", help="Binary topoly file", required=True
    )
    parser.add_argument(
        "-s_fmt", metavar="tpr", help="Format of topoly file", required=False
    )
    parser.add_argument(
        "-selQM",
        metavar="select_string",
        help="Selection command for the QM layer",
        default=None,
        required=True,
    )
    parser.add_argument(
        "-selMM",
        metavar="select_string",
        help="Selection command for the MM layer",
        default=None,
    )
    parser.add_argument(
        "-selPC",
        metavar="select_string",
        help="Selection command for the PC layer",
        default=None,
    )
    parser.add_argument(
        "-indQM", metavar="file.ndx", help="Index file for the QM layer", default=None
    )
    parser.add_argument(
        "-indMM", metavar="file.ndx", help="Index file for the MM layer", default=None
    )
    parser.add_argument(
        "-indPC", metavar="file.ndx", help="Index file for the PC layer", default=None
    )
    parser.add_argument(
        "-la",
        metavar="ATOM",
        help="Linking atom label (e.g. H-H-0.1)",
        default="H-H-0.1",
    )
    parser.add_argument(
        "-ob", help="Base name of generated Gaussian input", default="snap"
    )
    parser.add_argument(
        "-osfx", help="Suffix to be appended to Gaussian input name", default=""
    )
    parser.add_argument(
        "-nzero", help="Digits for the label with the step number", type=int, default=3
    )
    parser.add_argument(
        "-b",
        metavar="<time>",
        help="First frame (ps) to read from trajectory",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-e",
        metavar="<time>",
        help="Last frame (ps) to read from trajectory",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-dt",
        metavar="<time>",
        help="Only use frame when t MOD dt = first time (ps)",
        type=float,
        default=-1.0,
    )
    parser.add_argument(
        "-method",
        help='Whole route section (e.g. "#p hf/sto-3g")',
        default="#p hf/sto-3g",
    )
    parser.add_argument(
        "-chargemult",
        help='Charge and multiplicity line (e.g. "0 1 0 1")',
        default=None,
    )
    parser.add_argument(
        "-computeQ",
        action="store_true",
        help="Compute Charges in chargemult from topology",
        default=None,
    )
    parser.add_argument(
        "-FF",
        metavar="file.prm",
        help="FF file into be added to Gaussian input",
        default=None,
    )
    parser.add_argument(
        "-writeGRO",
        action="store_true",
        help="Write gro file to check layers",
        default=False,
    )
    parser.add_argument(
        "-nowrap",
        action="store_true",
        help="Skip wrapping the layers get a compact representation",
        default=False,
    )
    parser.add_argument(
        "-compact",
        action="store_true",
        help="Deprecated option to make system compact (now is the default)",
        default=False,
    )
    parser.add_argument(
        "-keep",
        action="store_true",
        help="Keep atom ordering from trajectory (instead of reordering QM then MM)",
        default=False,
    )
    parser.add_argument(
        "-fixnames",
        action="store_true",
        help="Try to convert atomnames into element names (WARNING: this might work unexpectedly)",
        default=False,
    )
    # Parse input
    args = parser.parse_args()

    if args.compact:
        warnings.warn("Argument -compact is deprecated and is *ignored*.")

    # Get topology and coordinates
    u = MDAnalysis.Universe(
        args.s, args.f, topology_format=args.s_fmt, format=args.f_fmt
    )
    
    # Check if we have box information. If not, wrap is not possible
    if (not args.nowrap) and np.linalg.norm(u.dimensions) == 0.:
        print('Note: no information about the box, setting -nowrap')
        args.nowrap = True

    # Selections
    # -- QM layer --
    if args.selQM:
        try:
            layerQM = u.select_atoms(args.selQM, updating=True)
        except:
            raise BaseException(
                "Error setting QM layer. Maybe due to missplells in selection keyword"
            )
    elif args.indQM:
        atoms_num = read_ndx(args.indQM)
        # Build selection kwd
        selection = "bynum "
        for num in atoms_num:
            selection += str(num) + " "
        # Apply selection (static)
        try:
            layerQM = u.select_atoms(selection, updating=False)
        except:
            raise BaseException("Error setting QM layer. Check index file")
    else:
        raise BaseException("No QM layer was set.")

    # -- MM layer --
    if args.selMM:
        try:
            # Ensure no overlap with QM layer
            # layerMM = exclusive_selection(u,args.selMM,layerQM)
            layerMM = u.select_atoms(args.selMM, updating=True)
        except:
            raise BaseException(
                "Error setting MM layer. Maybe due to missplells in selection keyword"
            )
    elif args.indMM:
        atoms_num = read_ndx(args.indMM)
        # Build selection kwd
        selection = "bynum "
        for num in atoms_num:
            selection += str(num) + " "
        # Apply selection (static)
        try:
            # Ensure no overlap with QM layer
            # layerMM = exclusive_selection(u,selection,layerQM)
            layerMM = u.select_atoms(selection, updating=False)
        except:
            raise BaseException("Error setting MM layer. Check index file")
    else:
        # return an empty AtomGroup (acts as None type within if clauses)
        layerMM = MDAnalysis.AtomGroup([], u)

    # -- Point charges layer --
    if args.selPC:
        try:
            # Ensure no overlap with QM/MM layers
            # layerPC = exclusive_selection(u,args.selPC,layerQM+layerMM)
            layerPC = u.select_atoms(args.selPC, updating=True)
        except:
            raise BaseException(
                "Error setting PC layer. Maybe due to missplells in selection keyword"
            )
    elif args.indPC:
        atoms_num = read_ndx(args.indPC)
        # Build selection kwd
        selection = "bynum "
        for num in atoms_num:
            selection += str(num) + " "
        # Apply selection (static)
        try:
            # Ensure no overlap with QM/MM layers
            # layerPC = exclusive_selection(u,selection,layerQM+layerMM)
            layerPC = u.select_atoms(selection, updating=False)
        except:
            raise BaseException("Error setting PC layer. Check index file")
    else:
        layerPC = MDAnalysis.AtomGroup([], u)

    # Fix names if requested (do on universe only once)
    if args.fixnames:
        for atom in u.atoms:
            atom.name = atname2element(atom.name, atom.mass)

    # Set some defaults to slice the traj
    if args.b == -1.0:
        tini = u.trajectory.time
    else:
        tini = args.b
    if args.e == -1.0:
        tfin = u.trajectory.totaltime
    else:
        tfin = args.e
    if args.dt == -1.0:
        dt = u.trajectory.dt
    else:
        dt = args.dt

    istp = 0
    for conf in u.trajectory:
        # Determine whether using the frame or not
        if u.trajectory.n_frames == 1:
            # Always use the frame when there is only one
            pass
        elif conf.time > tfin + 0.0001:
            sys.exit()
        elif conf.time % dt > 0.0001:
            continue
        elif conf.time < tini:
            continue

        # Default is trying to generate a compact representation by wrapping all layers
        if not args.nowrap:
            # Get box dimensions (assume rect box)
            box = u.dimensions[:3]
            compact_atoms(box, layerQM, layerMM, layerPC)

        # Write Gaussian input
        if args.nzero > 0:
            fmt = "%%s%%0%ig%%s.%%s" % args.nzero
            fname = fmt % (args.ob, istp, args.osfx, "com")
            chkname = fmt % (args.ob, istp, args.osfx, "chk")
        else:
            fmt = "%s%s.%s"
            fname = fmt % (args.ob, args.osfx, "com")
            chkname = fmt % (args.ob, args.osfx, "chk")
        f = open(fname, "w")
        # Manage defaults to set charge/mult
        if args.computeQ is None:
            if args.chargemult is not None:
                args.computeQ = False
            else:
                args.computeQ = True
        write_oniom(
            layerQM,
            layerMM - layerQM,
            layerPC - layerMM - layerQM,
            unit=f,
            title="Trajectory step %s (time=%s ps)" % (conf.frame, round(conf.time, 5)),
            chk=chkname,
            FFfile=args.FF,
            method=args.method,
            chargemult=args.chargemult,
            linking_atom=args.la,
            keep_traj_order=args.keep,
            use_computed_charges=args.computeQ,
        )
        f.close()
        files_gen = fname
        if args.writeGRO:
            if args.nzero > 0:
                grofile = fmt % (args.ob, istp, args.osfx, "gro")
            else:
                grofile = fmt % (args.ob, args.osfx, "gro")
            write_oniom_gro(
                layerQM,
                layerMM - layerQM,
                layerPC - layerMM - layerQM,
                grofile,
                vmd_visualization=True,
            )
            files_gen += " " + grofile + " viewLayers.vmd"
        # Indicate all files generated
        print("%s generated" % files_gen)
        # Prepare next iteration
        istp += 1
