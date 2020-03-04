#!/usr/bin/env python3

import MDAnalysis
import numpy as np
import sys
import re


def atname2element(name,resname):
    # This wont work for Nb...

    import string

    # Manage C-labels and CA, CL, CS
    if   resname.upper() == 'CA' and name.upper() == 'CA':
        name = 'Ca'
    elif resname.upper() == 'CL' and name.upper() == 'CL':
        name = 'Cl'
    elif resname.upper() == 'CS' and name.upper() == 'CS':
        name = 'Cs'
    elif name.upper()[0] == 'C':
        name = 'C'
    # Manage N-labels
    if   resname.upper() == 'NA' and name.upper() == 'NA':
        name = 'Na'
    elif name.upper()[0] == 'N':
        name = 'N'
    # Manage solvent
    if   resname.upper() == 'SOL' and name.upper() == 'HO':
        name = 'H'
    elif resname.upper() == 'HOH' and name.upper() == 'HO':
        name = 'H'
    elif resname.upper() == 'WAT' and name.upper() == 'HO':
        name = 'H'
    elif resname.upper() == 'MET' and name.upper() == 'HO':
        name = 'H'
    elif resname.upper() == 'ETH' and name.upper() == 'HO':
        name = 'H'
    if   resname.upper() == 'SOL' and name.upper() == 'OH':
        name = 'O'
    elif resname.upper() == 'HOH' and name.upper() == 'OH':
        name = 'O'
    elif resname.upper() == 'WAT' and name.upper() == 'OH':
        name = 'O'
    elif resname.upper() == 'MET' and name.upper() == 'OH':
        name = 'O'
    elif resname.upper() == 'ETH' and name.upper() == 'OH':
        name = 'O'
    # Unambiguous labels
    if name.upper() == 'HC':
        name = 'H'
    # Erase numbers and extrange chars
    # (from: https://bytes.com/topic/python/answers/850562-finding-all-numbers-string-replacing)
    name = re.sub('[%s]' % string.digits, '', name)

    return name

def read_ndx(fname,section=None):
    with open(fname) as f:
        if section:
            for line in f:
                if '[ %s ]'%section in line:
                    break
        else:
            # Reading first section
            line = f.readline()

        # Read all section (till next is encountered
        nums = []
        for line in f:
            if '[' in line:
                break
            nums += [int(c) for c in line.split()]

    return nums


def write_oniom(atomsQM,atomsMM,atomsPC,
                mem='2gb',nproc=16,chk='snap.chk',
                method="# hf/sto-3g",chargemult="0 1",unit=sys.stdout,
                title='Gaussian input from snapshot',FFfile=None):

    # Preliminary checks
    # Only ONIOM when atomsMM are specified
    if atomsMM:
        hl_flag = 'H'
    else:
        hl_flag = ''


    # Write header
    header="""%%mem=%s
%%nproc=%s
%%chk=%s

%s

%s

%s"""%(mem,nproc,chk,method,title,chargemult)
    print(header,file=unit)
    
    # Print High layer 
    for atom in atomsQM:
         atname=atname2element(atom.name,atom.resname)
         if atomsMM:
             atomlabel='%s'%(atname+'-'+atom.type+'-'+str(round(atom.charge,5)))
         else:
             atomlabel='%s'%(atname)
         print('%-20s %10.5f %10.5f %10.5f %s'%(atomlabel,*atom.position,hl_flag),file=unit) 
         
    if atomsMM:
        # Print Low layer
        for atom in atomsMM: 
             atname=atname2element(atom.name,atom.resname)
             print('%-20s %10.5f %10.5f %10.5f %s'%(atname+'-'+atom.type+'-'+str(round(atom.charge,5)),*atom.position,'L'),file=unit) 

        # Print connectivity
        # We first need a map between traj index and com index
        map_index = {}
        i = 0
        for atom in atomsQM.atoms:
            i += 1
            map_index[atom.index] = i
        for atom in atomsMM.atoms:
            i += 1
            map_index[atom.index] = i
        # And now generate connectivity
        print("",file=unit)
        for atom in atomsQM.atoms:
            i1 = map_index[atom.index]
            cnx_entry = str(i1)+' '
            for bonded_atom in atom.bonded_atoms: 
                i2 = map_index[bonded_atom.index] 
                cnx_entry += str(i2)+' 1.0 '
            print('%s'%cnx_entry,file=unit) 
        for atom in atomsMM.atoms: 
            i1 = map_index[atom.index]
            cnx_entry = str(i1)+' '
            for bonded_atom in atom.bonded_atoms: 
                i2 = map_index[bonded_atom.index] 
                cnx_entry += str(i2)+' 1.0 '
            print('%s'%cnx_entry,file=unit) 
            
        # Include FF
        if FFfile:
            print("",file=unit)
            with open(FFfile) as f:
                for line in f:
                    if len(line.lstrip()) != 0:
                        print(line,file=unit,end='')
    print("",file=unit)

    if atomsPC:
        # Print Point Charges
        for atom in atomsPC: 
             print('%10.5f %10.5f %10.5f %10.5f'%(*atom.position,round(atom.charge,5)),file=unit)
        print("",file=unit)
        
    return None
        
        
def write_oniom_gro(atomsQM,atomsMM,atomsPC,grofile,vmd_visualization=True):
    
    AllLayers = atomsQM.copy()
    # Setting this resnames leads to an error due to
    # a strange conflict with atomsPC (this is not
    # observed e.g. when atomsPC has updating=False)
    # We comment this for the moment, but it'd be better
    # to understand why the error occurs
    #AllLayers.residues.resnames = 'QML'
    
    if atomsMM:
        Layer = atomsMM.copy()
        Layer.residues.resnames = 'MML'
        AllLayers += Layer.copy()
        
    if atomsPC:
        Layer = atomsPC.copy()
        Layer.residues.resnames = 'PCL'
        AllLayers += Layer.copy()
    
    AllLayers.write(grofile)
    
    if vmd_visualization:
        with open('viewLayers.vmd','w') as f:
            print("""# Use this as:
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
mol representation Points 4.000000
mol color Name
mol selection {resname PCL}
mol material Opaque
mol addrep top""",file=f)
    
    return None


def compact_atoms(box,atomsQM,atomsMM,atomsPC):

    # Ref position is the center_of_geometry of QM layer
    Rref = atomsQM.center_of_geometry()
    
    # Iterate over all atoms
    for atom in atomsQM.atoms+atomsMM.atoms+atomsPC.atoms:
        for ix in range(3):
            d = atom.position[ix] - Rref[ix]
            if abs(d) > box[ix]/2:
                v = atom.position
                if (d>0):
                    v[ix] -= box[ix]
                else:
                    v[ix] += box[ix]
                atom.position =  v
                    
    return None
        


if __name__ == "__main__":
    import argparse

    # Input parser. Set flags
    parser = argparse.ArgumentParser(description='Generate ONIOM inputs from Gromacs trajectory.')
    parser.add_argument('-f',metavar='file.trr',help='Trajectory file',required=True)
    parser.add_argument('-s',metavar='file.tpr',help='Binary topoly file',required=True)
    parser.add_argument('-selQM',metavar='select_string',help='Selection command for the QM layer',default=None)
    parser.add_argument('-selMM',metavar='select_string',help='Selection command for the MM layer',default=None)
    parser.add_argument('-selPC',metavar='select_string',help='Selection command for the PC layer',default=None)
    parser.add_argument('-indQM',metavar='file.ndx',help='Index file for the QM layer',default=None)
    parser.add_argument('-indMM',metavar='file.ndx',help='Index file for the MM layer',default=None)
    parser.add_argument('-indPC',metavar='file.ndx',help='Index file for the PC layer',default=None)
    parser.add_argument('-ob',help='Base name of generated Gaussian input',default='snap')
    parser.add_argument('-osfx',help='Suffix to be appended to Gaussian input name',default='')
    parser.add_argument('-nzero',help='Digits for the label with the step number',type=int,default=3)
    parser.add_argument('-b',metavar='<time>',help='First frame (ps) to read from trajectory',type=float,default=-1.)
    parser.add_argument('-e',metavar='<time>',help='Last frame (ps) to read from trajectory',type=float,default=-1.)
    parser.add_argument('-dt',metavar='<time>',help='Only use frame when t MOD dt = first time (ps)',type=float,default=-1.)
    parser.add_argument('-method',help='Whole route section (e.g. "#p hf/sto-3g")',default='#p hf/sto-3g')
    parser.add_argument('-chargemult',help='Charge and multiplicity line (e.g. "0 1 0 1")',default='0 1')
    parser.add_argument('-FF',metavar='file.prm',help='FF file into be added to Gaussian input',default=None)
    parser.add_argument('-writeGRO',action='store_true',help='Write gro file to check layers',default=False)
    parser.add_argument('-compact',action='store_true',help='Process the snapshot to ensure compact representation',default=False)
    # Parse input
    args = parser.parse_args()

    # Get topology and coordinates
    u = MDAnalysis.Universe(args.s,args.f)

    # Selections
    # -- QM layer --
    if args.selQM:
        try:
            layerQM = u.select_atoms(args.selQM,updating=True)
        except:
            raise BaseException("Error setting QM layer. Maybe due to missplells in selection keyword")
    elif args.indQM:
        atoms_num = read_ndx(args.indQM)
        # Build selection kwd
        selection = 'bynum '
        for num in atoms_num:
            selection += str(num)+' '
        # Apply selection (static)
        try:
            layerQM = u.select_atoms(selection)
        except:
            raise BaseException("Error setting QM layer. Check index file")
    else:
        raise BaseException("No QM layer was set.")

    # -- MM layer --
    if args.selMM:
        try:
            # Ensure no overlap with QM layer
            selection='('+args.selMM+') and not group selQM'
            layerMM = u.select_atoms(selection,updating=True,selQM=layerQM,periodic=True)
        except:
            raise BaseException("Error setting MM layer. Maybe due to missplells in selection keyword")
    elif args.indMM:
        atoms_num = read_ndx(args.indMM)
        # Build selection kwd
        selection = 'bynum '
        for num in atoms_num:
            selection += str(num)+' '
        # Apply selection (static)
        try:
            # Ensure no overlap with QM layer
            selection='('+selection+') and not group selQM'
            layerMM = u.select_atoms(selection,selQM=layerQM)
        except:
            raise BaseException("Error setting MM layer. Check index file")
    else:
        # return an empty AtomGroup (acts as None type within if clauses)
        layerMM = MDAnalysis.AtomGroup([],u)

    # -- Point charges layer --
    if args.selPC:
        try:
            # Ensure no overlap with QM/MM layers
            selection = '('+args.selPC+') and not ((group selQM) or (group selMM))'
            layerPC = u.select_atoms(selection,updating=True,selQM=layerQM,selMM=layerMM,periodic=True)
        except:
            raise BaseException("Error setting PC layer. Maybe due to missplells in selection keyword")
    elif args.indPC:
        atoms_num = read_ndx(args.indPC)
        # Build selection kwd
        selection = 'bynum '
        for num in atoms_num:
            selection += str(num)+' '
        # Apply selection (static)
        try:
            # Ensure no overlap with QM/MM layers
            selection = '('+selection+') and not ((group selQM) or (group selMM))'
            layerPC = u.select_atoms(selection,selQM=layerQM,selMM=layerMM)
        except:
            raise BaseException("Error setting PC layer. Check index file")
    else:
        layerPC = MDAnalysis.AtomGroup([],u)

    
    # Set some defaults to slice the traj
    if args.b==-1.:
        tini = u.trajectory.time    
    else:
        tini = args.b
    if args.e==-1.:
        tfin = u.trajectory.totaltime
    else:
        tfin = args.e
    if args.dt==-1.:
        dt = u.trajectory.dt
    else:
        dt = args.dt

    istp=0
    for conf in u.trajectory:
        # Determine whether using the frame or not
        if u.trajectory.n_frames == 1:
            # Always use the frame when there is only one
            pass
        elif conf.time>tfin+0.0001:
            sys.exit()
        elif conf.time%dt > 0.0001:
            continue
        elif conf.time<tini:
            continue
        
        if args.compact:
            # Get box dimensions (assume rect box) 
            box = u.dimensions[:3]
            compact_atoms(box,layerQM,layerMM,layerPC)

        # Write Gaussian input
        fmt='%%s%%0%ig%%s.%%s'%args.nzero
        fname  = fmt%(args.ob,istp,args.osfx,'com')
        chkname= fmt%(args.ob,istp,args.osfx,'chk')
        f = open(fname,'w')
        write_oniom(layerQM,layerMM,layerPC,
                    unit=f,
                    title='Trajectory step %s (time=%s ps)'%(conf.frame,round(conf.time,5)),
                    chk=chkname,
                    FFfile=args.FF,
                    method=args.method,
                    chargemult=args.chargemult)
        f.close()
        files_gen = fname
        if args.writeGRO:
            grofile = fmt%(args.ob,istp,args.osfx,'gro')
            write_oniom_gro(layerQM,layerMM,layerPC,grofile,vmd_visualization=True)
            files_gen += ' '+grofile+' viewLayers.vmd'
        # Indicate all files generated
        print('%s generated'%files_gen)
        # Prepare next iteration
        istp += 1

