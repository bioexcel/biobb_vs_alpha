#!/usr/bin/env python

"""Python wrapper module for the GROMACS solvate module
"""
import sys
import json
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu

import os
import warnings
import shutil
from mmb_api import pdb
import Bio.PDB

import numpy
import fileinput
import re

from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as prot_one_letter


#from Bio.Struct.Geometry import center_of_mass
#from Bio import Struct


class Box(object):
    """Sets the center and the size of a rectangular parallelepiped box around a selection of residues found in a given PDB. The residue identifiers that compose the selection (i.e. binding site) are extracted from a second PDB.
    Args:
        input_pdb_path (str):  PDB protein structure for which the box will be build. Its size and center will be set around the 'resid_pdb_path' residues once mapped against this PDB.
        resid_pdb_path (str):  PDB file containing a selection of residue numbers mappable to 'input_pdb_path'.
        output_pdb_path (str): PDB protein structure coordinates including the annotation of the box center and size as REMARKs.
        offset:                Extra distance (Amstrongs) between the last residue atom and the box boundary. Default: 2
    """

    def __init__(self, input_pdb_path, resid_pdb_path, output_pdb_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        #i/o
        self.input_pdb_path       = input_pdb_path
        self.resid_pdb_path       = resid_pdb_path
        self.output_pdb_path      = output_pdb_path
        self.offset               = properties.get('offset', 2)

    def __align_sequences(self, seqA, seqB):
        """
        Performs a global pairwise alignment between two sequences using the Needleman-Wunsch algorithm as implemented in Biopython.
        Returns the alignment and the residue mapping between both original sequences.
        """

        #seq list to seq string
        sequence_A = ''.join([i[1] for i in seqA])
        sequence_B = ''.join([i[1] for i in seqB])

        # Do pairwaise alignment
        alns = Bio.pairwise2.align.globalds(sequence_A, sequence_B,self.matrix, self.gap_open, self.gap_extend, penalize_end_gaps=(False, False) )


        best_aln = alns[0]
        aligned_A, aligned_B, score, begin, end = best_aln

        # Equivalent residue numbering. Relative to reference
        mapping = {}
        aa_i_A, aa_i_B = 0, 0
        for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
            if aa_aln_A == '-':
                if aa_aln_B != '-':
                    aa_i_B += 1
            elif aa_aln_B == '-':
                if aa_aln_A != '-':
                    aa_i_A += 1
            else:
                assert seqA[aa_i_A][1] == aa_aln_A
                assert seqB[aa_i_B][1] == aa_aln_B
                mapping[seqA[aa_i_A][0]] = seqB[aa_i_B][0]
                aa_i_A += 1
                aa_i_B += 1

        return ((aligned_A, aligned_B), mapping)

    def __get_pdb_sequence(seft, structure):
        """
        Retrieves the AA sequence from a PDB structure.
        """

        aa = lambda r: (r.id[1], prot_one_letter.get(r.resname, 'X'))
        seq=[]
        for r in structure.get_residues():
            if Bio.PDB.Polypeptide.is_aa(r):
                seq.append(aa(r))
        return seq


    def launch(self):
        """Launches the pipeline to build a box around a selection of residues
        """
        #out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)

        ##
        ## Loading and parsing reference PDB structure

        parser = Bio.PDB.PDBParser()

        # Parse input structure
        print "Loading input PDB structure %s..." % self.input_pdb_path
        structure_name = os.path.basename(self.input_pdb_path.split('.')[0])
        structPDB      = parser.get_structure(structure_name,self.input_pdb_path)[0]
        structPDB_seq = [three_to_one(res.get_resname()) for res in structPDB.get_residues() if is_aa(res.get_resname(), standard=True)]
        print structPDB_seq

        # Parse residue structure
        print "Loading residue PDB selection %s..." % self.resid_pdb_path
        resid_name   = os.path.basename(self.resid_pdb_path.split('.')[0])
        residPDB     = parser.get_structure(resid_name,self.resid_pdb_path)[0]
        residPDB_seq = [three_to_one(res.get_resname()) for res in residPDB.get_residues() if is_aa(res.get_resname(), standard=True)]
        print residPDB_seq

        ##
        ## Mapping residue structure into input structure

        # Listing residues to be selected from the residue structure
        # residPDB_res_list = []
        # p = re.compile('H_|W_')
        # for residPDB_res in residPDB.get_residues():
        #     m_het = p.match(residPDB_res.get_id()[0])
        #     if not m_het:
        #         residPDB_res_list.append(residPDB_res.get_id())
        #         binding_site_CA_list.append(residPDB_res['CA'])
        

        # Aligning
        # alignments = pairwise2.align.localxx("".join(structPDB_seq), "".join(residPDB_seq))
        # print alignments[0]
        # mappings = Bio.PDB.StructureAlignment(alignments[0], structPDB, residPDB).get_maps()
        # print mappings
        # # Mapping selected residues to input structure
        # selection_res_list   = []
        # selection_atoms_num  = 0
        # for struct_chain in structPDB:
        #     for struct_res in struct_chain:
        #         if struct_res.get_id() in residPDB_res_list:
        #             selection_res_list.append(struct_res)
        #             selection_atoms_num += len(struct_res.get_list())

         # Get AA sequence
        clusterPDB_seq = self.__get_pdb_sequence(clusterPDB)

        # Pairwise align
        aln, residue_map = self.__align_sequences(structPDB_seq,clusterPDB_seq)
        print residue_map
        print "    Matching residues to input PDB structure. Alignment is:\n    %s" % aln[1]

        # Calculate (gapless) sequence identity
        seq_identity, gap_seq_identity = self.__calculate_alignment_identity(aln[0], aln[1])
        print "    Sequence identity (%%): %s" % seq_identity
        print "    Gap less identity (%%): %s" % gap_seq_identity


        ##
        ## Selecting aligned CA atoms from first model, first chain

        struct_atoms  = []
        cluster_atoms = []

        for struct_res in residue_map:
            try:
                cluster_atoms.append(clusterPDB[residue_map[struct_res]]['CA'])
                struct_atoms.append(structPDB[struct_res]['CA'])
            except KeyError:
                print "Cannot find CA atom for residue %s  (input PDB  %s)" % (structPDB[struct_res],struct_res)
                pass

        if len(cluster_atoms)==0:
            raise Exception('Cannot find CA atoms (1st model, 1st chain) in cluster member {1} when aligning against {2}. Ignoring this member.'.format(clusterPDB_path,structure_name))
        else:
            print "    Superimposing %s aligned protein residues" % len(cluster_atoms)

        # Align against input structure
        si = Bio.PDB.Superimposer()
        si.set_atoms(struct_atoms, cluster_atoms)
        si.apply(clusterPDB.get_atoms())
        print "    RMSD: %s" %si.rms

        # Save transformed structure (and ligand)
        clusterPDB_ligand_aligned = clusterPDB[clusterPDB_ligand.get_id()]
        print "    Saving transformed ligand coordinates"

        clusterPDB_ligands_aligned.append(clusterPDB_ligand_aligned)

        if len(selection_res_list) == 0:
            raise Exception('Cannot match any of the residues listed in %s into %s' % (self.resid_pdb_path,self.input_pdb_path) )
        elif len(selection_res_list) !=  len(residPDB_res_list):
            warnings.warn('Cannot match all the residues listed in %s into %s. Found %s out of %s'  % (self.resid_pdb_path,self.input_pdb_path,len(selection_res_list),len(residPDB_res_list)))
        else:
            print "Selection residues successfully matched"


        ##
        ## Compute binding site box size

        # compute box center
        selection_box_center = numpy.sum(atom.coord for res in selection_res_list for atom in res.get_atoms()) / selection_atoms_num
        print "Binding site center (Amstrongs): %8.3f%8.3f%8.3f" % (selection_box_center[1],selection_box_center[1],selection_box_center[2])

        # compute box size
        selection_coords_max = numpy.amax([atom.coord for res in selection_res_list for atom in res.get_atoms()],axis=0)
        selection_box_size   = selection_coords_max - selection_box_center
        if self.offset:
            selection_box_size = [c + self.offset for c in selection_box_size]
        print "Binding site size (Amstrongs):   %8.3f%8.3f%8.3f" % (selection_box_size[0],selection_box_size[1],selection_box_size[2])

        vol = numpy.prod(selection_box_size) * 2**3
        print "Volume (cubic Amstrongs): %.0f" % vol

        # add box details as PDB remarks
        #remarks  = "REMARK 900\nREMARK 900 RELATED  ENTRIES\nREMARK 900 RELATED ID:%s CHAIN:%s\n" % (self.pdb_code,self.pdb_chain)
        remarks = "REMARK BOX CENTER:%8.3f%8.3f%8.3f" % (selection_box_center[1],selection_box_center[1],selection_box_center[2])
        remarks += " SIZE:%8.3f%8.3f%8.3f" % (selection_box_size[0],selection_box_size[1],selection_box_size[2])

        # add (optional) box coordinates as 8 ATOM records
        #selection_box_coords_txt  = self.get_box_coordinates(selection_box_center,selection_box_size)
        selection_box_coords_txt   = ""

        # write output pdb

        shutil.copy2(self.input_pdb_path, self.output_pdb_path)

        with open(self.output_pdb_path, 'r+') as f:
            content = f.read()
            if "END" in content:
                content = content.replace("END", selection_box_coords_txt + "END")
            else:
                content += selection_box_coords_txt
            f.seek(0, 0)
            f.write(remarks.rstrip('\r\n') + '\n' + content)


        print "Output PDB file (with box setting annotations): %s" % self.output_pdb_path



    def get_box_coordinates(self, box_center, box_size, pdb_format=True):
        coords = [
            [box_center[0]-box_size[0],box_center[1]-box_size[1],box_center[2]-box_size[2]],
            [box_center[0]-box_size[0],box_center[1]-box_size[1],box_center[2]+box_size[2]],
            [box_center[0]-box_size[0],box_center[1]+box_size[1],box_center[2]-box_size[2]],
            [box_center[0]-box_size[0],box_center[1]+box_size[1],box_center[2]+box_size[2]],
            [box_center[0]+box_size[0],box_center[1]-box_size[1],box_center[2]-box_size[2]],
            [box_center[0]+box_size[0],box_center[1]-box_size[1],box_center[2]+box_size[2]],
            [box_center[0]+box_size[0],box_center[1]+box_size[1],box_center[2]-box_size[2]],
            [box_center[0]+box_size[0],box_center[1]+box_size[1],box_center[2]+box_size[2]]
        ]

        if pdb_format:
            coords_txt = ""
            at_num = 10000
            at_nam = "ZN"
            re_nam = "ZN"
            chain  = "Z"
            res_num= 9999
            occ    = 1
            bfact  = 50
            elem   = "ZN"
            for coord in coords:
                coords_txt += "HETATM%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (at_num,at_nam,re_nam,chain,res_num,coord[0],coord[1],coord[2],occ,bfact,elem)
                at_num +=1
            return coords_txt
        else:
            return coords



#Creating a main function to be compatible with CWL
def main():
    system         = sys.argv[1]
    step           = sys.argv[2]
    properties_file= sys.argv[3]
    prop           = settings.YamlReader(properties_file, system).get_prop_dic()[step]

    BindingSite(input_pdb_path   = sys.argv[4],
                resid_pdb_path   = sys.argv[5],
                output_pdb_path  = sys.argv[6],
                properties           = prop    ).launch()


if __name__ == '__main__':
    main()
