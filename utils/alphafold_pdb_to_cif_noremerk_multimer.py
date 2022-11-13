#!/usr/bin/python3
#
# File:    alphafold_pdb_to_cif.py
# Author:  Brinda Vallat
# Date:    Nov 29, 2021
# Version: 0.001
# Adapted from Ben Webb's modbase_pdb_to_cif.py available at https://github.com/salilab/modbase_utils
# Note: Only works for single chain models
# Input: PDB file with HEADER, TITLE, and REMARKS, and AlphaFold Pickle File
# Output: mmCIF file
# Command line arguments: Input AF pickle filename, Input PDB filename, Output mmCIF filename

import re
import collections
import pickle

# Reference sequence database
SequenceDB = collections.namedtuple(
    "SequenceDB", ["name", "code", "accession"])


# Mapping between one-letter codes and PDB names
three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

one_to_three = {val: key for key, val in three_to_one.items()}



class AFPickle:
    """Data from an AlphaFold Pickle File"""

    def __init__(self, fname):
        pickle_off = open(fname, "rb")
        pdata = pickle.load(pickle_off)
        self.plddt = pdata['plddt']
        self.pae = pdata['predicted_aligned_error']

class CifLoop:
    """Helper class to write an mmCIF loop construct"""

    def __init__(self, fh, category, keys):
        self.fh, self.category, self.keys = fh, category, keys
        self._empty_loop = True

    def write(self, line):
        f = self.fh
        if self._empty_loop:
            f.write("#\nloop_\n")
            for k in self.keys:
                f.write("_%s.%s\n" % (self.category, k))
            self._empty_loop = False
        print(line, file=f)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if not self._empty_loop:
            self.fh.write("#\n")


class CifWriter:
    def __init__(self, fh, pkl):
        class CifID:
            pass

        self.fh, self.pkl = fh, pkl 
        # Assign consecutive CIF numeric IDs from 1
        self.target = CifID()
        self.coord = CifID()
        self.alignment = CifID()
        self.target.entity_id = 1
        self.target.data_id = 1
        self.alignment.data_id, self.coord.data_id = 2, 3

    def print(self, s):
        print(s, file=self.fh)

    def loop(self, category, keys):
        return CifLoop(self.fh, category, keys)

    def write_header(self, model_id, title):
        self.print("data_ %s" % model_id)
#        comment-out because these need additional REMARK files
#        self.print("_entry.id %s" % model_id)
#        self.print("_struct.entry_id %s" % model_id)
#        if title:
#            self.print("_struct.title '%s'" % title)

    def write_exptl(self, model_id, expdta):
        if expdta.startswith('THEORETICAL MODEL'):
            details = "Computational model generated using AlphaFold2"
            self.print("#\n_exptl.entry_id %s" % model_id)
            self.print("_exptl.method 'THEORETICAL MODEL'")
            self.print("_exptl.details '%s'" % details)

    def write_audit_conform(self):
        with self.loop("audit_conform", ["dict_name", "dict_version", "dict_location"]) as lp:
            lp.write("mmcif_ma.dic  1.3.3  https://github.com/ihmwg/MA-dictionary/blob/master/mmcif_ma.dic")
            lp.write("mmcif_pdbx.dic  5.342  http://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic")

    def write_audit_author(self):
        with self.loop("audit_author", ["name", "pdbx_ordinal"]) as lp:
            lp.write("'Doe J' 1")
            lp.write("'Doe R' 2")

    def write_citation(self):
        with self.loop(
                "citation",
                ["id", "title", "journal_abbrev", "journal_volume",
                 "page_first", "page_last", "year", "pdbx_database_id_PubMed",
                 "pdbx_database_id_DOI"]) as lp:
            # Primary citation
            lp.write("primary . . . . . . . .")
            # AlphaFold citation
            lp.write("2 'Highly accurate protein structure prediction with AlphaFold"
                     "'\n'Nature' 596 583 589 2021 "
                     "34265844 10.1038/s41586-021-03819-2")
        with self.loop(
                "citation_author", ["citation_id", "name", "ordinal"]) as lp:
            lp.write("2 'Jumper J' 1")
            lp.write("2 'Evans R' 2")
            lp.write("2 'Pritzel A' 3")
            lp.write("2 'Green T' 4")
            lp.write("2 'Figurnov M' 5")
            lp.write("2 'Ronneberger O' 6")
            lp.write("2 'Tunyasuvunakool K' 7")
            lp.write("2 'Bates R' 8")
            lp.write("2 'Zidek A' 9")
            lp.write("2 'Potapenko A' 10")
            lp.write("2 'Bridgland A' 11")
            lp.write("2 'Meyer C' 12")
            lp.write("2 'Kohl SAA' 13")
            lp.write("2 'Ballard AJ' 14")
            lp.write("2 'Cowie A' 15")
            lp.write("2 'Romera-Paredes  B' 16")
            lp.write("2 'Nikolov S' 17")
            lp.write("2 'Jain R' 18")
            lp.write("2 'Adler J' 19")
            lp.write("2 'Back T' 20")
            lp.write("2 'Petersen S' 21")
            lp.write("2 'Reiman D' 22")
            lp.write("2 'Clancy E' 23")
            lp.write("2 'Zielinski M' 24")
            lp.write("2 'Steinegger M' 25")
            lp.write("2 'Pacholska M' 26")
            lp.write("2 'Berghammer T' 27")
            lp.write("2 'Bodenstein S' 28")
            lp.write("2 'Silver D' 29")
            lp.write("2 'Vinyals O' 30")
            lp.write("2 'Senior AW' 31")
            lp.write("2 'Kavukcuoglu K' 32")
            lp.write("2 'Kohli P' 33")
            lp.write("2 'Hassabis D' 34")
            lp.write("primary 'Doe J' 35")
            lp.write("primary 'Doe R' 36")

    def write_software(self):
        alphafold_version = "2.0.1"
        with self.loop(
                "software",
                ["pdbx_ordinal", "name", "classification",
                 "version", "type", "location", "citation_id"]) as lp:
            lp.write("1 AlphaFold 'model building' %s program "
                     "https://github.com/deepmind/alphafold 2" % alphafold_version)

        # Put each piece of software in its own group
        with self.loop(
                "ma_software_group",
                ["ordinal_id", "group_id", "software_id"]) as lp:
            for i in range(1, 2):
                lp.write("%d %d %d" % (i, i, i))

    def write_chem_comp(self):
        # Just assume all 20 standard amino acids are in the model
        with self.loop(
                "chem_comp",
                ["id", "type", "name", "formula", "formula_weight"]) as lp:
            lp.write("""ALA 'L-peptide linking' ALANINE 'C3 H7 N O2' 89.094
ARG 'L-peptide linking' ARGININE 'C6 H15 N4 O2 1' 175.212
ASN 'L-peptide linking' ASPARAGINE 'C4 H8 N2 O3' 132.119
ASP 'L-peptide linking' 'ASPARTIC ACID' 'C4 H7 N O4' 133.103
CYS 'L-peptide linking' CYSTEINE 'C3 H7 N O2 S' 121.154
GLN 'L-peptide linking' GLUTAMINE 'C5 H10 N2 O3' 146.146
GLU 'L-peptide linking' 'GLUTAMIC ACID' 'C5 H9 N O4' 147.130
GLY 'peptide linking' GLYCINE 'C2 H5 N O2' 75.067
HIS 'L-peptide linking' HISTIDINE 'C6 H10 N3 O2 1' 156.165
ILE 'L-peptide linking' ISOLEUCINE 'C6 H13 N O2' 131.175
LEU 'L-peptide linking' LEUCINE 'C6 H13 N O2' 131.175
LYS 'L-peptide linking' LYSINE 'C6 H15 N2 O2 1' 147.198
MET 'L-peptide linking' METHIONINE 'C5 H11 N O2 S' 149.208
PHE 'L-peptide linking' PHENYLALANINE 'C9 H11 N O2' 165.192
PRO 'L-peptide linking' PROLINE 'C5 H9 N O2' 115.132
SER 'L-peptide linking' SERINE 'C3 H7 N O3' 105.093
THR 'L-peptide linking' THREONINE 'C4 H9 N O3' 119.120
TRP 'L-peptide linking' TRYPTOPHAN 'C11 H12 N2 O2' 204.229
TYR 'L-peptide linking' TYROSINE 'C9 H11 N O3' 181.191
VAL 'L-peptide linking' VALINE 'C5 H11 N O2' 117.148""")

    def write_entity_details(self, sequence3_chain, chain_name_list, chain_ids, genes):
        # entities for target
        with self.loop(
                "entity",
                ["id", "type", "src_method", "pdbx_description"]) as lp:
            for i in range(len(genes)):
                lp.write("%d polymer man %s" % (i+1, genes[i]))

        with self.loop(
                "entity_poly",
                ["entity_id", "type", "nstd_linkage",
                 "pdbx_seq_one_letter_code",
                 "pdbx_seq_one_letter_code_can"]) as lp:
            i=0
            for chain_id in chain_ids:
                target_primary = "".join(three_to_one[x] for x in sequence3_chain[i])
                lp.write("%d polypeptide(L) no %s %s"
                         % (i+1, target_primary, target_primary))
                i=i+1

        with self.loop(
                "entity_poly_seq",
                ["entity_id", "num", "mon_id", "hetero"]) as lp:
            for j in range(len(chain_ids)):
                for i, s in enumerate(sequence3_chain[j]):
                    lp.write("%d %d %s ." % (j+1, i+1, s))

    def write_target_details(self, chain_ids, sequence3, seqdb, target_begin, target_end):
        with self.loop(
                "ma_target_entity", ["entity_id", "data_id", "origin"]) as lp:
            di = 1
            for chain_id in chain_ids:
                lp.write("%d %d ." % (di, di))
                di+=1

        with self.loop(
                "ma_target_entity_instance",
                ["asym_id", "entity_id", "details"]) as lp:
            ei=1
            for chain_id in chain_ids:
                lp.write("%s %d ." % (chain_id, ei))
                ei+=1

        with self.loop(
                "ma_target_ref_db_details",
                ["target_entity_id", "db_name", "db_name_other_details",
                 "db_code", "db_accession", "seq_db_isoform",
                 "seq_db_align_begin", "seq_db_align_end"]) as lp:
            for db in seqdb:
                if db.name == 'UniProt':
                    lp.write("%d UNP . %s %s ? %s %s"
                             % (self.target.entity_id, db.code, db.accession,
                                target_begin, target_end))
                elif db.name == 'RefSeq':
                    lp.write("%d %s . %s %s ? %s %s"
                             % (self.target.entity_id, db.name, db.code,
                                db.accession, target_begin, target_end))

    def write_assembly(self, chain_ids, sequence3_chain):
        with self.loop(
                'ma_struct_assembly',
                ['ordinal_id', 'assembly_id', 'entity_id', 'asym_id',
                 'seq_id_begin', 'seq_id_end']) as lp:
            ei = 0
            for chain_id in chain_ids:                        
                l = len(sequence3_chain[ei])
                ei += 1
                lp.write("1 1 %d %s 1 %d"
                         % (ei, chain_id, l))

    def write_data(self):
        with self.loop("ma_data", ["id", "name", "content_type"]) as lp:
            lp.write("%d 'Target Sequence' target" % self.target.data_id)
            lp.write("%d 'Coevolution MSA' "
                     "'coevolution MSA'" % self.alignment.data_id)
            lp.write("%d 'Target Structure' 'model coordinates'"
                     % self.coord.data_id)

        # Put each data item in its own group
        with self.loop(
                "ma_data_group", ["ordinal_id", "group_id", "data_id"]) as lp:
            for i in range(1, 4):
                lp.write("%d %d %d" % (i, i, i))

    def write_protocol(self):
        with self.loop(
                'ma_protocol_step',
                ['ordinal_id', 'protocol_id', 'step_id', 'method_type',
                 'step_name', 'software_group_id']) as lp:
            lp.write("1 1 1 'coevolution MSA' 'AlphaFold' 1")
            lp.write("2 1 2 'modeling' 'AlphaFold' 1")
            lp.write("3 1 3 'model selection' 'AlphaFold' 1")

    def write_model_list(self):
        with self.loop(
                'ma_model_list',
                ['ordinal_id', 'model_id', 'model_group_id', 'model_name',
                 'model_group_name', 'assembly_id', 'data_id',
                 'model_type']) as lp:
            lp.write("1 1 1 'Top scoring model' . 1 %d 'Other'" % self.coord.data_id)

    def write_asym(self, chain_ids):
        with self.loop('struct_asym', ['id', 'entity_id', 'details']) as lp:
            ei = 1
            for chain_id in chain_ids:
                lp.write("%s %d ?" % (chain_id, ei))
                ei += 1

    def write_seq_scheme(self, chain_ids, sequence3, tgtbeg, tgtend):
        assert len(sequence3) == tgtend - tgtbeg + 1
        entity_id = self.target.entity_id
        with self.loop(
            'pdbx_poly_seq_scheme',
            ['asym_id', 'entity_id', 'seq_id', 'mon_id', 'pdb_seq_num',
             'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
             'pdb_strand_id']) as lp:
            for i, s in enumerate(sequence3):
                seqid = 1 + i
                auth_seqid = tgtbeg + i
                lp.write("%s %d %d %s %d %d %s %s %s" % (chain_ids[seqid-1], entity_id, seqid, s, auth_seqid, auth_seqid, s, s, chain_ids[seqid-1]))

    def write_atom_site(self, chain_id, atoms, resnum_begin, resnum_end):
        elements = set()
        auth_seqid = resnum_begin
        entity_id = 0
        seqid = 1
        ordinal = 1
        model_num = 1
        pdb_resnum = None
        this_chain = None
        chain=None
        with self.loop(
            'atom_site',
            ['group_PDB', 'type_symbol', 'label_atom_id',
             'label_comp_id', 'label_asym_id', 'label_seq_id',
             'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_asym_id',
             'label_alt_id', 'Cartn_x', 'Cartn_y', 'Cartn_z',
             'occupancy', 'B_iso_or_equiv', 'label_entity_id', 'pdbx_PDB_model_num', 'id']) as lp:
            for a in atoms:
                # Detect new residue if PDB resnum changed
                pdb_this_resnum = a[22:26]
                if pdb_resnum is not None and pdb_this_resnum != pdb_resnum:
                    auth_seqid += 1    
                    seqid += 1

                this_chain = a[21]
                if this_chain is not None and this_chain!=chain:
                    seqid = 1
                    entity_id += 1
                        
                chain = this_chain
                pdb_resnum = pdb_this_resnum
                inscode = a[26:27].strip() or '?'
                group_pdb = a[:6]
                element = a[76:78].strip() or '?'
                elements.add(element)
                atmnam = a[12:16]
                resnam = a[17:20]
                x = a[30:38]
                y = a[38:46]
                z = a[46:54]
                occ = a[54:60]
                tfac = a[60:66]
                cid = a[21]
                lp.write("%s %s %s %s %s %d %d %s %s . %s %s %s %s %s %d %d %d"
                         % (group_pdb, element, atmnam, resnam, cid,
                            seqid, auth_seqid, inscode, cid, x, y, z,
                            occ, tfac, entity_id, model_num, ordinal))
                ordinal += 1

        with CifLoop(fh, 'atom_type', ['symbol']) as lp:
            for element in sorted(elements):
                lp.write(element)

    def write_af_scores(self, chain_name_list, sequence3):
        if not self.pkl:
            return

        assert len(sequence3) == len(self.pkl.plddt), "Length of sequence in the PDB and pLDDT in AlphaFold Pickle file do not match"

        paeLenExp = int(len(sequence3)**2)
        paeLen = int(len(self.pkl.pae)*len(self.pkl.pae[0]))

        assert paeLenExp == paeLen, "Length of sequence in the PDB and PAE in AlphaFold Pickle file do not match"

        mode = ('global', 'local', 'local-pairwise')
        name = ('pLDDT', 'pLDDT', 'PAE')
        type1 = ('pLDDT', 'pLDDT', 'PAE')
        software_group_id = 1

        with self.loop(
            'ma_qa_metric',
            ['id', 'mode', 'name', 'software_group_id', 'type']) as lp:
            for id1 in range(3):
                lp.write("%d %s %s %d %s"
                        % (id1+1, mode[id1], name[id1], software_group_id, type1[id1]))
            
        ordinal_id = 1
        model_id = 1
        glddtsum = 0

        with self.loop(
            'ma_qa_metric_local',
            ['ordinal_id', 'model_id', 'label_asym_id', 'label_seq_id', 'label_comp_id', 'metric_id', 'metric_value']) as lp:
            for seq in sequence3:
                label_comp_id = seq
                label_seq_id = ordinal_id
                metric_id = 2
                metric_value = self.pkl.plddt[ordinal_id-1]
                glddtsum = glddtsum + metric_value
                lp.write("%d %d %s %d %s %d %.2f"
                         % (ordinal_id, model_id, chain_name_list[ordinal_id-1], label_seq_id, label_comp_id, metric_id, metric_value))
                ordinal_id += 1

        glddt = glddtsum/len(self.pkl.plddt)
        metric_id = 1
        ordinal_id = 1

        with self.loop(
            'ma_qa_metric_global',
            ['ordinal_id', 'model_id', 'metric_id', 'metric_value']) as lp:
            lp.write("%d %d %d %.2f"
                    % (ordinal_id, model_id, metric_id, glddt))

        ordinal_id = 1

        with self.loop(
            'ma_qa_metric_local_pairwise',
            ['ordinal_id', 'model_id', 'label_asym_id_1', 'label_seq_id_1', 'label_comp_id_1', 'label_asym_id_2', 'label_seq_id_2', 'label_comp_id_2', 'metric_id', 'metric_value']) as lp:
            nseq1 = 0
            for seq1 in sequence3:
                nseq1 += 1
                nseq2 = 0
                for seq2 in sequence3:
                    nseq2 += 1
                    label_comp_id_1 = seq1
                    label_comp_id_2 = seq2
                    label_seq_id_1 = nseq1
                    label_seq_id_2 = nseq2
                    metric_id = 3
                    metric_value = self.pkl.pae[nseq1-1][nseq2-1]
                    lp.write("%d %d %s %d %s %s %d %s %d %.2f"
                             % (ordinal_id, model_id, chain_name_list[nseq1-1], label_seq_id_1, label_comp_id_1, chain_name_list[nseq2-1], label_seq_id_2, label_comp_id_2, metric_id, metric_value))
                    ordinal_id +=1  
        
class Structure:
    """Handle read of PDB structure and write of mmCIF"""

    def _read_pdb(self, fh):
        self.remarks = {}
        self.seqdb = []
        self.target = []
        self.expdta = None
        self.title = None
        self.atoms = []

        for line in fh:
            # Handle PDB headers
            if line.startswith('REMARK 220 SEQDB:'):
                val = [x.strip() for x in line[17:].split()]
                if len(val) == 3 and val[0] == 'UniProt':
                    self.seqdb.append(SequenceDB(
                        name=val[0], accession=val[1], code=val[2]))
                elif len(val) == 2 and val[0] == 'RefSeq' and '.' in val[1]:
                    self.seqdb.append(SequenceDB(
                        name=val[0], accession=val[1].split('.', 1)[0],
                        code=val[1]))
            elif line.startswith('REMARK') and line.count(':') == 1:
                key, val = [x.strip() for x in line[11:].split(':')]
                #print(key,val)
                self.remarks[key] = val
            elif line.startswith('TITLE     '):
                self.title = line[10:].strip()
            elif line.startswith('EXPDTA    '):
                self.expdta = line[10:].strip()
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                self.atoms.append(line)
        # Assume that all models are single chain
        if self.atoms:
            self.chain_ids = list(set(self.atoms[:][21].strip()))
            #print(self.chain_ids)

    def get_sequence3(self):
        """Get PDB sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                yield a[17:20].strip()  # residue name
                resnum = this_resnum

    def get_sequence3_chain(self,chain_id):
        """Get PDB sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                if a[21] == chain_id:
                    yield a[17:20].strip()  # residue name
                    resnum = this_resnum

    def get_chains(self):
        """Get PDB sequence as a sequence of 3-letter residue names"""
        resnum = None
        for a in self.atoms:
            this_resnum = a[22:26]  # residue sequence number
            if this_resnum != resnum:
                yield a[21].strip()  # residue name
                resnum = this_resnum

    def write_mmcif(self, fh, pkl):
        """Write current structure out to a mmCIF file handle"""
        # mmCIF models must always have a chain ID
        sequence3 = list(self.get_sequence3())
        chain_name_list= list(self.get_chains())
        chain_ids = sorted(list(set(chain_name_list)))
        sequence3_chain=list()
        for chain_id in chain_ids:
            sequence3_chain.append(list(self.get_sequence3_chain(chain_id)))

        pkl = AFPickle(pkl)

        c = CifWriter(fh, pkl)
        c.write_header("", self.title)
        c.write_audit_conform()
        c.write_audit_author()
        c.write_citation()
        c.write_software()
        c.write_chem_comp()
        # pad description by chain_ids
        c.write_entity_details(sequence3_chain, chain_name_list,chain_ids,chain_ids)
        tgtbeg = 1
        tgtend = len(sequence3)
        c.write_target_details(chain_ids, sequence3, self.seqdb, tgtbeg, tgtend)
        c.write_assembly(chain_ids, sequence3_chain)
        c.write_data()
        c.write_protocol()
        c.write_model_list()
        c.write_asym(chain_ids)
        c.write_seq_scheme(chain_name_list, sequence3, tgtbeg, tgtend)
        c.write_atom_site(chain_ids, self.atoms, tgtbeg, tgtend)
        c.write_af_scores(chain_name_list, sequence3)

def read_pdb(fh):
    """Read PDB file from filehandle and return a new Structure"""
    s = Structure()
    s._read_pdb(fh)
    return s

if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser(
            description="Utility to convert AlphaFold PDB files to mmCIF",
            epilog="""
Convert a PDB file, generated by AlphaFold, to mmCIF format. This requires
the PDB file with HEADER, TITLE, and REMARKS, and the corresponding pickle
file with pLDDT scores (key: 'plddt') of shape (N,) and PAE scores 
(key: 'predicted_aligned_error') of shape (N,N), where N is the total number 
of residues in the corresponding PDB file. 
""")
    a.add_argument("pkl", help="Input AF pickle file")
    a.add_argument("pdb", help="Input PDB file")
    a.add_argument("mmcif", help="Output mmCIF file")
    args = a.parse_args()

    with open(args.pdb) as fh:
        s = read_pdb(fh)
    with open(args.mmcif, 'w') as fh:
        s.write_mmcif(fh, args.pkl)
