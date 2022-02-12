# modelcif_utils

This repository contains scripts that support the [ModelCIF](https://github.com/ihmwg/ModelCIF) 
data representation. 

[ModelCIF](https://github.com/ihmwg/ModelCIF) is an extension of the [PDBx/mmCIF](http://mmcif.wwpdb.org) 
dictionary. 

## Scripts
 - `utils/alphafold_pdb_to_cif.py` takes a PDB file and an AlphaFold pickle file (with pLDDT and PAE scores) to
   generate an mmCIF file compliant with [ModelCIF](https://github.com/ihmwg/ModelCIF). The script relies
   on metadata information (Title, reference sequence information etc.) in the header section of the PDB file. 
   For an example of PDB input file header information, see [here](utils/P69905_HEADER.pdb). 
