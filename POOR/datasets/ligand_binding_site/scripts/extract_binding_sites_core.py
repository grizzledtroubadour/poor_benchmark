#!/usr/bin/env python3
"""
Core functions for protein-ligand binding site extraction.

This module contains the parsing and binding-site detection logic, shared by
01_extract_binding_sites.py and the batch drivers (03/05); the former resume
driver 02_resume_extract_binding_sites.py is archived under output/scripts_archive/.
"""

import gzip
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# Standard amino acid 3-letter codes
AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}

# DNA/RNA residue names to skip when looking for protein chains
NUCLEIC_NAMES = {
    "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
    "RA", "RC", "RG", "RU", "DI", "CI", "GI", "TI", "UI",
}

# Water / deuterated water
WATER_NAMES = {"HOH", "DOD", "WAT", "H2O", "OH2", "OHH"}

# Common modified amino acids that are part of the protein backbone.
MODIFIED_AA_NAMES = {
    "MSE", "SEP", "TPO", "PTR", "MEN", "FME", "CME", "CSD",
    "CSW", "OCY", "KCX", "MLY", "M3L", "DM0", "DMH", "HYP",
    "NLN", "ALC", "SAC", "PCA", "MLZ", "HLY", "CSO", "LLP",
    "CXM", "FME", "MVA", "NLE", "ORN", "PFF", "TPQ", "TYS",
    "SME", "SNC", "YCM",
}

# Mapping from modified amino acid to one-letter standard amino acid
MODIFIED_AA_1LETTER = {
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y",
    "KCX": "K", "MLY": "K", "M3L": "K", "MLZ": "K",
    "HLY": "K", "FME": "M", "CME": "C", "CSD": "C",
    "CSW": "C", "OCY": "C", "CSO": "C", "HYP": "P",
    "NLN": "L", "NLE": "L", "ALC": "A", "SAC": "S",
    "PCA": "E", "CXM": "M", "MVA": "V", "ORN": "R",
    "PFF": "F", "TPQ": "Y", "TYS": "Y", "SME": "M",
    "SNC": "C", "YCM": "C", "MEN": "N", "DM0": "K",
    "DMH": "K", "LLP": "K",
}

# Common crystallization additives / non-biological small molecules
CRYSTALLIZATION_ADDITIVE_NAMES = {
    # Alcohols / polyols
    "GOL", "EDO", "MPD", "PG4", "PGE", "PEG", "PE3", "PE4", "PE5", "PE6",
    "PE7", "PE8", "P33", "P34", "P35", "P36", "P37", "P38", "P39", "P40",
    "P6G", "PGT", "XLT", "1PE", "2PE", "3PE", "4PE", "SOG", "BOG", "BNG",
    "BMA", "MAN", "GLC", "GAL", "FUC", "RIB", "XYS", "FRU", "SUC", "TRE",
    "MAL", "CEL", "LCT", "NAG", "NDG", "BDP",
    # Acids / anions
    "SO4", "SUL", "PO4", "PHO", "ACT", "ACE", "ACA", "ACI", "CIT", "TAR",
    "MLA", "OXL", "FOR", "FMT", "BEN", "SUC", "MAL", "GLU", "ASP",
    # Buffers
    "TRS", "MES", "HEP", "EPP", "MPO", "TBU", "TAM", "TEA", "TFA", "HEZ",
    # Solvents
    "DMS", "DMSO", "DMF", "DIO", "EOH", "ETOH", "IPA", "MOH", "ACN", "THF",
    "BME", "BET", "DTT", "BUA", "BU1", "PDO",
    # Ions
    "CL", "NA", "K", "CA", "MG", "ZN", "FE", "CO", "NI", "CU", "MN",
    "CD", "HG", "PB", "SR", "BA", "RB", "CS", "AL", "CR", "MO", "W",
    "SE", "BR", "IOD",
    # Other common additives
    "AZI", "CYN", "SCN", "PER", "NO3", "NO2", "NH2", "NH4", "AMO", "AMM",
    "UNX", "UNL", "UNK",
}

HEAVY_ATOM_NAMES = {"H", "D"}


def is_heavy_atom(atom):
    """Return True if atom is not H/D."""
    elem = atom.element.strip().upper() if atom.element else ""
    if elem:
        return elem not in HEAVY_ATOM_NAMES
    name = atom.name.strip()
    return not (name.startswith("H") or name.startswith("D"))


def is_protein_residue(residue):
    """Return True if residue is part of a protein chain."""
    resname = residue.resname.strip()
    return resname in AA_NAMES or resname in MODIFIED_AA_NAMES


def is_ligand_residue(residue):
    """Return True if HETATM residue is a ligand candidate."""
    hetflag = residue.id[0]
    if hetflag == " ":
        return False
    resname = residue.resname.strip()
    if resname in WATER_NAMES:
        return False
    if resname in MODIFIED_AA_NAMES:
        return False
    if resname in AA_NAMES or resname in NUCLEIC_NAMES:
        return False
    if resname in CRYSTALLIZATION_ADDITIVE_NAMES:
        return False
    return True


def get_ligand_residue_id(residue):
    """Return a tuple uniquely identifying a ligand residue."""
    hetflag, resseq, icode = residue.id
    chain_id = residue.get_parent().id if residue.get_parent() else ""
    return (residue.resname.strip(), chain_id, resseq, icode.strip())


def count_ligand_heavy_atoms(residue):
    """Count heavy atoms in a ligand residue."""
    count = 0
    for atom in residue.get_atoms():
        if is_heavy_atom(atom):
            count += 1
    return count


def get_chain_sequence_and_atoms(chain):
    """
    Extract observed protein sequence and per-residue heavy atom coordinates.
    Returns None if chain has no protein residues/atoms.
    """
    seq_chars = []
    residue_list = []
    coord_rows = []
    atom_to_residue_idx = []

    aa_1letter = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    aa_1letter.update(MODIFIED_AA_1LETTER)

    for residue in chain.get_residues():
        if not is_protein_residue(residue):
            continue
        resname = residue.resname.strip()
        if resname not in aa_1letter:
            continue

        start_idx = len(coord_rows)
        has_atom = False
        for atom in residue.get_atoms():
            if not is_heavy_atom(atom):
                continue
            coord = atom.get_coord()
            if coord is None or len(coord) != 3:
                continue
            coord_rows.append(coord)
            atom_to_residue_idx.append(len(residue_list))
            has_atom = True

        if has_atom:
            seq_chars.append(aa_1letter[resname])
            residue_list.append(residue)

    if not coord_rows:
        return None

    return {
        "seq": "".join(seq_chars),
        "residue_list": residue_list,
        "coords": np.array(coord_rows, dtype=float),
        "atom_to_residue_idx": atom_to_residue_idx,
    }


def get_ligand_atoms(ligand_residue):
    """Return heavy atom coordinates for a ligand residue."""
    coords = []
    for atom in ligand_residue.get_atoms():
        if not is_heavy_atom(atom):
            continue
        coord = atom.get_coord()
        if coord is None or len(coord) != 3:
            continue
        coords.append(coord)
    return np.array(coords, dtype=float) if coords else None


def process_structure_file(filepath, distance_threshold, min_ligand_atoms, min_binding_residues):
    """
    Process a single structure file and return binding site records.

    Returns dict with keys:
        records: list of record dicts
        error: str or None
        n_chains: int
        n_ligands: int
    """
    filepath = Path(filepath)
    pdb_id = filepath.stem.split(".")[0].upper()
    suffix = filepath.suffixes

    if ".pdb" in suffix and ".gz" in suffix:
        file_type = "pdb"
        parser = PDBParser(QUIET=True)
        opener = lambda: gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    elif ".cif" in suffix and ".gz" in suffix:
        file_type = "cif"
        parser = MMCIFParser(QUIET=True)
        opener = lambda: gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    elif ".pdb" in suffix:
        file_type = "pdb"
        parser = PDBParser(QUIET=True)
        opener = lambda: open(filepath, "r", encoding="utf-8", errors="replace")
    elif ".cif" in suffix:
        file_type = "cif"
        parser = MMCIFParser(QUIET=True)
        opener = lambda: open(filepath, "r", encoding="utf-8", errors="replace")
    else:
        return {"records": [], "error": f"Unsupported file extension: {filepath.name}", "n_chains": 0, "n_ligands": 0}

    records = []
    try:
        with opener() as fh:
            structure = parser.get_structure(pdb_id, fh)

        models = list(structure.get_models())
        if not models:
            return {"records": [], "error": "No models found", "n_chains": 0, "n_ligands": 0}
        model = models[0]

        # Collect ligands
        ligands = []
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if not is_ligand_residue(residue):
                    continue
                n_heavy = count_ligand_heavy_atoms(residue)
                if n_heavy < min_ligand_atoms:
                    continue
                coords = get_ligand_atoms(residue)
                if coords is None or len(coords) == 0:
                    continue
                ligand_id = get_ligand_residue_id(residue)
                ligands.append((chain.id.strip(), ligand_id, residue, coords, n_heavy))

        n_ligands_total = len(ligands)
        n_protein_chains = 0

        for chain in model.get_chains():
            chain_id = chain.id.strip()
            chain_data = get_chain_sequence_and_atoms(chain)
            if chain_data is None or len(chain_data["seq"]) == 0:
                continue
            n_protein_chains += 1

            protein_coords = chain_data["coords"]
            atom_to_residue_idx = chain_data["atom_to_residue_idx"]
            seq = chain_data["seq"]
            n_residues = len(seq)

            for lig_chain_id, ligand_id, ligand_residue, lig_coords, n_heavy in ligands:
                tree = cKDTree(lig_coords)
                hit_indices = tree.query_ball_point(protein_coords, r=distance_threshold)

                binding_residue_set = set()
                for atom_idx, neighbors in enumerate(hit_indices):
                    if neighbors:
                        binding_residue_set.add(atom_to_residue_idx[atom_idx])

                if len(binding_residue_set) <= min_binding_residues:
                    continue

                binding_mask = [0] * n_residues
                for res_idx in binding_residue_set:
                    binding_mask[res_idx] = 1

                lig_name, lig_chain, lig_resnum, lig_icode = ligand_id
                ligand_id_str = f"{lig_name}_{lig_chain}_{lig_resnum}{lig_icode}"
                unique_id = f"{pdb_id}_{chain_id}_{ligand_id_str}"

                records.append({
                    "unique_id": unique_id,
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "ligand_name": lig_name,
                    "ligand_chain": lig_chain,
                    "ligand_resnum": lig_resnum,
                    "ligand_icode": lig_icode,
                    "ligand_id": ligand_id_str,
                    "file_type": file_type,
                    "aa_seq": seq,
                    "labels": str(binding_mask),
                    "num_binding_residues": len(binding_residue_set),
                    "num_ligand_heavy_atoms": n_heavy,
                    "protein_length": n_residues,
                })

        return {
            "records": records,
            "error": None,
            "n_chains": n_protein_chains,
            "n_ligands": n_ligands_total,
        }

    except Exception as e:
        return {
            "records": [],
            "error": f"{type(e).__name__}: {str(e)}",
            "n_chains": 0,
            "n_ligands": 0,
        }


def process_file_wrapper(args):
    """Wrapper for parallel execution that returns (filepath, result)."""
    filepath, distance_threshold, min_ligand_atoms, min_binding_residues = args
    result = process_structure_file(filepath, distance_threshold, min_ligand_atoms, min_binding_residues)
    return filepath, result
