"""Ligand resolution shared by the prediction input builders."""

import pytest
from biotite.structure import BondList
from biotite.structure.info import residue

from blocks.networks.prediction.helpers import _resolve_ligand


def _custom_ligand(bonds="ccd"):
    ligand = residue("EOH")
    ligand = ligand[ligand.element != "H"]
    ligand.res_name[:] = "L:0"
    if bonds == "empty":
        ligand.bonds = BondList(len(ligand))
    elif bonds is None:
        ligand.bonds = None
    return ligand


def test_bonded_custom_ligand_resolves_to_its_smiles():
    assert _resolve_ligand("L:0", _custom_ligand()) == ("smiles", "CCO")


@pytest.mark.parametrize("bonds", ["empty", None])
def test_bondless_custom_ligand_needs_specified_smiles(bonds):
    with pytest.raises(ValueError, match="specified_smiles"):
        _resolve_ligand("L:0", _custom_ligand(bonds))
    assert _resolve_ligand(
        "L:0", _custom_ligand(bonds), specified_smiles={"L:0": "CCO"}
    ) == ("smiles", "CCO")


def test_single_atom_without_bonds_still_resolves():
    ion = residue("ZN")
    ion.res_name[:] = "L:1"
    ion.bonds = BondList(1)
    kind, value = _resolve_ligand("L:1", ion)
    assert kind == "smiles" and "Zn" in value
