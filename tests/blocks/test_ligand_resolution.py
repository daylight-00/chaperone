"""Ligand resolution shared by the prediction input builders."""

import pytest
from biotite.structure import BondList
from biotite.structure.info import residue

from blocks.networks.prediction.helpers import _resolve_ligand


def _ethanol():
    ethanol = residue("EOH")
    ethanol = ethanol[ethanol.element != "H"]
    ethanol.res_name[:] = "L:0"
    return ethanol


def test_bonded_custom_ligand_resolves_to_its_smiles():
    assert _resolve_ligand("L:0", _ethanol()) == ("smiles", "CCO")


@pytest.mark.parametrize("bonds", [BondList(3), None])
def test_bondless_custom_ligand_needs_specified_smiles(bonds):
    ethanol = _ethanol()
    ethanol.bonds = bonds
    with pytest.raises(ValueError, match="specified_smiles"):
        _resolve_ligand("L:0", ethanol)
    assert _resolve_ligand("L:0", ethanol, specified_smiles={"L:0": "CCO"}) == (
        "smiles",
        "CCO",
    )


def test_single_atom_without_bonds_still_resolves():
    ion = residue("ZN")
    ion.res_name[:] = "L:1"
    ion.bonds = BondList(1)
    kind, value = _resolve_ligand("L:1", ion)
    assert kind == "smiles" and "Zn" in value
