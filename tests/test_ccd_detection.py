"""Portable CCD detection for ligand-containing predictor inputs."""

import pytest
from biotite.structure import Atom, array
from biotite.structure.info import residue

from blocks.networks.prediction.helpers import _resolve_ligand
from blocks.utils.oracles import is_ccd_name


def test_bundled_ccd_without_site_mirror(monkeypatch):
    monkeypatch.setattr("chaperone.get_setting_field", lambda name, default=None: None)
    assert is_ccd_name("MN")
    assert is_ccd_name("ALA")
    assert not is_ccd_name("NOT_A_COMPONENT")
    assert not is_ccd_name("")
    assert not is_ccd_name(None)


def test_mirror_extensions_and_bundled_fallback(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "chaperone.get_setting_field", lambda name, default=None: str(tmp_path)
    )
    entry = tmp_path / "N" / "NEWCCD" / "NEWCCD.cif"
    entry.parent.mkdir(parents=True)
    entry.write_text("data_NEWCCD\n")
    assert is_ccd_name("NEWCCD")
    assert is_ccd_name("MN")
    assert not is_ccd_name("UNKNOWN_CCD")


@pytest.fixture
def bundled_ccd(monkeypatch):
    monkeypatch.setattr("chaperone.get_setting_field", lambda name, default=None: None)


def test_ccd_ligand_with_unresolved_atoms_keeps_its_code(bundled_ccd):
    heme = residue("HEM")
    heme = heme[heme.element != "H"][:-4]
    assert _resolve_ligand("HEM", heme) == ("ccd", "HEM")


def test_ligand_named_like_a_ccd_code_with_foreign_atoms_is_refused(bundled_ccd):
    ethanol = residue("EOH")
    ethanol.res_name[:] = "LIG"
    with pytest.raises(ValueError, match="specified_smiles"):
        _resolve_ligand("LIG", ethanol)
    assert _resolve_ligand("LIG", ethanol, specified_smiles={"LIG": "CCO"}) == (
        "smiles",
        "CCO",
    )


def test_rfd3_origin_marker_is_not_taken_as_ccd_ori(bundled_ccd):
    origin = array(
        [Atom([0.0, 0.0, 0.0], res_name="ORI", atom_name="ORI", element="X")]
    )
    with pytest.raises(ValueError, match="'ORI'"):
        _resolve_ligand("ORI", origin)
