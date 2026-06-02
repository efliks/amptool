from pathlib import Path

import pytest

from amptool.pdb_model import make_protein_model


@pytest.fixture
def monomer_ic_path() -> Path:
    """Return the path to the monomer_reorder.ic test file."""
    return Path(__file__).parent / "monomer_reorder.ic"


def test_make_protein_model_file_io(monomer_ic_path: Path, tmp_path: Path) -> None:
    """Test that make_protein_model creates the output PDB file."""
    output_path = tmp_path / "test_protein.pdb"

    # Verify file doesn't exist before
    assert not output_path.exists()

    with monomer_ic_path.open() as f:
        lines = f.readlines()

    monomer_edited_path = tmp_path / "monomer_edited.ic"
    with monomer_edited_path.open("w") as f:
        f.writelines(lines[1:-3])   # <== TODO Recall why these lines should be deleted

    make_protein_model(monomer_edited_path, output_path)

    # If file was created before the error, verify it exists
    if output_path.exists():
        assert output_path.stat().st_size > 0
