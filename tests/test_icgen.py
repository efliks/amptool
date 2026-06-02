from pathlib import Path

import pytest

from amptool.icgen import process_monomer


@pytest.fixture
def monomer_pdb_path() -> Path:
    """Return the path to the monomer.pdb test file."""
    return Path(__file__).parent / "monomer.pdb"


@pytest.fixture
def expected_ic_path() -> Path:
    """Return the path to the expected monomer_reorder.ic file."""
    return Path(__file__).parent / "monomer_reorder.ic"


def test_process_monomer(
    monomer_pdb_path: Path, expected_ic_path: Path, tmp_path: Path
) -> None:
    """Test that process_monomer generates correct internal coordinates from monomer.pdb."""
    # Create a temporary output path
    output_path = tmp_path / "test_output.ic"

    # Process the monomer
    process_monomer(monomer_pdb_path, output_path)

    # Read the generated output
    with output_path.open() as f:
        generated_content = f.read()

    # Read the expected output
    with expected_ic_path.open() as f:
        expected_content = f.read()

    # Compare the outputs
    assert generated_content == expected_content


def test_process_monomer_file_exists(monomer_pdb_path: Path, tmp_path: Path) -> None:
    """Test that process_monomer creates the output file."""
    output_path = tmp_path / "output.ic"

    # Verify file doesn't exist before
    assert not output_path.exists()

    # Process the monomer
    process_monomer(monomer_pdb_path, output_path)

    # Verify file was created
    assert output_path.exists()
    assert output_path.stat().st_size > 0

    # Uncomment to regenerate baseline
    # result_path = Path(__file__).parent / "monomer_reorder.ic"
    #
    # with output_path.open() as f:
    #     lines = f.readlines()
    #     with result_path.open("w") as f:
    #         f.writelines(lines)


def test_process_monomer_output_format(monomer_pdb_path: Path, tmp_path: Path) -> None:
    """Test that process_monomer output has the correct format."""
    output_path = tmp_path / "format_check.ic"

    process_monomer(monomer_pdb_path, output_path)

    with output_path.open() as f:
        lines = f.readlines()

    # Each line should have 7 space-separated values
    for line in lines:
        parts = line.strip().split()
        assert len(parts) == 7, f"Expected 7 values per line, got {len(parts)}: {line}"

        # First 4 should be atom labels (strings)
        for i in range(4):
            assert parts[i].isalnum(), (
                f"Expected atom label at position {i}, got {parts[i]}"
            )

        # Last 3 should be numeric values (bond, angle, dihedral)
        for i in range(4, 7):
            assert (
                parts[i].replace(".", "").replace("-", "").isdigit() or "." in parts[i]
            ), f"Expected numeric value at position {i}, got {parts[i]}"
