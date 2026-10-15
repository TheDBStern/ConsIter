import pytest

from consiter.fasta import low_coverage_intervals, read_dict_names, read_fasta_names, rename_headers


def write(path, text):
    path.write_text(text)
    return path


def test_read_dict_names(tmp_path):
    d = write(
        tmp_path / "ref.dict",
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\tM5:abc\n@SQ\tLN:50\tSN:chr2\n",
    )
    assert read_dict_names(d) == ["chr1", "chr2"]


def test_read_fasta_names(tmp_path):
    fa = write(tmp_path / "x.fa", ">seq1 some description\nACGT\n>seq2\nAC\nGT\n")
    assert read_fasta_names(fa) == ["seq1", "seq2"]


def test_rename_single_record(tmp_path):
    # FastaAlternateReferenceMaker writes numbered headers like this
    src = write(tmp_path / "in.fa", ">1 MN908947.3:1-29903\nACGT\nAC\n")
    rename_headers(src, tmp_path / "out.fa", ["MN908947.3"])
    assert (tmp_path / "out.fa").read_text() == ">MN908947.3\nACGT\nAC\n"


def test_rename_with_sample_prefix(tmp_path):
    src = write(tmp_path / "in.fa", ">MN908947.3\nACGT\n")
    rename_headers(src, tmp_path / "out.fa", ["MN908947.3"], prefix="S1")
    assert (tmp_path / "out.fa").read_text() == ">S1|MN908947.3\nACGT\n"


def test_rename_many_records_keeps_order(tmp_path):
    # >1 must not match >10, >11, ... (bug in the original script)
    names = [f"seg{i}" for i in range(1, 13)]
    src = write(tmp_path / "in.fa", "".join(f">{i} x\nAC\n" for i in range(1, 13)))
    rename_headers(src, tmp_path / "out.fa", names, prefix="S")
    assert read_fasta_names(tmp_path / "out.fa") == [f"S|{n}" for n in names]


@pytest.mark.parametrize("names", [["a"], ["a", "b", "c"]])
def test_rename_count_mismatch(tmp_path, names):
    src = write(tmp_path / "in.fa", ">1\nA\n>2\nC\n")
    with pytest.raises(ValueError, match="records"):
        rename_headers(src, tmp_path / "out.fa", names)


def test_low_coverage_intervals():
    bedgraph = ["c1\t0\t10\t0\n", "c1\t10\t20\t19\n", "c1\t20\t30\t20\n", "c2\t0\t5\t3\n"]
    assert list(low_coverage_intervals(bedgraph, 20)) == [
        "c1\t0\t10\n",
        "c1\t10\t20\n",
        "c2\t0\t5\n",
    ]
    assert list(low_coverage_intervals(bedgraph, 0)) == []
