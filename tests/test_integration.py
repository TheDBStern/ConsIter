"""End-to-end runs of the ``consiter`` command on simulated data.

Run with ``pytest -m integration``; needs the tools from environment.yml.
"""

from __future__ import annotations

import random
import subprocess

import pytest

from consiter.tools import missing_programs
from simulate import mutate, random_seq, simulate_reads, write_fasta

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(bool(missing_programs()), reason="external tools not on PATH"),
]


def read_fasta(path):
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = []
        else:
            records[name].append(line.strip())
    return {k: "".join(v) for k, v in records.items()}


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    rng = random.Random(42)
    d = tmp_path_factory.mktemp("data")
    ref = {"contigA": random_seq(rng, 3000), "contigB": random_seq(rng, 2000)}
    truth = {
        "contigA": mutate(rng, ref["contigA"], 20, [(1000, "GTA"), (2000, 4)]),
        "contigB": mutate(rng, ref["contigB"], 10, []),
    }
    write_fasta(d / "ref.fa", ref)
    simulate_reads(rng, truth, d / "R1.fq", d / "R2.fq")
    return d, ref, truth


def run_consiter(d, out, *extra):
    cmd = [
        "consiter", "-r", d / "ref.fa", "-n", "S1", "-1", d / "R1.fq", "-2", d / "R2.fq",
        "-o", out, "-t", "2", "--xmx", "2", "-c", "5", *extra,
    ]  # fmt: skip
    subprocess.run([str(c) for c in cmd], check=True)


def bam_lengths(bam):
    header = subprocess.run(
        ["samtools", "view", "-H", str(bam)], check=True, capture_output=True, text=True
    ).stdout
    lengths = {}
    for line in header.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            lengths[fields["SN"]] = int(fields["LN"])
    return lengths


def matches_outside_mask(consensus, truth):
    return len(consensus) == len(truth) and all(
        c == t for c, t in zip(consensus, truth, strict=True) if c != "N"
    )


def test_recovers_snps_and_indels(dataset, tmp_path):
    d, ref, truth = dataset
    out = tmp_path / "out"
    run_consiter(d, out)

    consensus = read_fasta(out / "S1.consensus.fa")
    assert list(consensus) == ["contigA", "contigB"]
    for name in truth:
        assert matches_outside_mask(consensus[name], truth[name]), name
        # ends are expected to be masked, but most of the genome should not be
        assert consensus[name].count("N") < 0.05 * len(truth[name])

    assert list(read_fasta(out / "S1.consensus.sample_name.fa")) == ["S1|contigA", "S1|contigB"]
    assert bam_lengths(out / "S1.bt2.rmdup.bam") == {k: len(v) for k, v in consensus.items()}
    assert (out / "S1.bt2.rmdup.bai").is_file()
    assert not (out / "tmp").exists()
    assert not (d / "ref.fa.fai").exists()  # no index files next to the input


def test_no_indel_keeps_coordinates(dataset, tmp_path):
    d, ref, truth = dataset
    out = tmp_path / "out"
    run_consiter(d, out, "--no-indel", "--keep-tmp")

    consensus = read_fasta(out / "S1.consensus.fa")
    assert {k: len(v) for k, v in consensus.items()} == {k: len(v) for k, v in ref.items()}
    # contigB has only SNPs, so it must be recovered exactly
    assert matches_outside_mask(consensus["contigB"], truth["contigB"])
    assert (out / "tmp" / "logs").is_dir()


def test_tool_failure_is_reported(dataset, tmp_path):
    d, _, _ = dataset
    bad = tmp_path / "bad_R2.fq"
    bad.write_text("not a fastq\n")
    cmd = ["consiter", "-r", d / "ref.fa", "-n", "S1", "-1", d / "R1.fq", "-2", bad,
           "-o", tmp_path / "out", "-t", "2", "--xmx", "2"]  # fmt: skip
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "bowtie2" in proc.stderr
    assert (tmp_path / "out" / "tmp").is_dir()
