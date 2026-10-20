"""Iteration logic, with every external tool replaced by a fake."""

import csv
import gzip
import logging

import pytest

from consiter import pipeline, tools
from consiter.cli import Config

# A MarkDuplicates metrics file, cut down to the columns that are read.
METRICS = "## METRICS CLASS\nLIBRARY\tREAD_PAIRS_EXAMINED\tPERCENT_DUPLICATION\n1\t1000\t{}\n"


class FakeTools:
    """Records calls and fakes outputs; ``rates`` are returned by successive mappings."""

    def __init__(self, monkeypatch, rates, *, consensus=None, duplicates=0.1, low_bases=7):
        self.rates = list(rates)
        self.consensus = list(consensus or [])
        """Sequence written by successive consensus rounds; ``ACGT`` repeats by default."""
        self.duplicates = duplicates
        self.low_bases = low_bases
        self.mapped_refs = []
        self.called_refs = []
        self.bcftools_calls = []
        self.gatk_gates = []
        self.indexed = []
        self.selected_snps = 0
        self.coverage_depths = []
        self.masked_with = None
        mp = monkeypatch
        mp.setattr(tools, "bowtie2_build", lambda ref, log: None)
        mp.setattr(tools, "bowtie2_map", self.bowtie2_map)
        mp.setattr(tools, "mark_duplicates", self.mark_duplicates)
        mp.setattr(tools, "index_bam", self.index_bam)
        mp.setattr(tools, "coverage_summary", self.coverage_summary)
        mp.setattr(tools, "prepare_reference", lambda ref, log: None)
        mp.setattr(tools, "haplotype_caller", self.haplotype_caller)
        mp.setattr(tools, "variant_filtration", self.variant_filtration)
        mp.setattr(tools, "select_snps", self.select_snps)
        mp.setattr(tools, "alternate_reference", self.alternate_reference)
        mp.setattr(tools, "bcftools_consensus", self.bcftools_consensus)
        mp.setattr(tools, "mask_fasta", self.mask_fasta)

    def bowtie2_map(self, ref, left, right, name, threads, out_bam, stats, log):
        self.mapped_refs.append(ref.name)
        out_bam.write_text(f"aligned to {ref.name}")
        return self.rates.pop(0)

    def mark_duplicates(self, in_bam, out_bam, metrics, xmx, log):
        out_bam.write_text(in_bam.read_text())
        out_bam.with_suffix(".bai").write_text("index")
        metrics.write_text(METRICS.format(self.duplicates))

    def index_bam(self, bam, log):
        self.indexed.append(bam.name)
        bam.with_suffix(".bai").write_text("index")

    def coverage_summary(self, bam, min_cov, out_bed, log):
        self.coverage_depths.append(min_cov)
        out_bed.write_text("segA\t0\t1\n")
        return tools.Coverage(low_bases=self.low_bases, total_bases=8, mean_depth=30.0)

    def haplotype_caller(self, ref, bam, vcf, ploidy, xmx, log):
        self.called_refs.append(ref.name)
        assert bam.read_text() == f"aligned to {ref.name}"
        vcf.touch()

    def variant_filtration(self, ref, vcf, out, min_depth, xmx, log):
        self.gatk_gates.append(min_depth)
        out.touch()

    def select_snps(self, ref, vcf, out, xmx, log):
        self.selected_snps += 1
        out.touch()

    def _sequence(self):
        return self.consensus.pop(0) if self.consensus else "ACGT"

    def alternate_reference(self, ref, vcf, out, xmx, log):
        # mimic GATK's numbered headers
        n = sum(line.startswith(">") for line in ref.read_text().splitlines())
        seq = self._sequence()
        out.write_text("".join(f">{i} x:1-4\n{seq}\n" for i in range(1, n + 1)))

    def bcftools_consensus(self, ref, bam, out, min_depth, ploidy, log, *, snps_only=False):
        self.bcftools_calls.append((ref.name, min_depth, ploidy, snps_only))
        names = [line for line in ref.read_text().splitlines() if line.startswith(">")]
        seq = self._sequence()
        out.write_text("".join(f"{name}\n{seq}\n" for name in names))

    def mask_fasta(self, ref, bed, out_fasta, log):
        self.masked_with = (ref.name, bed.name)
        out_fasta.write_text(ref.read_text())


@pytest.fixture
def cfg(tmp_path):
    (tmp_path / "ref.fa").write_text(">segA desc\nACGT\n>segB\nACGT\n")
    (tmp_path / "R1.fq").touch()
    (tmp_path / "R2.fq").touch()
    return Config(
        ref=tmp_path / "ref.fa", name="S1", left=tmp_path / "R1.fq", right=tmp_path / "R2.fq",
        outdir=tmp_path / "out", max_iter=5,
    )  # fmt: skip


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_stops_when_rate_does_not_improve(monkeypatch, cfg):
    # distinct consensus sequences, so neither the fixed point nor the cycle check fires
    fake = FakeTools(monkeypatch, [90.0, 95.0, 95.0], consensus=["ACGA", "ACGG"])
    out = pipeline.run(cfg)

    assert fake.mapped_refs == ["ref.iter0.fa", "ref.iter1.fa", "ref.iter2.fa"]
    assert fake.called_refs == ["ref.iter0.fa", "ref.iter1.fa"]
    # latest consensus is masked with the reads mapped to it
    assert fake.masked_with == ("ref.iter2.fa", "iter2.lowcov.bed")
    assert out.iterations == 3
    assert out.stop_reason == "rate"
    assert out.alignment_rates == [90.0, 95.0, 95.0]
    assert out.masked_bases == 7
    assert (out.masked_fraction, out.breadth) == (7 / 8, 1 / 8)
    assert out.bam.read_text() == "aligned to ref.iter2.fa"
    assert out.bai.is_file()
    assert out.consensus.read_text() == ">segA\nACGG\n>segB\nACGG\n"
    assert out.consensus_sample_name.read_text() == ">S1|segA\nACGG\n>S1|segB\nACGG\n"
    assert not (cfg.outdir / "tmp").exists()


def test_first_round_always_makes_a_consensus(monkeypatch, cfg):
    fake = FakeTools(monkeypatch, [99.0, 98.0], consensus=["ACGA"])
    pipeline.run(cfg)
    assert fake.called_refs == ["ref.iter0.fa"]
    assert fake.masked_with[0] == "ref.iter1.fa"


def test_stops_at_max_iter(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "max_iter": 2})
    fake = FakeTools(monkeypatch, [80.0, 85.0, 90.0], consensus=["ACGA", "ACGG"])
    out = pipeline.run(cfg)
    assert fake.called_refs == ["ref.iter0.fa", "ref.iter1.fa"]  # two consensus rounds
    assert fake.mapped_refs[-1] == "ref.iter2.fa"  # final consensus is mapped once more
    assert fake.masked_with == ("ref.iter2.fa", "iter2.lowcov.bed")
    assert (out.iterations, out.stop_reason) == (3, "max_iter")


def test_stops_at_a_fixed_point_while_the_rate_is_still_rising(monkeypatch, cfg):
    """The caller applying nothing is the real signal that the loop is finished."""
    fake = FakeTools(monkeypatch, [80.0, 85.0, 90.0], consensus=["ACGA", "ACGA"])
    out = pipeline.run(cfg)
    assert fake.mapped_refs == ["ref.iter0.fa", "ref.iter1.fa"]
    assert (out.iterations, out.stop_reason) == (2, "fixed_point")
    assert out.consensus.read_text() == ">segA\nACGA\n>segB\nACGA\n"


def test_stops_on_a_cycle(monkeypatch, cfg):
    """A caller that oscillates between two consensuses never reaches a fixed point.

    Here round 2 returns the starting reference, which no fixed-point check catches.
    """
    fake = FakeTools(monkeypatch, [80.0, 85.0, 90.0], consensus=["ACGA", "ACGT"])
    out = pipeline.run(cfg)
    assert fake.mapped_refs == ["ref.iter0.fa", "ref.iter1.fa"]
    assert (out.iterations, out.stop_reason) == (2, "cycle")


def test_stop_rate_ignores_a_fixed_point(monkeypatch, cfg):
    """With --stop rate the repeated consensus is remapped, and the flat rate stops it."""
    cfg = Config(**{**cfg.__dict__, "stop": "rate"})
    fake = FakeTools(monkeypatch, [80.0, 85.0, 85.0], consensus=["ACGA", "ACGA"])
    out = pipeline.run(cfg)
    assert fake.mapped_refs == ["ref.iter0.fa", "ref.iter1.fa", "ref.iter2.fa"]
    assert (out.iterations, out.stop_reason) == (3, "rate")


def test_stop_consensus_ignores_a_falling_rate(monkeypatch, cfg):
    """With --stop consensus the loop keeps going while the caller still changes bases."""
    cfg = Config(**{**cfg.__dict__, "stop": "consensus"})
    fake = FakeTools(monkeypatch, [90.0, 80.0, 70.0], consensus=["ACGA", "ACGG", "ACGG"])
    out = pipeline.run(cfg)
    assert fake.called_refs == ["ref.iter0.fa", "ref.iter1.fa", "ref.iter2.fa"]
    assert (out.iterations, out.stop_reason) == (3, "fixed_point")
    assert out.alignment_rates == [90.0, 80.0, 70.0]


def test_the_call_gate_is_separate_from_the_masking_threshold(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "min_cov": 20, "call_min_cov": 2})
    fake = FakeTools(monkeypatch, [90.0, 90.0], consensus=["ACGA"])
    pipeline.run(cfg)
    assert fake.gatk_gates == [2]
    assert fake.coverage_depths == [20, 20]


def test_bcftools_caller(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "caller": "bcftools", "no_indel": True})
    fake = FakeTools(monkeypatch, [90.0, 95.0, 95.0], consensus=["ACGA", "ACGG"])
    out = pipeline.run(cfg)
    assert fake.bcftools_calls == [("ref.iter0.fa", 2, 1, True), ("ref.iter1.fa", 2, 1, True)]
    assert fake.called_refs == []  # GATK not used at all
    assert out.consensus.read_text() == ">segA\nACGG\n>segB\nACGG\n"


def test_no_dedup_keeps_duplicates_and_still_indexes(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "no_dedup": True})
    fake = FakeTools(monkeypatch, [90.0, 90.0], consensus=["ACGA"])

    def fail(*args, **kwargs):
        raise AssertionError("MarkDuplicates must not run with --no-dedup")

    monkeypatch.setattr(tools, "mark_duplicates", fail)
    out = pipeline.run(cfg)
    assert fake.indexed == ["iter0.rmdup.bam", "iter1.rmdup.bam"]
    assert out.bai.is_file()
    assert read_tsv(out.iterations_tsv)[0]["duplicate_fraction"] == ""


def test_warns_about_an_implausible_duplicate_fraction(monkeypatch, cfg, caplog):
    FakeTools(monkeypatch, [90.0, 90.0], duplicates=0.991)
    with caplog.at_level(logging.WARNING, logger="consiter"):
        out = pipeline.run(cfg)
    assert "99.1% of reads were removed as duplicates" in caplog.text
    assert "--no-dedup" in caplog.text
    assert read_tsv(out.iterations_tsv)[0]["duplicate_fraction"] == "0.9910"


def test_warns_when_the_reference_has_no_foothold(monkeypatch, cfg, caplog):
    FakeTools(monkeypatch, [90.0, 90.0], low_bases=8)  # breadth 0.0
    with caplog.at_level(logging.WARNING, logger="consiter"):
        pipeline.run(cfg)
    assert "Only 0.0% of the reference reaches 20x" in caplog.text
    assert "100.0% of the consensus is masked" in caplog.text


def test_no_warning_when_coverage_is_fine(monkeypatch, cfg, caplog):
    FakeTools(monkeypatch, [90.0, 90.0], low_bases=0)
    with caplog.at_level(logging.WARNING, logger="consiter"):
        pipeline.run(cfg)
    assert caplog.text == ""


def test_iterations_tsv(monkeypatch, cfg):
    fake = FakeTools(monkeypatch, [90.0, 95.0, 95.0], consensus=["ACGA", "ACGG"])
    out = pipeline.run(cfg)
    rows = read_tsv(out.iterations_tsv)
    assert [r["iteration"] for r in rows] == ["0", "1", "2"]
    assert [r["alignment_rate"] for r in rows] == ["90.00", "95.00", "95.00"]
    assert [r["changed"] for r in rows] == ["1", "1", ""]  # no round after the last iteration
    assert [r["stop_reason"] for r in rows] == ["", "", "rate"]
    assert rows[0]["breadth"] == "0.1250"
    assert rows[0]["masked_bases"] == "7"
    assert rows[0]["duplicate_fraction"] == "0.1000"
    assert rows[0]["ref_hash"] != rows[1]["ref_hash"]
    assert fake.mapped_refs  # the fakes really ran


def test_iterations_tsv_survives_a_failure(monkeypatch, cfg):
    FakeTools(monkeypatch, [90.0])
    monkeypatch.setattr(tools, "haplotype_caller", lambda *a: 1 / 0)
    with pytest.raises(ZeroDivisionError):
        pipeline.run(cfg)
    rows = read_tsv(cfg.outdir / "S1.iterations.tsv")
    assert [r["iteration"] for r in rows] == ["0"]


def test_no_indel_selects_snps(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "no_indel": True})
    fake = FakeTools(monkeypatch, [90.0, 90.0])
    pipeline.run(cfg)
    assert fake.selected_snps == 1


def test_keep_tmp(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "keep_tmp": True})
    FakeTools(monkeypatch, [90.0, 90.0])
    pipeline.run(cfg)
    assert (cfg.outdir / "tmp" / "ref.iter1.fa").is_file()


def test_gzipped_reference(monkeypatch, cfg):
    with gzip.open(cfg.ref.with_suffix(".fa.gz"), "wt") as fh:
        fh.write(">segA\nACGT\n>segB\nACGT\n")
    cfg = Config(**{**cfg.__dict__, "ref": cfg.ref.with_suffix(".fa.gz")})
    FakeTools(monkeypatch, [90.0, 90.0])
    out = pipeline.run(cfg)
    assert out.consensus.read_text().startswith(">segA\n")


def test_tool_failure_keeps_tmp(monkeypatch, cfg):
    FakeTools(monkeypatch, [90.0])

    def fail(*args):
        raise tools.ToolError("HaplotypeCaller exited with status 1")

    monkeypatch.setattr(tools, "haplotype_caller", fail)
    with pytest.raises(tools.ToolError):
        pipeline.run(cfg)
    assert (cfg.outdir / "tmp").is_dir()


def test_refuses_existing_tmp(monkeypatch, cfg):
    FakeTools(monkeypatch, [90.0, 90.0])
    unrelated = cfg.outdir / "tmp" / "precious.txt"
    unrelated.parent.mkdir(parents=True)
    unrelated.write_text("keep me")
    with pytest.raises(ValueError, match="already exists"):
        pipeline.run(cfg)
    assert unrelated.read_text() == "keep me"
