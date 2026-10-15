"""Iteration logic, with every external tool replaced by a fake."""

import gzip

import pytest

from consiter import pipeline, tools
from consiter.cli import Config


class FakeTools:
    """Records calls and fakes outputs; ``rates`` are returned by successive mappings."""

    def __init__(self, monkeypatch, rates):
        self.rates = list(rates)
        self.mapped_refs = []
        self.called_refs = []
        self.selected_snps = 0
        self.masked_with = None
        mp = monkeypatch
        mp.setattr(tools, "bowtie2_build", lambda ref, log: None)
        mp.setattr(tools, "bowtie2_map", self.bowtie2_map)
        mp.setattr(tools, "mark_duplicates", self.mark_duplicates)
        mp.setattr(tools, "prepare_reference", lambda ref, log: None)
        mp.setattr(tools, "haplotype_caller", self.haplotype_caller)
        mp.setattr(
            tools, "variant_filtration", lambda ref, vcf, out, min_depth, xmx, log: out.touch()
        )
        mp.setattr(tools, "select_snps", self.select_snps)
        mp.setattr(tools, "alternate_reference", self.alternate_reference)
        mp.setattr(tools, "mask_low_coverage", self.mask_low_coverage)

    def bowtie2_map(self, ref, left, right, name, threads, out_bam, stats, log):
        self.mapped_refs.append(ref.name)
        out_bam.write_text(f"aligned to {ref.name}")
        return self.rates.pop(0)

    def mark_duplicates(self, in_bam, out_bam, metrics, xmx, log):
        out_bam.write_text(in_bam.read_text())
        out_bam.with_suffix(".bai").write_text("index")

    def haplotype_caller(self, ref, bam, vcf, ploidy, xmx, log):
        self.called_refs.append(ref.name)
        assert bam.read_text() == f"aligned to {ref.name}"
        vcf.touch()

    def select_snps(self, ref, vcf, out, xmx, log):
        self.selected_snps += 1
        out.touch()

    def alternate_reference(self, ref, vcf, out, xmx, log):
        # mimic GATK's numbered headers
        n = sum(line.startswith(">") for line in ref.read_text().splitlines())
        out.write_text("".join(f">{i} x:1-4\nACGT\n" for i in range(1, n + 1)))

    def mask_low_coverage(self, ref, bam, min_cov, out_fasta, out_bed, log):
        self.masked_with = (ref.name, bam.read_text())
        out_fasta.write_text(ref.read_text())
        return 7


@pytest.fixture
def cfg(tmp_path):
    (tmp_path / "ref.fa").write_text(">segA desc\nACGT\n>segB\nACGT\n")
    (tmp_path / "R1.fq").touch()
    (tmp_path / "R2.fq").touch()
    return Config(
        ref=tmp_path / "ref.fa", name="S1", left=tmp_path / "R1.fq", right=tmp_path / "R2.fq",
        outdir=tmp_path / "out", max_iter=5,
    )  # fmt: skip


def test_stops_when_rate_does_not_improve(monkeypatch, cfg):
    fake = FakeTools(monkeypatch, [90.0, 95.0, 95.0])
    out = pipeline.run(cfg)

    assert fake.mapped_refs == ["ref.iter0.fa", "ref.iter1.fa", "ref.iter2.fa"]
    assert fake.called_refs == ["ref.iter0.fa", "ref.iter1.fa"]
    # latest consensus is masked with the reads mapped to it
    assert fake.masked_with == ("ref.iter2.fa", "aligned to ref.iter2.fa")
    assert out.iterations == 3
    assert out.alignment_rates == [90.0, 95.0, 95.0]
    assert out.masked_bases == 7
    assert out.bam.read_text() == "aligned to ref.iter2.fa"
    assert out.bai.is_file()
    assert out.consensus.read_text() == ">segA\nACGT\n>segB\nACGT\n"
    assert out.consensus_sample_name.read_text() == ">S1|segA\nACGT\n>S1|segB\nACGT\n"
    assert not (cfg.outdir / "tmp").exists()


def test_first_round_always_makes_a_consensus(monkeypatch, cfg):
    fake = FakeTools(monkeypatch, [99.0, 98.0])
    pipeline.run(cfg)
    assert fake.called_refs == ["ref.iter0.fa"]
    assert fake.masked_with[0] == "ref.iter1.fa"


def test_stops_at_max_iter(monkeypatch, cfg):
    cfg = Config(**{**cfg.__dict__, "max_iter": 2})
    fake = FakeTools(monkeypatch, [80.0, 85.0, 90.0])
    out = pipeline.run(cfg)
    assert fake.called_refs == ["ref.iter0.fa", "ref.iter1.fa"]  # two consensus rounds
    assert fake.mapped_refs[-1] == "ref.iter2.fa"  # final consensus is mapped once more
    assert fake.masked_with == ("ref.iter2.fa", "aligned to ref.iter2.fa")
    assert out.iterations == 3


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
