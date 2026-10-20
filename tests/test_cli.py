from pathlib import Path

import pytest
import tyro

from consiter import cli


@pytest.fixture
def inputs(tmp_path):
    for f in ("ref.fa", "R1.fq", "R2.fq"):
        (tmp_path / f).write_text("x\n")
    return tmp_path


def parse(args):
    return tyro.cli(cli.Config, args=args, config=(tyro.conf.FlagCreatePairsOff,))


def test_short_aliases_and_defaults(inputs):
    cfg = parse(["-r", str(inputs / "ref.fa"), "-n", "S1",
                 "-1", str(inputs / "R1.fq"), "-2", str(inputs / "R2.fq")])  # fmt: skip
    assert cfg == cli.Config(
        ref=inputs / "ref.fa", name="S1", left=inputs / "R1.fq", right=inputs / "R2.fq"
    )
    assert cfg.outdir == Path(".")
    assert (cfg.max_iter, cfg.threads, cfg.ploidy, cfg.min_cov, cfg.xmx) == (10, 16, 1, 20, 50)
    assert (cfg.call_min_cov, cfg.caller, cfg.stop) == (2, "gatk", "both")
    assert not cfg.keep_tmp and not cfg.no_indel and not cfg.no_dedup
    assert cli.validate(cfg) == []


def test_long_options(inputs):
    cfg = parse([
        "--ref", str(inputs / "ref.fa"), "--name", "S1",
        "--left", str(inputs / "R1.fq"), "--right", str(inputs / "R2.fq"),
        "--outdir", "out", "--max-iter", "3", "--threads", "4", "--ploidy", "2",
        "--min-cov", "10", "--call-min-cov", "3", "--caller", "bcftools",
        "--stop", "consensus", "--xmx", "8", "--keep-tmp", "--no-indel", "--no-dedup",
    ])  # fmt: skip
    assert (cfg.max_iter, cfg.threads, cfg.ploidy, cfg.min_cov, cfg.xmx) == (3, 4, 2, 10, 8)
    assert (cfg.call_min_cov, cfg.caller, cfg.stop) == (3, "bcftools", "consensus")
    assert cfg.outdir == Path("out")
    assert cfg.keep_tmp and cfg.no_indel and cfg.no_dedup


def test_required_arguments():
    with pytest.raises(SystemExit):
        parse(["-n", "S1"])


def test_validate_reports_all_problems(tmp_path):
    cfg = cli.Config(
        ref=tmp_path / "missing.fa", name="bad name", left=tmp_path / "a", right=tmp_path / "b",
        max_iter=0, threads=0, min_cov=-1, call_min_cov=-1,
    )  # fmt: skip
    errors = cli.validate(cfg)
    assert any("reference" in e for e in errors)
    assert any("left reads" in e for e in errors)
    assert any("right reads" in e for e in errors)
    assert any("--max-iter" in e for e in errors)
    assert any("--threads" in e for e in errors)
    assert any(e.startswith("--min-cov") for e in errors)
    assert any(e.startswith("--call-min-cov") for e in errors)
    assert any("--name" in e for e in errors)


def test_validate_rejects_a_ploidy_bcftools_cannot_set(inputs):
    kwargs = dict(ref=inputs / "ref.fa", name="S1", left=inputs / "R1.fq",
                  right=inputs / "R2.fq", outdir=inputs / "out", ploidy=3)  # fmt: skip
    assert cli.validate(cli.Config(**kwargs)) == []  # GATK takes any ploidy
    errors = cli.validate(cli.Config(**kwargs, caller="bcftools"))
    assert any("--ploidy 1 or 2" in e for e in errors)


def test_main_checks_for_bcftools_only_when_it_is_the_caller(inputs, monkeypatch, capsys):
    monkeypatch.setattr(cli.tools, "missing_programs", lambda programs: list(programs))
    args = ["-r", str(inputs / "ref.fa"), "-n", "S1",
            "-1", str(inputs / "R1.fq"), "-2", str(inputs / "R2.fq")]  # fmt: skip
    with pytest.raises(SystemExit) as gatk:
        cli.main(args)
    with pytest.raises(SystemExit) as bcf:
        cli.main([*args, "--caller", "bcftools"])
    assert "bcftools" not in str(gatk.value.code)
    assert "bcftools" in str(bcf.value.code)


def test_main_exits_on_invalid_input(tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["-r", str(tmp_path / "nope.fa"), "-n", "S1",
                  "-1", str(tmp_path / "a"), "-2", str(tmp_path / "b")])  # fmt: skip
    assert "Cannot find reference" in str(exc.value.code)


def test_validate_rejects_existing_tmp(inputs):
    (inputs / "tmp").mkdir()
    cfg = cli.Config(ref=inputs / "ref.fa", name="S1", left=inputs / "R1.fq",
                     right=inputs / "R2.fq", outdir=inputs)  # fmt: skip
    assert any("already exists" in e for e in cli.validate(cfg))
