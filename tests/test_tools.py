import os
import sys

import pytest

from consiter import tools

BOWTIE2_SUMMARY = """\
100000 reads; of these:
  100000 (100.00%) were paired; of these:
    3402 (3.40%) aligned concordantly 0 times
    95911 (95.91%) aligned concordantly exactly 1 time
    687 (0.69%) aligned concordantly >1 times
    ----
    3402 pairs aligned concordantly 0 times; of these:
      211 (6.20%) aligned discordantly 1 time
97.53% overall alignment rate
"""


def test_parse_alignment_rate():
    assert tools.parse_alignment_rate(BOWTIE2_SUMMARY) == pytest.approx(97.53)
    assert tools.parse_alignment_rate("100.00% overall alignment rate\n") == 100.0
    assert tools.parse_alignment_rate("0.00% overall alignment rate\n") == 0.0


def test_parse_alignment_rate_missing():
    with pytest.raises(ValueError, match="overall alignment rate"):
        tools.parse_alignment_rate("Error: reads file does not look like a FASTQ file\n")


def test_run_logs_command_and_output(tmp_path):
    log = tmp_path / "logs" / "step.log"
    tools.run([sys.executable, "-c", "print('hello')"], log)
    text = log.read_text()
    assert text.startswith("$ ")
    assert "hello" in text


def test_run_stdout_redirect(tmp_path):
    out = tmp_path / "out.txt"
    tools.run([sys.executable, "-c", "print('data')"], tmp_path / "step.log", stdout=out)
    assert out.read_text() == "data\n"


def test_run_failure_raises_with_log_tail(tmp_path):
    log = tmp_path / "step.log"
    code = "import sys; print('boom', file=sys.stderr); sys.exit(3)"
    with pytest.raises(tools.ToolError, match="status 3") as exc:
        tools.run([sys.executable, "-c", code], log)
    assert "boom" in str(exc.value)
    assert str(log) in str(exc.value)


def test_run_handles_spaces_in_paths(tmp_path):
    out = tmp_path / "dir with space" / "out.txt"
    out.parent.mkdir()
    tools.run([sys.executable, "-c", "print(1)"], out.parent / "a log.log", stdout=out)
    assert out.read_text() == "1\n"


def test_missing_programs():
    assert tools.missing_programs(["definitely-not-a-real-program-xyz"]) == [
        "definitely-not-a-real-program-xyz"
    ]
    assert tools.missing_programs([]) == []


def test_variant_filtration_uses_min_depth_and_no_sor(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(tools, "run", lambda cmd, log, stdout=None: calls.append(cmd))
    tools.variant_filtration(tmp_path / "r.fa", tmp_path / "in.vcf", tmp_path / "out.vcf", 7, 2,
                             tmp_path / "f.log")  # fmt: skip
    (cmd,) = calls
    exprs = [cmd[i + 1] for i, a in enumerate(cmd) if a == "--filter-expression"]
    assert "DP < 7" in exprs
    assert not any("SOR" in e for e in exprs)
    assert any("FS" in e for e in exprs)


def test_bowtie2_map_reports_downstream_failure(monkeypatch, tmp_path):
    # samtools sort dies (e.g. disk full); bowtie2 and samtools view then die of SIGPIPE.
    # The error must name samtools sort and show its message, not just bowtie2's.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "bowtie2").write_text("#!/bin/sh\nyes '@SQ\tSN:x'\n")
    (bin_dir / "samtools").write_text(
        '#!/bin/sh\nif [ "$1" = sort ]; then echo "No space left on device" >&2; exit 1; fi\ncat\n'
    )
    for f in bin_dir.iterdir():
        f.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")

    with pytest.raises(tools.ToolError) as exc:
        tools.bowtie2_map(
            tmp_path / "ref.fa", tmp_path / "R1.fq", tmp_path / "R2.fq", "S1", 1,
            tmp_path / "out.bam", tmp_path / "stats.txt", tmp_path / "map.log",
        )  # fmt: skip
    msg = str(exc.value)
    assert "samtools sort exited with status 1" in msg
    assert "No space left on device" in msg
