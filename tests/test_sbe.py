import types
from pathlib import Path
import suncat.sbe.run as sbe


def test_run_module_builds_the_exact_sbe_command(monkeypatch, tmp_path):
    work = tmp_path / "work"; work.mkdir()
    captured = {}

    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        (work / "x_filt.cnv").write_text("x")   # pretend SBE wrote the output
        return types.SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(sbe.subprocess, "run", fake_run)
    out = sbe.run_module("FilterW.exe", Path("C:/psa/05_filter.psa"), Path("C:/in/x_al.cnv"),
                         "x_filt", sbe_bin=Path("C:/SBE"), work_dir=work, echo=True, warn=True)
    cmd = captured["cmd"]
    # exactly the step01 pattern: exe /p<psa> /i<in> /o<work> /f<stem> /s
    assert cmd[0].endswith("FilterW.exe")
    assert cmd[1].startswith("/p") and cmd[1].endswith("05_filter.psa")
    assert cmd[2].startswith("/i") and cmd[2].endswith("x_al.cnv")
    assert cmd[3].startswith("/o")
    assert cmd[4] == "/fx_filt"
    assert cmd[5] == "/s"
    assert out == work / "x_filt.cnv"


def test_run_module_returns_none_when_sbe_writes_nothing(monkeypatch, tmp_path):
    def fake_run(cmd, **kw):
        return types.SimpleNamespace(returncode=1, stderr="boom")   # failed, no output file
    monkeypatch.setattr(sbe.subprocess, "run", fake_run)
    out = sbe.run_module("LoopEditW.exe", Path("p.psa"), Path("in.cnv"), "stem",
                         sbe_bin=Path("C:/SBE"), work_dir=tmp_path, echo=False, warn=True)
    assert out is None
