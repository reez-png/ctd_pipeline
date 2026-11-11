import suncat
from suncat import config


def test_version():
    assert isinstance(suncat.__version__, str) and suncat.__version__


def test_config_paths(tmp_path):
    c = config.for_cruise("P45_06", root=tmp_path)
    assert c.cruise_dir == tmp_path / "cruises" / "P45_06"
    assert c.align_dir.name == "Align_CTD"
    assert c.L3.parts[-2:] == ("L3", "CTD")


def test_root_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("SUNCAT_CTD_ROOT", str(tmp_path))
    assert config.for_cruise("P45_06").ctd_root == tmp_path.resolve()
