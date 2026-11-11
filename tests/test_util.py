import pandas as pd
from suncat.util import pick


def test_pick_exact_then_substring():
    df = pd.DataFrame({"prdM": [1.0], "tv290C": [2.0]})
    assert pick(df, "prdM", "prDM") == "prdM"       # exact match wins
    assert pick(df, "nope", "290") == "tv290C"      # case-insensitive substring
    assert pick(df, "xxx") is None


def test_pick_accepts_tuple_or_varargs():
    df = pd.DataFrame({"t090C": [1.0]})
    assert pick(df, ("tv290C", "t090C")) == "t090C"   # single tuple (step00a form)
    assert pick(df, "tv290C", "t090C") == "t090C"     # varargs (step00b/step00c form)


def test_pick_none_df():
    assert pick(None, "a") is None
