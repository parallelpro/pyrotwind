"""Shared fixtures for the pyrotwind test suite.

Track/namelist/golden-output fixtures live in tests/fixtures/ and are
generated ahead of time by generate_fixtures.py from a real MESA model
grid -- see that script's docstring. Tests only ever read these small,
committed files; they don't need the external grid at run time.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"

# pyrotwind.cpython-*.so lives directly in the repo root (HERE.parent),
# not inside a further-nested "pyrotwind" package -- so this is a bare
# `import pyrotwind`, matching compile_notes' own documented usage from
# within this repo (as opposed to the outer sg-rotation project, which
# imports this repo as a subdirectory package via `from pyrotwind import
# pyrotwind`).
sys.path.insert(0, str(HERE.parent))
import pyrotwind  # noqa: E402

OUT_NAMES = ["swe", "swc", "swcore", "sj", "sje", "sjc", "sjcore", "sprot", "faccen"]
NML_NAMES = ["nml_solid", "nml_twozone", "nml_threezone"]
TRACK_NAMES = ["dense_m1.0", "dense_m1.5"]


def load_track(name):
    return np.load(FIXTURES / f"track_{name}.npz")


def setup_namelist(nml_name):
    """(Re-)read the given namelist. Must be called before rotwind() --
    it sets lsolid/lthreezone/iwind/... and pdisk/tlock for every track
    processed until the next call."""
    nml = str(FIXTURES / f"{nml_name}.nml")
    pyrotwind.constm.init_const(nml)
    pyrotwind.params.read_wind_params(nml)


def run_track(track_name):
    """Call rotwind() on a fixture track under whatever namelist config
    was last set up via setup_namelist(). Returns rotwind()'s raw tuple."""
    d = load_track(track_name)
    nm = len(d["sage"])
    return pyrotwind.rotwind(
        d["taucouple"], d["taucouplecm"], d["sage"], d["sl"], d["sr"],
        d["steffl"], d["smcz"], d["srcz"], d["sxc"], d["si"], d["sie"],
        d["sic"], d["sicore"], d["staucz"], d["shec"], d["sp"], d["sm"], nm,
    )


def run(nml_name, track_name):
    """setup_namelist(nml_name) then run_track(track_name), for the
    common case of one track under one namelist in isolation."""
    setup_namelist(nml_name)
    return run_track(track_name)


@pytest.fixture(scope="session")
def golden():
    return np.load(FIXTURES / "golden_outputs.npz")
