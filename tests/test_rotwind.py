"""Regression and sanity tests for pyrotwind.rotwind(), driven by real
MESA tracks (see generate_fixtures.py / tests/fixtures/).
"""
import numpy as np
import pytest

from conftest import OUT_NAMES, NML_NAMES, TRACK_NAMES, load_track, run, run_track, setup_namelist, pyrotwind


def test_import_and_init():
    """constm.init_const / params.read_wind_params should run cleanly and
    populate the expected module state."""
    setup_namelist("nml_threezone")
    assert pyrotwind.constm.solm == pytest.approx(1.9891e33)
    assert pyrotwind.params.lthreezone == True  # noqa: E712 (Fortran logical, not a Python bool)
    assert pyrotwind.params.lsolid == False  # noqa: E712
    assert pyrotwind.params.pdisk == pytest.approx(4.0)


@pytest.mark.parametrize("track_name", TRACK_NAMES)
@pytest.mark.parametrize("nml_name", NML_NAMES)
def test_rotwind_matches_golden(nml_name, track_name, golden):
    """rotwind() output for each (namelist, track) fixture combination
    must match the recorded golden reference exactly. This is the main
    regression guard -- any change to the physics or numerics, intended
    or not, shows up here. Each case is set up from a fresh
    read_wind_params call, independent of test execution order."""
    out = run(nml_name, track_name)
    key = f"{nml_name}__{track_name}"
    for out_name, arr in zip(OUT_NAMES, out[:-1]):
        expected = golden[f"{key}__{out_name}"]
        np.testing.assert_array_equal(
            arr, expected, err_msg=f"{key}__{out_name} differs from golden reference"
        )
    assert int(out[-1]) == int(golden[f"{key}__iermsg"])


def test_no_state_leak_between_tracks(golden):
    """Regression guard for a real bug: the rossby-cutoff flag (lrocrit)
    used to live in the params module and was only reset once per *run*
    (in read_wind_params), not once per *track* -- so a track that
    crossed the critical rossby number could leak that state into
    whichever track ran next in the same session. It's now local to
    solidevol/drevol, threaded into int1zone as an explicit argument.

    Running dense_m1.0 immediately before dense_m1.5 in the same
    read_wind_params session (mirroring the `for track in tracks:
    rotwind(...)` loop real callers use) must give dense_m1.5 the exact
    same result as running it in isolation -- which is how the golden
    reference itself was generated, so comparing against it here is
    exactly the regression check for this bug class.
    """
    setup_namelist("nml_solid")
    run_track("dense_m1.0")  # exercises int1zone's rossby-cutoff branch
    out = run_track("dense_m1.5")
    key = "nml_solid__dense_m1.5"
    for out_name, arr in zip(OUT_NAMES, out[:-1]):
        np.testing.assert_array_equal(arr, golden[f"{key}__{out_name}"])
    assert int(out[-1]) == int(golden[f"{key}__iermsg"])


def test_deterministic_repeat():
    """Calling rotwind() twice with identical inputs in the same process
    must give identical output -- a general regression guard against any
    future state accidentally persisting across calls."""
    out1 = run("nml_twozone", "dense_m1.0")
    out2 = run_track("dense_m1.0")
    for a1, a2 in zip(out1[:-1], out2[:-1]):
        np.testing.assert_array_equal(a1, a2)
    assert out1[-1] == out2[-1]


def test_disk_locked_initial_period():
    """Every track should start disk-locked at pdisk (days), regardless
    of physics model."""
    for nml_name in NML_NAMES:
        out = run(nml_name, "dense_m1.0")
        sprot = out[7]
        assert sprot[0] == pytest.approx(4.0), nml_name


def test_threezone_core_rotates_independently():
    """Sanity check that the three-zone path actually exercises distinct
    core dynamics on this fixture, not silently falling back to the
    two-zone/solid behavior (swcore corotating with swc) throughout --
    i.e. that this fixture is a meaningful test of int3zone, not just of
    threezoneevol's fully-convective/no-core-yet fallback paths."""
    out = run("nml_threezone", "dense_m1.0")
    swc, swcore = out[1], out[2]
    d = load_track("dense_m1.0")
    has_core = d["sicore"] > 0
    assert np.any(has_core)
    assert not np.allclose(swc[has_core], swcore[has_core])
