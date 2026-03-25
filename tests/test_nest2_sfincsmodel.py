# -*- coding: utf-8 -*-
"""
Tests for nest2 functions that use a hydromt_sfincs SfincsModel as the detail model.

Covered combinations
--------------------
* nest2_sfincsmodel_in_sfincsmodel  — overall=SfincsModel,  detail=SfincsModel
* nest2_delft3dfm_in_sfincsmodel    — overall=Delft3DFM,    detail=SfincsModel
* nest2_hurrywave_in_sfincsmodel    — overall=HurryWave,    detail=SfincsModel
* nest2_beware_in_sfincsmodel       — overall=BEWARE,       detail=SfincsModel

Each combination is tested for:
  - Basic "happy path" (BCs are set, no exception raised)
  - Optional parameter variants (bc_path, boundary_water_level_correction, etc.)
  - Combination-specific variants (filter_incoming, return_maximum)
  - Error path: station not found / no stations in bounding box

After calling set_timeseries(), boundary conditions are stored in the component's
``.data`` attribute as an ``xr.Dataset``.  Water-level BCs use variable ``"bzs"``
(dims ``["time", "index"]``); SnapWave BCs use variables ``["hs", "tp", "wd", "ds"]``.

Fixtures are provided by ``conftest.py``.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from conftest import (
    DETAIL_POINT_NAMES,
    N_POINTS,
    N_TIMES,
    OBS_PREFIX,
    make_mock_beware,
    make_mock_delft3dfm,
    make_mock_hurrywave,
    make_mock_overall_sfincs,
)

from cht_nesting.nest2.nest2_beware_in_sfincsmodel import nest2_beware_in_sfincsmodel
from cht_nesting.nest2.nest2_delft3dfm_in_sfincsmodel import (
    nest2_delft3dfm_in_sfincsmodel,
)
from cht_nesting.nest2.nest2_hurrywave_in_sfincsmodel import (
    nest2_hurrywave_in_sfincsmodel,
)
from cht_nesting.nest2.nest2_sfincsmodel_in_sfincsmodel import (
    nest2_sfincsmodel_in_sfincsmodel,
)

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------


def _assert_wl_timeseries_set(
    detail_model, n_times: int = N_TIMES, n_points: int = N_POINTS
) -> None:
    """Assert that water level timeseries (variable ``bzs``) have been set."""
    data = detail_model.water_level.data
    assert data is not None, (
        "water_level.data should not be None after set_timeseries()"
    )
    assert "bzs" in data, (
        f"Expected variable 'bzs' in water_level.data; got {list(data)}"
    )
    shape = data["bzs"].shape  # (time, index)
    assert shape == (n_times, n_points), (
        f"Expected shape ({n_times}, {n_points}), got {shape}"
    )
    assert isinstance(data.indexes["time"], pd.DatetimeIndex), (
        "time coordinate must be a DatetimeIndex"
    )


def _assert_sw_timeseries_set(
    detail_model, n_times: int = N_TIMES, n_points: int = N_POINTS
) -> None:
    """Assert that SnapWave timeseries (variables hs/tp/wd/ds) have been set."""
    data = detail_model.snapwave_boundary_conditions.data
    assert data is not None, (
        "snapwave_boundary_conditions.data should not be None after set_timeseries()"
    )
    for var in ("hs", "tp", "wd", "ds"):
        assert var in data, (
            f"Expected variable '{var}' in snapwave_boundary_conditions.data"
        )
    shape = data["hs"].shape  # (time, index)
    assert shape == (n_times, n_points), (
        f"Expected shape ({n_times}, {n_points}), got {shape}"
    )


# ---------------------------------------------------------------------------
# nest2_sfincsmodel_in_sfincsmodel
# ---------------------------------------------------------------------------


class TestNest2SfincsmodelInSfincsmodel:
    """Tests for overall=SfincsModel, detail=SfincsModel."""

    def test_basic(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """Basic nesting sets water level BCs without error."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        result = nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(sfincs_his_path.parent),
        )
        # Default return value is the water_level component
        assert result is detail_model_wl.water_level
        _assert_wl_timeseries_set(detail_model_wl)

    def test_correction(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
        wl_values: np.ndarray,
    ) -> None:
        """Boundary water level correction is applied uniformly to all time series."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        correction = 0.25
        nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(sfincs_his_path.parent),
            boundary_water_level_correction=correction,
        )
        bzs = detail_model_wl.water_level.data["bzs"].values  # (time, index)
        np.testing.assert_allclose(bzs[:, 0], wl_values[:, 0] + correction, rtol=1e-5)

    def test_return_maximum(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """return_maximum=True returns a pd.Series of peak water levels."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        result = nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(sfincs_his_path.parent),
            return_maximum=True,
        )
        assert isinstance(result, pd.Series), (
            "return_maximum=True should return pd.Series"
        )
        assert len(result) == N_TIMES

    def test_filter_incoming(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """filter_incoming=True runs without error (u/v=0 → no distortion)."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(sfincs_his_path.parent),
            filter_incoming=True,
        )
        _assert_wl_timeseries_set(detail_model_wl)

    def test_bc_path_triggers_write(
        self,
        tmp_path: Path,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """Passing bc_path calls detail.water_level.write() without error."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(sfincs_his_path.parent),
            bc_path=str(tmp_path),
        )
        _assert_wl_timeseries_set(detail_model_wl)

    def test_station_not_found_raises(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """ValueError is raised when an observation point is missing from the his-file."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        with pytest.raises(ValueError, match="not found"):
            nest2_sfincsmodel_in_sfincsmodel(
                overall,
                detail_model_wl,
                obs_point_prefix="wrong_prefix",
                output_path=str(sfincs_his_path.parent),
            )

    def test_default_output_path_from_root(
        self,
        sfincs_his_path: Path,
        detail_model_wl,
    ) -> None:
        """output_path is inferred from overall.root.path when not supplied."""
        overall = make_mock_overall_sfincs(sfincs_his_path.parent)
        # overall.root.path must be a Path (or str) pointing at the output dir
        overall.root.path = sfincs_his_path.parent
        nest2_sfincsmodel_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
        )
        _assert_wl_timeseries_set(detail_model_wl)


# ---------------------------------------------------------------------------
# nest2_delft3dfm_in_sfincsmodel
# ---------------------------------------------------------------------------


class TestNest2SfincsmodelInDelft3dfm:
    """Tests for overall=Delft3DFM, detail=SfincsModel (SfincsModel nested in Delft3DFM)."""

    def _make_bzs(self, sim_times: pd.DatetimeIndex) -> pd.DataFrame:
        """Synthetic water level DataFrame mimicking read_timeseries_output output."""
        rng = np.random.default_rng(7)
        data = rng.uniform(-0.5, 2.0, size=(len(sim_times), N_POINTS))
        return pd.DataFrame(data, index=sim_times)

    def test_basic(
        self,
        tmp_path: Path,
        sim_times: pd.DatetimeIndex,
        detail_model_wl,
    ) -> None:
        """Basic nesting sets water level BCs without error."""
        bzs = self._make_bzs(sim_times)
        overall = make_mock_delft3dfm(tmp_path, bzs)
        nest2_delft3dfm_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(tmp_path),
        )
        _assert_wl_timeseries_set(detail_model_wl)

    def test_correction(
        self,
        tmp_path: Path,
        sim_times: pd.DatetimeIndex,
        detail_model_wl,
    ) -> None:
        """Boundary water level correction is applied uniformly to all time series."""
        bzs = self._make_bzs(sim_times)
        overall = make_mock_delft3dfm(tmp_path, bzs)
        correction = -0.10
        nest2_delft3dfm_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(tmp_path),
            boundary_water_level_correction=correction,
        )
        result = detail_model_wl.water_level.data["bzs"].values  # (time, index)
        np.testing.assert_allclose(
            result[:, 0], bzs.iloc[:, 0].values + correction, rtol=1e-5
        )

    def test_bc_path_triggers_write(
        self,
        tmp_path: Path,
        sim_times: pd.DatetimeIndex,
        detail_model_wl,
    ) -> None:
        """Passing bc_path calls detail.water_level.write() without error."""
        bzs = self._make_bzs(sim_times)
        overall = make_mock_delft3dfm(tmp_path, bzs)
        nest2_delft3dfm_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(tmp_path),
            bc_path=str(tmp_path),
        )
        _assert_wl_timeseries_set(detail_model_wl)

    def test_read_timeseries_called_with_correct_names(
        self,
        tmp_path: Path,
        sim_times: pd.DatetimeIndex,
        detail_model_wl,
    ) -> None:
        """read_timeseries_output is called with the expected station name list."""
        bzs = self._make_bzs(sim_times)
        overall = make_mock_delft3dfm(tmp_path, bzs)
        nest2_delft3dfm_in_sfincsmodel(
            overall,
            detail_model_wl,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(tmp_path),
        )
        call_kwargs = overall.read_timeseries_output.call_args
        # Support both positional and keyword argument call styles.
        name_list = (
            call_kwargs.kwargs.get("name_list")
            if call_kwargs.kwargs.get("name_list") is not None
            else call_kwargs.args[0]
        )
        expected = [f"{OBS_PREFIX}_{n}" for n in DETAIL_POINT_NAMES]
        assert name_list == expected


# ---------------------------------------------------------------------------
# nest2_hurrywave_in_sfincsmodel
# ---------------------------------------------------------------------------


class TestNest2SfincsmodelInHurrywave:
    """Tests for overall=HurryWave, detail=SfincsModel (SfincsModel nested in HurryWave)."""

    def test_basic(
        self,
        hurrywave_his_path: Path,
        detail_model_sw,
    ) -> None:
        """Basic nesting sets SnapWave BCs without error."""
        overall = make_mock_hurrywave(hurrywave_his_path.parent)
        nest2_hurrywave_in_sfincsmodel(
            overall,
            detail_model_sw,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(hurrywave_his_path.parent),
        )
        _assert_sw_timeseries_set(detail_model_sw)

    def test_wave_values_match_fixture(
        self,
        hurrywave_his_path: Path,
        detail_model_sw,
    ) -> None:
        """Wave parameter values in the BCs match the fixture constants."""
        overall = make_mock_hurrywave(hurrywave_his_path.parent)
        nest2_hurrywave_in_sfincsmodel(
            overall,
            detail_model_sw,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(hurrywave_his_path.parent),
        )
        data = detail_model_sw.snapwave_boundary_conditions.data
        np.testing.assert_allclose(data["hs"].values, 1.5, rtol=1e-5)
        np.testing.assert_allclose(data["tp"].values, 8.0, rtol=1e-5)
        np.testing.assert_allclose(data["wd"].values, 270.0, rtol=1e-5)
        np.testing.assert_allclose(data["ds"].values, 20.0, rtol=1e-5)

    def test_default_output_path(
        self,
        hurrywave_his_path: Path,
        detail_model_sw,
    ) -> None:
        """output_path defaults to overall.path."""
        overall = make_mock_hurrywave(hurrywave_his_path.parent)
        overall.path = str(hurrywave_his_path.parent)
        nest2_hurrywave_in_sfincsmodel(
            overall,
            detail_model_sw,
            obs_point_prefix=OBS_PREFIX,
        )
        _assert_sw_timeseries_set(detail_model_sw)

    def test_station_not_found_raises(
        self,
        hurrywave_his_path: Path,
        detail_model_sw,
    ) -> None:
        """ValueError is raised when an observation point is missing from the his-file."""
        overall = make_mock_hurrywave(hurrywave_his_path.parent)
        with pytest.raises(ValueError, match="not found"):
            nest2_hurrywave_in_sfincsmodel(
                overall,
                detail_model_sw,
                obs_point_prefix="wrong_prefix",
                output_path=str(hurrywave_his_path.parent),
            )

    def test_bc_path_triggers_write(
        self,
        tmp_path: Path,
        hurrywave_his_path: Path,
        detail_model_sw,
    ) -> None:
        """Passing bc_path calls detail.snapwave_boundary_conditions.write()."""
        overall = make_mock_hurrywave(hurrywave_his_path.parent)
        nest2_hurrywave_in_sfincsmodel(
            overall,
            detail_model_sw,
            obs_point_prefix=OBS_PREFIX,
            output_path=str(hurrywave_his_path.parent),
            bc_path=str(tmp_path),
        )
        _assert_sw_timeseries_set(detail_model_sw)


# ---------------------------------------------------------------------------
# nest2_beware_in_sfincsmodel
# ---------------------------------------------------------------------------


class TestNest2SfincsmodelInBeware:
    """Tests for overall=BEWARE, detail=SfincsModel (SfincsModel nested in BEWARE)."""

    def test_basic(
        self,
        beware_his_path: Path,
        detail_model_wl,
    ) -> None:
        """Basic nesting sets water level BCs without error."""
        overall = make_mock_beware(beware_his_path.parent)
        nest2_beware_in_sfincsmodel(
            overall,
            detail_model_wl,
            output_path=str(beware_his_path.parent),
        )
        assert detail_model_wl.water_level.data is not None
        assert "bzs" in detail_model_wl.water_level.data

    def test_correction(
        self,
        beware_his_path: Path,
        detail_model_wl,
    ) -> None:
        """Boundary water level correction is applied uniformly to all BCs.

        The fixture uses WL=0.5 and R2_setup=0.1, so the expected value is 0.6 + correction.
        """
        overall = make_mock_beware(beware_his_path.parent)
        correction = 0.3
        nest2_beware_in_sfincsmodel(
            overall,
            detail_model_wl,
            output_path=str(beware_his_path.parent),
            boundary_water_level_correction=correction,
        )
        bzs = detail_model_wl.water_level.data["bzs"].values  # (time, index)
        np.testing.assert_allclose(bzs[:, 0], 0.6 + correction, rtol=1e-5)

    def test_default_output_path(
        self,
        beware_his_path: Path,
        detail_model_wl,
    ) -> None:
        """output_path defaults to overall.path."""
        overall = make_mock_beware(beware_his_path.parent)
        overall.path = str(beware_his_path.parent)
        nest2_beware_in_sfincsmodel(
            overall,
            detail_model_wl,
        )
        assert "bzs" in detail_model_wl.water_level.data

    def test_option_wave_raises(
        self,
        beware_his_path: Path,
        detail_model_wl,
    ) -> None:
        """NotImplementedError is raised when option != 'flow'."""
        overall = make_mock_beware(beware_his_path.parent)
        with pytest.raises(NotImplementedError):
            nest2_beware_in_sfincsmodel(
                overall,
                detail_model_wl,
                output_path=str(beware_his_path.parent),
                option="wave",
            )

    def test_no_stations_in_bbox_raises(
        self,
        beware_his_path: Path,
        detail_model_wl,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """ValueError is raised when no BEWARE stations fall within the detail bbox."""
        overall = make_mock_beware(beware_his_path.parent)
        # Patch detail.bounds to a region far from the BEWARE fixture points.
        monkeypatch.setattr(
            type(detail_model_wl),
            "bounds",
            property(lambda self: (1_000_000.0, 1_000_000.0, 2_000_000.0, 2_000_000.0)),
        )
        with pytest.raises(ValueError, match="No BEWARE output locations"):
            nest2_beware_in_sfincsmodel(
                overall,
                detail_model_wl,
                output_path=str(beware_his_path.parent),
            )

    def test_bc_path_triggers_write(
        self,
        tmp_path: Path,
        beware_his_path: Path,
        detail_model_wl,
    ) -> None:
        """Passing bc_path calls detail.water_level.write() without error."""
        overall = make_mock_beware(beware_his_path.parent)
        nest2_beware_in_sfincsmodel(
            overall,
            detail_model_wl,
            output_path=str(beware_his_path.parent),
            bc_path=str(tmp_path),
        )
        assert "bzs" in detail_model_wl.water_level.data
