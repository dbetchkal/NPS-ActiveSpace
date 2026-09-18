"""Tests for omni source ladder selection."""

from nps_active_space.utils.helpers import get_omni_sources


class TestGetOmniSources:
    def test_default_half_db_step(self) -> None:
        sources = get_omni_sources(0.0, 2.0)
        stems = [s.rsplit("/", 1)[-1].replace(".src", "") for s in sources]
        assert stems == ["O_+000", "O_+005", "O_+010", "O_+015", "O_+020"]

    def test_five_db_step_full_nmsim_span(self) -> None:
        sources = get_omni_sources(-10.0, 40.0, step_db=5.0)
        stems = [s.rsplit("/", 1)[-1].replace(".src", "") for s in sources]
        assert stems == [
            "O_-100",
            "O_-050",
            "O_+000",
            "O_+050",
            "O_+100",
            "O_+150",
            "O_+200",
            "O_+250",
            "O_+300",
            "O_+350",
            "O_+400",
        ]
        assert len(stems) == 11
