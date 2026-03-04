import os
from pathlib import Path
import pytest

from mmalignments.models.qc.mosdepth import Mosdepth


def test_mosdepth_builds_command(tmp_path, monkeypatch):
    # Create fake bam
    bam = tmp_path / "sample.bam"
    bam.write_text("fakebam")

    mos = Mosdepth(primary_binary="/usr/bin/mosdepth")

    called = {}

    def fake_run(*args, **kwargs):
        # emulate External.run returning a zero-arg callable
        called['args'] = args
        called['kwargs'] = kwargs

        def _inner():
            # emulate running the command
            called['ran'] = True

        return _inner

    monkeypatch.setattr(mos, "run", fake_run)

    runner = mos.run_mosdepth(prefix=tmp_path / "outprefix", bam_file=bam, threads=2, by=tmp_path / "regions.bed")

    assert callable(runner)
    runner()
    assert called.get('ran', False) is True

    # ensure that run() was called with extra_arguments containing prefix and bam
    extra = called['kwargs'].get('extra_arguments')
    assert extra is not None
    assert str(tmp_path / "outprefix") in extra[0]
    assert str(bam) in extra[-1]
