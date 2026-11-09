import json
from pathlib import Path

import pytest
from analyze import Member, overlap, recover_cycles
from shapely.geometry import LineString, box


def test_recovery_deduplicates_rotations_and_checks_arrows():
    pair: list[Member] = [{"cell": [0, 0], "raw": 8}, {"cell": [0, 1], "raw": 4}]
    assert len(recover_cycles([pair, pair[::-1]], expected=1)) == 1
    with pytest.raises(ValueError, match="closure"):
        recover_cycles([[{"cell": [0, 0], "raw": 8}, {"cell": [0, 1], "raw": 8}]], expected=1)
    with pytest.raises(ValueError, match="population"):
        recover_cycles([pair], expected=294)


def test_footprint_boundary_touch_is_not_every_member():
    cells = [box(0, 0, 30, 30), box(30, 0, 60, 30)]
    river = LineString([(0, -10), (0, 40)])
    result = overlap(cells, [river])
    assert result == [0.0, 30.0]
    assert sum(distance == 0 for distance in result) == 1


def test_historical_inventory_population():
    source = json.loads((Path(__file__).parent / "results/recovered-input.json").read_text())
    assert len(source["initial"]) == 292
    assert len(source["continuations"]) == 5
    initial = recover_cycles(source["initial"], expected=292)
    cycles = recover_cycles(source["initial"] + source["continuations"])
    assert len(cycles) - len(initial) == 2
    assert sum(map(len, cycles)) == 666
    assert any(m["cell"] == [147888, 234986] and m["raw"] == 5 for cycle in cycles for m in cycle)
