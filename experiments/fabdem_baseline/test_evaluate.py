import numpy as np
from evaluate import trace_d8
from numpy.testing import assert_array_equal


def test_inclusive_counts_and_connected_cells():
    # Four sources meet the center, which exits east to a terminal cell.
    directions = np.array([[2, 4, 8], [1, 1, 0], [128, 64, 32]], dtype=np.uint8)
    result = trace_d8(directions, 5)
    assert result["cycle_cells"] == 0
    assert result["invalid_cells"] == 0
    assert_array_equal(result["accumulation"], [[1, 1, 1], [1, 8, 9], [1, 1, 1]])
    assert result["connected"].all()


def test_cycle_is_not_a_terminal_and_invalid_is_reported():
    directions = np.array([[1, 16, 3], [64, 64, 0]], dtype=np.uint8)
    result = trace_d8(directions, 5)
    assert result["cycle_cells"] == 2
    assert result["invalid_cells"] == 1
    assert_array_equal(result["connected"], [[False, False, False], [False, False, True]])


def test_nodata_and_off_grid_directions_are_not_cycles():
    directions = np.array([[247, 1], [64, 0]], dtype=np.uint8)
    result = trace_d8(directions, 3)
    assert result["cycle_cells"] == 0
    assert result["invalid_cells"] == 0
    assert_array_equal(result["connected"], [[False, False], [False, True]])
