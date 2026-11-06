import numpy as np
from numpy.testing import assert_array_equal
from process import route


def test_depression_filled_to_unique_spill_and_counts_conserved():
    dem = np.array([[9, 9, 9], [9, 0, 1], [9, 9, 9]], dtype=np.float32)
    filled, directions, counts = route(dem)
    assert filled[1, 1] == 1
    assert directions[1, 1] == 1  # east to the unique lowest edge
    assert counts[1, 2] == 9
    assert_array_equal(dem, [[9, 9, 9], [9, 0, 1], [9, 9, 9]])


def test_aligned_reprojection_overlap_is_identical():
    from affine import Affine
    from process import project_elevation

    source = np.arange(10000, dtype=np.float32).reshape(100, 100) % 17
    native = Affine(0.0002777777778, 0, -89.34, 0, -0.0002777777778, 31.66)
    small_transform = Affine(30, 0, 279000, 0, -30, 3505000)
    large_transform = small_transform @ Affine.translation(-5, -5)
    small = project_elevation(source, native, "EPSG:4326", -9999, small_transform, (20, 20))
    large = project_elevation(source, native, "EPSG:4326", -9999, large_transform, (30, 30))
    assert_array_equal(small, large[5:25, 5:25])
    assert (small != -9999).all()
