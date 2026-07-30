"""Dependency-light regression test for OpenPlumes mesh extraction.

QGIS is not available in ordinary Python test environments. This test loads
only the pure NumPy constants and functions from the plugin source, preserving
the current single-file plugin architecture while exercising real source code.
"""

import ast
from pathlib import Path
import unittest

import numpy as np


SOURCE_FILE = (
    Path(__file__).resolve().parents[1]
    / "openplumes"
    / "openplumes_algorithm.py"
)
PURE_NAMES = {
    "TETRAHEDRA",
    "_interpolate_edge",
    "_orient_triangle",
    "_triangles_for_tetrahedron",
    "extract_isosurface_marching_tetrahedra",
    "_run_mesh_smoke_test",
}


def load_pure_mesh_namespace():
    tree = ast.parse(SOURCE_FILE.read_text(encoding="utf-8"))
    selected = []
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.FunctionDef)):
            if isinstance(node, ast.Assign):
                names = {
                    target.id for target in node.targets
                    if isinstance(target, ast.Name)
                }
            else:
                names = {node.name}
            if names & PURE_NAMES:
                selected.append(node)

    module = ast.Module(body=selected, type_ignores=[])
    namespace = {"np": np}
    exec(compile(module, str(SOURCE_FILE), "exec"), namespace)
    return namespace


class MarchingTetrahedraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh = load_pure_mesh_namespace()

    def test_smoke_surface_honors_threshold(self):
        stats = self.mesh["_run_mesh_smoke_test"]()
        self.assertGreater(stats["vertices"], 0)
        self.assertGreater(stats["triangles"], 0)

    def test_rejects_threshold_outside_volume_range(self):
        x, y, z = np.mgrid[0:1:2j, 0:1:2j, 0:1:2j]
        with self.assertRaisesRegex(ValueError, "strictly between"):
            self.mesh["extract_isosurface_marching_tetrahedra"](
                x + y + z, x, y, z, iso_value=10.0,
            )

    def test_rejects_mismatched_grid_shape(self):
        x, y, z = np.mgrid[0:1:2j, 0:1:2j, 0:1:2j]
        with self.assertRaisesRegex(ValueError, "equal shapes"):
            self.mesh["extract_isosurface_marching_tetrahedra"](
                x + y + z, x[:, :, :1], y, z, iso_value=1.5,
            )


if __name__ == "__main__":
    unittest.main()
