# -*- coding: utf-8 -*-

"""OpenPlumes QGIS Processing algorithm.

The first working pipeline reads monitoring-well samples, interpolates a 3D
concentration field with SciPy, and extracts a georeferenced triangular
isosurface with a serial, NumPy-only marching-tetrahedra implementation.
"""

import time

import numpy as np
from scipy.interpolate import RBFInterpolator

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPoint,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterEnum,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsLineString,
    QgsPolygon,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsPointXY,
    QgsProcessingException,
    QgsWkbTypes,
)


# Every voxel is divided around the same body diagonal (corner 0 to corner 6).
# Using the same decomposition in every voxel makes the face diagonals agree
# across neighboring voxels, preventing cracks caused by inconsistent splits.
TETRAHEDRA = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)


def _interpolate_edge(point_a, point_b, value_a, value_b, iso_value):
    """Return the XYZ location where an edge crosses ``iso_value``."""

    difference = value_b - value_a
    if np.isclose(difference, 0.0):
        return (point_a + point_b) / 2.0

    fraction = (iso_value - value_a) / difference
    return point_a + fraction * (point_b - point_a)


def _orient_triangle(triangle, inside_points, outside_points):
    """Orient a triangle from higher concentration toward lower concentration."""

    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
    outward = outside_points.mean(axis=0) - inside_points.mean(axis=0)
    if np.dot(normal, outward) < 0:
        triangle[[1, 2]] = triangle[[2, 1]]
    return triangle


def _triangles_for_tetrahedron(points, values, iso_value):
    """Extract zero, one, or two triangles from one tetrahedron.

    Values equal to the threshold are classified as inside. Degenerate
    triangles are filtered by the caller after vertices are deduplicated.
    """

    inside = np.flatnonzero(values >= iso_value)
    outside = np.flatnonzero(values < iso_value)

    if len(inside) == 0 or len(inside) == 4:
        return []

    inside_points = points[inside]
    outside_points = points[outside]

    if len(inside) in (1, 3):
        # One corner on one side and three on the other creates one triangle.
        if len(inside) == 1:
            anchor = inside[0]
            others = outside
        else:
            anchor = outside[0]
            others = inside

        triangle = np.array([
            _interpolate_edge(
                points[anchor], points[other], values[anchor], values[other],
                iso_value,
            )
            for other in others
        ])
        return [_orient_triangle(triangle, inside_points, outside_points)]

    # Two corners on each side create a quadrilateral. Split it into two
    # triangles using a stable edge order around the tetrahedron.
    inside_a, inside_b = inside
    outside_a, outside_b = outside
    p_ac = _interpolate_edge(
        points[inside_a], points[outside_a], values[inside_a],
        values[outside_a], iso_value,
    )
    p_ad = _interpolate_edge(
        points[inside_a], points[outside_b], values[inside_a],
        values[outside_b], iso_value,
    )
    p_bc = _interpolate_edge(
        points[inside_b], points[outside_a], values[inside_b],
        values[outside_a], iso_value,
    )
    p_bd = _interpolate_edge(
        points[inside_b], points[outside_b], values[inside_b],
        values[outside_b], iso_value,
    )

    triangles = [
        np.array([p_ac, p_ad, p_bd]),
        np.array([p_ac, p_bd, p_bc]),
    ]
    return [
        _orient_triangle(triangle, inside_points, outside_points)
        for triangle in triangles
    ]


def extract_isosurface_marching_tetrahedra(
        concentration_volume, grid_x, grid_y, grid_z, iso_value,
        is_canceled=None):
    """Extract a serial triangular isosurface from a vertex-centered grid.

    Parameters are plain NumPy arrays so this function remains independently
    testable. All four arrays must have the same 3D shape. ``grid_x``,
    ``grid_y``, and ``grid_z`` contain real-world coordinates, so returned
    vertices are already georeferenced.

    Returns ``(vertices, faces, stats)``. Faces contain indices into the
    deduplicated vertex array.
    """

    arrays = (concentration_volume, grid_x, grid_y, grid_z)
    if any(array.ndim != 3 for array in arrays):
        raise ValueError("Concentration and coordinate grids must be 3D arrays.")
    if any(array.shape != concentration_volume.shape for array in arrays[1:]):
        raise ValueError("Concentration and coordinate grids must have equal shapes.")
    if any(size < 2 for size in concentration_volume.shape):
        raise ValueError("Each grid dimension must contain at least two points.")
    if not np.isfinite(concentration_volume).all():
        raise ValueError("Concentration volume contains NaN or infinite values.")
    if not np.isfinite(iso_value):
        raise ValueError("Isosurface threshold must be finite.")

    volume_min = float(concentration_volume.min())
    volume_max = float(concentration_volume.max())
    if not volume_min < iso_value < volume_max:
        raise ValueError(
            f"Isosurface threshold {iso_value:g} must be strictly between "
            f"the modeled minimum {volume_min:g} and maximum {volume_max:g}."
        )

    vertices = []
    faces = []
    vertex_lookup = {}
    active_voxels = 0
    tetrahedra_tested = 0

    def vertex_index(point):
        # Adjacent voxels independently calculate the same edge intersection.
        # Rounding provides a stable key so they share one mesh vertex.
        key = tuple(np.round(point, decimals=10))
        index = vertex_lookup.get(key)
        if index is None:
            index = len(vertices)
            vertex_lookup[key] = index
            vertices.append(np.asarray(point, dtype=float))
        return index

    nx, ny, nz = concentration_volume.shape
    for i in range(nx - 1):
        if is_canceled is not None and is_canceled():
            break
        for j in range(ny - 1):
            for k in range(nz - 1):
                corner_indices = (
                    (i, j, k),
                    (i + 1, j, k),
                    (i + 1, j + 1, k),
                    (i, j + 1, k),
                    (i, j, k + 1),
                    (i + 1, j, k + 1),
                    (i + 1, j + 1, k + 1),
                    (i, j + 1, k + 1),
                )
                corner_values = np.array([
                    concentration_volume[index] for index in corner_indices
                ], dtype=float)

                # Fast rejection: a threshold surface exists only if corner
                # values occur on both sides of the requested threshold.
                if corner_values.min() > iso_value:
                    continue
                if corner_values.max() < iso_value:
                    continue

                active_voxels += 1
                corner_points = np.array([
                    [grid_x[index], grid_y[index], grid_z[index]]
                    for index in corner_indices
                ], dtype=float)

                for tetrahedron in TETRAHEDRA:
                    tetrahedra_tested += 1
                    tetrahedron = np.asarray(tetrahedron)
                    triangles = _triangles_for_tetrahedron(
                        corner_points[tetrahedron],
                        corner_values[tetrahedron],
                        iso_value,
                    )
                    for triangle in triangles:
                        face = tuple(vertex_index(point) for point in triangle)
                        if len(set(face)) == 3:
                            faces.append(face)

    vertex_array = np.asarray(vertices, dtype=float).reshape((-1, 3))
    face_array = np.asarray(faces, dtype=np.int64).reshape((-1, 3))
    stats = {
        "active_voxels": active_voxels,
        "tetrahedra_tested": tetrahedra_tested,
        "vertices": len(vertex_array),
        "triangles": len(face_array),
    }
    return vertex_array, face_array, stats


def _run_mesh_smoke_test():
    """Minimal deterministic debug path for the pure geometry extractor."""

    x, y, z = np.mgrid[0:1:2j, 0:1:2j, 0:1:2j]
    volume = x + y + z
    vertices, faces, stats = extract_isosurface_marching_tetrahedra(
        volume, x, y, z, iso_value=1.5,
    )
    if not len(vertices) or not len(faces):
        raise AssertionError("Marching-tetrahedra smoke test produced no surface.")
    if not np.allclose(vertices.sum(axis=1), 1.5, atol=1e-9):
        raise AssertionError("Smoke-test vertices do not honor the threshold.")
    return stats


class OpenPlumesAlgorithm(QgsProcessingAlgorithm):
    """Create a georeferenced 3D contaminant-threshold shell."""

    INPUT = "INPUT"
    CONTAMINANT_ATTRIBUTE = "CONTAMINANT_ATTRIBUTE"
    INTERPOLATOR = "INTERPOLATOR"
    MODEL_BOUNDARY = "MODEL_BOUNDARY"
    ISO_VALUE = "ISO_VALUE"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT,
            self.tr("Monitoring well layer"),
            [QgsProcessing.TypeVectorPoint],
        ))
        self.addParameter(QgsProcessingParameterField(
            self.CONTAMINANT_ATTRIBUTE,
            self.tr("Contaminant"),
            parentLayerParameterName=self.INPUT,
            type=QgsProcessingParameterField.Numeric,
        ))
        self.addParameter(QgsProcessingParameterEnum(
            self.INTERPOLATOR,
            self.tr("Interpolator"),
            options=(
                "RBF (Radial Basis Function)",
                "Kriging (planned)",
                "IDW (planned)",
            ),
            defaultValue=0,
        ))
        self.addParameter(QgsProcessingParameterExtent(
            self.MODEL_BOUNDARY,
            self.tr("Model boundary"),
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.ISO_VALUE,
            self.tr("Concentration threshold"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=100.0,
        ))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT,
            self.tr("Plume shell triangles"),
            type=QgsProcessing.TypeVectorPolygon,
        ))

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException("A monitoring well layer is required.")

        contaminant = self.parameterAsString(
            parameters, self.CONTAMINANT_ATTRIBUTE, context,
        )
        interpolator = self.parameterAsEnum(
            parameters, self.INTERPOLATOR, context,
        )
        model_extent = self.parameterAsExtent(
            parameters, self.MODEL_BOUNDARY, context,
        )
        iso_value = self.parameterAsDouble(parameters, self.ISO_VALUE, context)

        if interpolator != 0:
            raise QgsProcessingException(
                "This prototype currently implements RBF interpolation only."
            )
        if model_extent.isEmpty():
            raise QgsProcessingException("Model boundary must not be empty.")

        field_names = set(source.fields().names())
        required_fields = {"Z", "Depth", contaminant}
        missing_fields = sorted(required_fields - field_names)
        if missing_fields:
            raise QgsProcessingException(
                "Missing required input fields: " + ", ".join(missing_fields)
            )

        # Internal modeling CRS for the prototype.
        # UTM Zone 13N uses meters and is appropriate for the current dummy dataset.
        model_crs = QgsCoordinateReferenceSystem("EPSG:32613")

        source_crs = source.sourceCrs()

        if not source_crs.isValid():
            raise QgsProcessingException("Input layer has an invalid CRS.")

        to_model_crs = QgsCoordinateTransform(
            source_crs,
            model_crs,
            context.transformContext(),
        )

        feedback.pushInfo(f"Input CRS: {source_crs.authid()}")
        feedback.pushInfo(f"Modeling CRS: {model_crs.authid()}")


        interpolation_points = []
        for current, feature in enumerate(source.getFeatures()):
            if feedback.isCanceled():
                return {}

            geometry = feature.geometry()
            if geometry.isNull() or geometry.isEmpty():
                feedback.pushWarning(
                    f"Skipped feature {feature.id()}: missing point geometry."
                )
                continue

            point = geometry.asPoint()

            projected_point = to_model_crs.transform(
                QgsPointXY(point.x(), point.y())
            )
            feedback.pushInfo(
                f"{feature['Well_ID']}: "
                f"original=({point.x():.6f}, {point.y():.6f}) "
                f"projected=({projected_point.x():.2f}, {projected_point.y():.2f})"
            )
            row = np.array([
                projected_point.x(),
                projected_point.y(),
                feature["Z"] - feature["Depth"],
                feature[contaminant],
            ], dtype=float)

            if not np.isfinite(row).all():
                feedback.pushWarning(
                    f"Skipped feature {feature.id()}: non-finite sample data."
                )
                continue
            interpolation_points.append(row)

            if source.featureCount():
                feedback.setProgress(int(20 * (current + 1) / source.featureCount()))

        sample_data = np.asarray(interpolation_points, dtype=float)
        if sample_data.ndim != 2 or sample_data.shape[1:] != (4,):
            raise QgsProcessingException("Sample data must form an N by 4 XYZC array.")
        if len(sample_data) < 4:
            raise QgsProcessingException(
                "RBF interpolation in 3D requires at least four valid samples."
            )

        unique_points = np.unique(sample_data[:, :3], axis=0)
        if len(unique_points) != len(sample_data):
            raise QgsProcessingException(
                "Duplicate XYZ sample locations must be resolved before modeling."
            )

        x_min = model_extent.xMinimum()
        x_max = model_extent.xMaximum()
        y_min = model_extent.yMinimum()
        y_max = model_extent.yMaximum()

        # Prototype grid settings. These remain explicit and reproducible until
        # model-resolution and vertical-boundary controls are added to the UI.
        z_bottom = 0.0
        z_top = 100.0
        grid_x, grid_y, grid_z = np.mgrid[
            x_min:x_max:50j,
            y_min:y_max:50j,
            z_bottom:z_top:50j,
        ]
        prediction_points = np.column_stack((
            grid_x.ravel(), grid_y.ravel(), grid_z.ravel(),
        ))

        feedback.pushInfo(
            f"Sample Z range: "
            f"{sample_data[:, 2].min():.2f} to "
            f"{sample_data[:, 2].max():.2f}"
        )

        feedback.pushInfo(
            f"Grid Z range: "
            f"{grid_z.min():.2f} to "
            f"{grid_z.max():.2f}"
        )

        feedback.pushInfo(f"Valid samples (XYZC): {sample_data.shape}")
        feedback.pushInfo(f"Prediction grid: {grid_x.shape}")
        feedback.pushInfo(f"Prediction locations: {prediction_points.shape}")

        try:
            neighbor_num = 9
            epsilon_num = 1
            rbf = RBFInterpolator(
                sample_data[:, :3], sample_data[:, 3],
                kernel="thin_plate_spline", neighbors=8
                #kernel='linear', neighbors=neighbor_num
                #kernel='cubic', neighbors=neighbor_num
                #kernel='quintic', neighbors=neighbor_num
                #kernel='multiquadric', neighbors=neighbor_num, epsilon=epsilon_num
                #kernel='inverse_multiquadric', neighbors=neighbor_num, epsilon=epsilon_num
                #kernel='inverse_quadratic', neighbors=neighbor_num, epsilon=epsilon_num
                #kernel='gaussian', neighbors=neighbor_num
            )



# ---------------------------------------------------------
# DEBUG: Profile between MW-01 and MW-02
# ---------------------------------------------------------

            t = np.linspace(0.0, 1.0, 25)

            start = sample_data[0, :3]
            end = sample_data[1, :3]

            profile_points = start + t[:, None] * (end - start)

            profile_values = rbf(profile_points)

            feedback.pushInfo("MW-01 → MW-02 profile:")

            for fraction, value in zip(t, profile_values):
                feedback.pushInfo(
                    f"{fraction:5.2f}  {value:8.2f}"
                )

            feedback.pushInfo(
                f"Profile range: "
                f"{profile_values.min():.2f} to "
                f"{profile_values.max():.2f}"
            )

            predicted_values = rbf(prediction_points)
            predicted_at_samples = rbf(sample_data[:, :3])

            feedback.pushInfo(
                f"Max interpolation error: "
                f"{np.max(np.abs(predicted_at_samples - sample_data[:, 3])):.6f}"
            )

        except (ValueError, np.linalg.LinAlgError) as error:
            raise QgsProcessingException(
                f"RBF interpolation failed: {error}"
            ) from error

        concentration_volume = predicted_values.reshape(grid_x.shape)
        volume_min = float(concentration_volume.min())
        volume_max = float(concentration_volume.max())
        feedback.pushInfo(
            f"Modeled concentration range: {volume_min:.3f} to {volume_max:.3f}"
        )

        started = time.perf_counter()
        try:
            vertices, faces, mesh_stats = extract_isosurface_marching_tetrahedra(
                concentration_volume,
                grid_x,
                grid_y,
                grid_z,
                iso_value,
                is_canceled=feedback.isCanceled,
            )
        except ValueError as error:
            raise QgsProcessingException(str(error)) from error

        elapsed = time.perf_counter() - started
        if feedback.isCanceled():
            return {}
        if len(faces) == 0:
            raise QgsProcessingException(
                "No plume shell triangles were produced at this threshold."
            )

        output_fields = QgsFields()
        output_fields.append(QgsField("triangle_id", QVariant.Int))
        output_fields.append(QgsField("threshold", QVariant.Double))
        output_fields.append(QgsField("contaminant", QVariant.String))

        sink, destination_id = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            output_fields,
            QgsWkbTypes.PolygonZ,
            model_crs,
        )
        if sink is None:
            raise QgsProcessingException("Could not create plume shell output.")

        for triangle_id, face in enumerate(faces):
            triangle = vertices[face]

            ring = [QgsPoint(*point) for point in triangle]
            ring.append(QgsPoint(*triangle[0]))

            polygon = QgsPolygon()
            polygon.setExteriorRing(QgsLineString(ring))

            output_feature = QgsFeature(output_fields)
            output_feature.setGeometry(QgsGeometry(polygon))
            output_feature.setAttributes([
                triangle_id,
                iso_value,
                contaminant,
            ])

            sink.addFeature(
                output_feature,
                QgsFeatureSink.FastInsert,
            )

            if triangle_id % 250 == 0:
                feedback.setProgress(
                    80 + int(20 * triangle_id / max(len(faces), 1))
                )      

        feedback.setProgress(100)
        feedback.pushInfo(f"Isosurface threshold: {iso_value:g}")
        feedback.pushInfo(f"Active voxels: {mesh_stats['active_voxels']}")
        feedback.pushInfo(f"Tetrahedra tested: {mesh_stats['tetrahedra_tested']}")
        feedback.pushInfo(f"Unique mesh vertices: {mesh_stats['vertices']}")
        feedback.pushInfo(f"Output triangles: {mesh_stats['triangles']}")
        feedback.pushInfo(f"Serial extraction time: {elapsed:.3f} seconds")

        return {self.OUTPUT: destination_id}

    def name(self):
        return "create_plume"

    def displayName(self):
        return self.tr("Create Plume")

    def group(self):
        return self.tr("Plume Modeling")

    def groupId(self):
        return "plume_modeling"

    def tr(self, string):
        return QCoreApplication.translate("Processing", string)

    def createInstance(self):
        return OpenPlumesAlgorithm()
