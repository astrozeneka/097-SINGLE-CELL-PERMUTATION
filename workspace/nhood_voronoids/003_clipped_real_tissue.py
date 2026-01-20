import argparse
import json
import os

import numpy as np
import pandas as pd
import freud
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon as ShapelyMultiPolygon
from shapely.ops import unary_union

parser = argparse.ArgumentParser(description="Generate clipped Voronoi tessellation for neighborhood motifs.")
parser.add_argument('--geojson', type=str, help='Input GeoJSON file path.', default="geojson/SP19-4321C1_Cycle1_qptiff_ome_annotations_20260120_104656.geojson")
parser.add_argument('--sample-name', type=str, help='Image ID.', default="XX4321C1")
parser.add_argument("--sample-name-column", type=str, help="Column name for sample name in input CSV.", default="imageid")
parser.add_argument("--input", type=str, help="Input file path.", default="geojson/nhood_input_data_683fa287_clustered.csv")
parser.add_argument('--output', type=str, help='Output file path.', required=False)

parser.add_argument('--region-column', type=str, help="Column name for region name in CSV file", default="Parent Region")
parser.add_argument("--region", type=str, help="Region name to plot", default="Stroma01")

parser.add_argument('--x-coord-column', type=str, default='Centroid X', help="Column name for X coordinate")
parser.add_argument('--y-coord-column', type=str, default='Centroid Y', help="Column name for Y coordinate")
parser.add_argument('--motif-column', type=str, default='motif', help="Column name for motif")
args = parser.parse_args()


def extract_parent_shape(geojson_data, region_name):
    """Extract the parent region polygon from GeoJSON by matching properties.name."""
    parent_region = None
    for feature in geojson_data:
        name = feature['properties'].get('name', None)
        if name == region_name:
            parent_region = feature
            break

    if parent_region is None:
        raise ValueError(f"Region '{region_name}' not found in GeoJSON")

    # Convert to Shapely polygon(s)
    parent_polygons = []
    geometry = parent_region['geometry']
    if geometry['type'] == 'MultiPolygon':
        for polygon in geometry['coordinates']:
            exterior = polygon[0]
            holes = polygon[1:] if len(polygon) > 1 else []
            parent_polygons.append(Polygon(exterior, holes))
    elif geometry['type'] == 'Polygon':
        polygon = geometry['coordinates']
        exterior = polygon[0]
        holes = polygon[1:] if len(polygon) > 1 else []
        parent_polygons.append(Polygon(exterior, holes))

    return unary_union(parent_polygons) if len(parent_polygons) > 1 else parent_polygons[0]


def compute_clipped_voronoi(cells_df, x_col, y_col, motif_col, parent_shape):
    """Compute Voronoi tessellation and clip cells to parent region."""
    # Build points array for freud (3D with z=0)
    points = np.zeros((len(cells_df), 3))
    points[:, 0] = cells_df[x_col].values
    points[:, 1] = cells_df[y_col].values
    points[:, 2] = 0

    motifs = cells_df[motif_col].values

    # Center points for freud (required for periodic box)
    x_center = (points[:, 0].min() + points[:, 0].max()) / 2
    y_center = (points[:, 1].min() + points[:, 1].max()) / 2
    points_centered = points.copy()
    points_centered[:, 0] -= x_center
    points_centered[:, 1] -= y_center

    # Create bounding box with some padding
    x_range = points_centered[:, 0].max() - points_centered[:, 0].min()
    y_range = points_centered[:, 1].max() - points_centered[:, 1].min()
    box = freud.box.Box(Lx=x_range * 1.2, Ly=y_range * 1.2, is2D=True)

    # Compute Voronoi
    voro = freud.locality.Voronoi()
    voro.compute((box, points_centered))

    # Clip each Voronoi cell to the parent region
    clipped_cells = []
    clipped_motifs = []

    for i, cell in enumerate(voro.polytopes):
        # Transform back to original coordinates
        cell_coords = cell[:, :2].copy()
        cell_coords[:, 0] += x_center
        cell_coords[:, 1] += y_center

        # Clip to parent shape
        voronoi_poly = Polygon(cell_coords)
        clipped = voronoi_poly.intersection(parent_shape)

        if not clipped.is_empty:
            clipped_cells.append(clipped)
            clipped_motifs.append(motifs[i])

    return clipped_cells, clipped_motifs


def plot_voronoi(cells, motifs, parent_shape, output_path, title=""):
    """Plot the clipped Voronoi cells colored by motif."""
    fig, ax = plt.subplots(figsize=(12, 12))

    # Build color map using default colormap
    unique_motifs = sorted(set(motifs))
    cmap = plt.cm.get_cmap('tab20', len(unique_motifs))
    motif_colors = {motif: cmap(i) for i, motif in enumerate(unique_motifs)}

    # Draw clipped Voronoi cells
    for cell, motif in zip(cells, motifs):
        color = motif_colors[motif]
        if isinstance(cell, Polygon):
            x, y = cell.exterior.xy
            ax.fill(x, y, edgecolor='black', facecolor=color, alpha=0.7, linewidth=0.1)
        elif isinstance(cell, ShapelyMultiPolygon):
            for poly in cell.geoms:
                x, y = poly.exterior.xy
                ax.fill(x, y, edgecolor='black', facecolor=color, alpha=0.7, linewidth=0.1)

    # Draw parent region boundary
    if isinstance(parent_shape, Polygon):
        x, y = parent_shape.exterior.xy
        ax.plot(x, y, 'gray', linewidth=1.5)
    elif isinstance(parent_shape, ShapelyMultiPolygon):
        for poly in parent_shape.geoms:
            x, y = poly.exterior.xy
            ax.plot(x, y, 'gray', linewidth=1.5)

    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Add scale bar
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    scale_length = 10 ** np.floor(np.log10(x_range / 5))
    if x_range / scale_length > 10:
        scale_length *= 2
    bar_x = xlim[0] + x_range * 0.05
    bar_y = ylim[0] + (ylim[1] - ylim[0]) * 0.05
    ax.plot([bar_x, bar_x + scale_length], [bar_y, bar_y], 'k-', linewidth=3)
    scale_length_um = scale_length * 0.5  # 1 pixel = 0.5 µm
    ax.text(bar_x + scale_length / 2, bar_y + (ylim[1] - ylim[0]) * 0.02,
            f'{int(scale_length_um)} μm', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Add legend
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=motif_colors[m], edgecolor='black', label=m)
                       for m in unique_motifs]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


if __name__ == '__main__':
    # Load GeoJSON
    with open(args.geojson) as f:
        data = json.load(f)

    # List available region names for user reference
    available_regions = [f['properties'].get('name', None) for f in data]
    available_regions = [r for r in available_regions if r is not None]
    print(f"Available regions: {available_regions}")

    # Extract parent shape
    parent_shape = extract_parent_shape(data, args.region)
    print(f"Found region '{args.region}'")

    # Load cell data
    nhood_df = pd.read_csv(args.input)
    nhood_df = nhood_df[nhood_df[args.sample_name_column] == args.sample_name]
    print(f"Loaded {len(nhood_df)} cells for sample '{args.sample_name}'")
    nhood_df = nhood_df[nhood_df[args.region_column] == args.region]
    print(f"{len(nhood_df)} cells remain after filtering for region '{args.region}'")

    if len(nhood_df) == 0:
        print(f"No cells found for sample '{args.sample_name}'. Available samples:")
        all_samples = pd.read_csv(args.input)[args.sample_name_column].unique()
        print(all_samples[:20])
        exit(1)

    # Compute clipped Voronoi
    clipped_cells, clipped_motifs = compute_clipped_voronoi(
        nhood_df,
        args.x_coord_column,
        args.y_coord_column,
        args.motif_column,
        parent_shape
    )
    print(f"Computed {len(clipped_cells)} clipped Voronoi cells")

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        os.makedirs("voronoids", exist_ok=True)
        output_path = f"voronoids/{args.sample_name}_{args.region}.png"

    # Plot and save
    plot_voronoi(
        clipped_cells,
        clipped_motifs,
        parent_shape,
        output_path,
        title=f"{args.sample_name} - {args.region}"
    )
