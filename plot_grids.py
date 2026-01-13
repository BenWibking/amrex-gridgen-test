#!/usr/bin/env python3
"""
Plot a 2D projection (x-y) of AMReX grid boxes that intersect a z-plane.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

KPC_IN_CM = 3.08567758e21


def parse_domain_line(line: str) -> tuple[list[int], list[int], int]:
    ints = [int(value) for value in re.findall(r"-?\d+", line)]
    if len(ints) < 7:
        raise ValueError(f"Malformed domain line: {line}")
    lo = ints[0:3]
    hi = ints[3:6]
    ngrid = ints[-1]
    return lo, hi, ngrid


def parse_grid_line(line: str) -> tuple[list[int], list[int]]:
    ints = [int(value) for value in re.findall(r"-?\d+", line)]
    if len(ints) < 6:
        raise ValueError(f"Malformed grid line: {line}")
    lo = ints[0:3]
    hi = ints[3:6]
    return lo, hi


def read_grids(path: Path) -> list[dict[str, object]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"No data found in {path}")
    nlevels = int(lines[0].split()[0])
    levels: list[dict[str, object]] = []
    idx = 1
    for _ in range(nlevels):
        domain_line = lines[idx]
        idx += 1
        domain_lo, domain_hi, ngrid = parse_domain_line(domain_line)
        grids: list[tuple[list[int], list[int]]] = []
        for _ in range(ngrid):
            grid_line = lines[idx]
            idx += 1
            grids.append(parse_grid_line(grid_line))
        levels.append({"domain": (domain_lo, domain_hi), "grids": grids})
    return levels


def read_prob_bounds(path: Path) -> tuple[list[float], list[float]]:
    text = path.read_text()
    lo_match = re.search(r"^\s*geometry\.prob_lo\s*=\s*([^\n#]+)", text, re.MULTILINE)
    hi_match = re.search(r"^\s*geometry\.prob_hi\s*=\s*([^\n#]+)", text, re.MULTILINE)
    if lo_match is None or hi_match is None:
        raise ValueError(f"Missing geometry.prob_lo/prob_hi in {path}")
    prob_lo = [float(value) for value in lo_match.group(1).split()]
    prob_hi = [float(value) for value in hi_match.group(1).split()]
    if len(prob_lo) != 3 or len(prob_hi) != 3:
        raise ValueError(f"Expected 3 values for prob_lo/prob_hi in {path}")
    return prob_lo, prob_hi


def add_box(ax, lo, hi, color, linestyle="-", linewidth=1.0) -> None:
    width = hi[0] - lo[0]
    height = hi[1] - lo[1]
    rect = Rectangle((lo[0], lo[1]), width, height, fill=False, edgecolor=color, linewidth=linewidth, linestyle=linestyle)
    ax.add_patch(rect)


def plot_grids(
    levels: list[dict[str, object]],
    prob_lo: list[float],
    prob_hi: list[float],
    z_phys: float,
    length_scale: float,
    length_label: str,
    show_domain: bool,
    output: Path | None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    color_cycle = plt.cm.tab10.colors
    legend_handles = []
    plotted_count = 0
    min_x = float("inf")
    max_x = float("-inf")
    min_y = float("inf")
    max_y = float("-inf")

    for level_index, level in enumerate(levels):
        color = color_cycle[level_index % len(color_cycle)]
        domain_lo, domain_hi = level["domain"]
        nx = domain_hi[0] - domain_lo[0] + 1
        ny = domain_hi[1] - domain_lo[1] + 1
        nz = domain_hi[2] - domain_lo[2] + 1
        dx = (prob_hi[0] - prob_lo[0]) / nx
        dy = (prob_hi[1] - prob_lo[1]) / ny
        dz = (prob_hi[2] - prob_lo[2]) / nz
        grids = level["grids"]
        plotted = False
        for lo, hi in grids:
            z_lo = (prob_lo[2] + dz * lo[2]) * length_scale
            z_hi = (prob_lo[2] + dz * (hi[2] + 1)) * length_scale
            if z_lo <= z_phys <= z_hi:
                x_lo = (prob_lo[0] + dx * lo[0]) * length_scale
                x_hi = (prob_lo[0] + dx * (hi[0] + 1)) * length_scale
                y_lo = (prob_lo[1] + dy * lo[1]) * length_scale
                y_hi = (prob_lo[1] + dy * (hi[1] + 1)) * length_scale
                add_box(ax, (x_lo, y_lo), (x_hi, y_hi), color=color, linewidth=1.0)
                min_x = min(min_x, x_lo)
                max_x = max(max_x, x_hi)
                min_y = min(min_y, y_lo)
                max_y = max(max_y, y_hi)
                plotted_count += 1
                plotted = True
        if show_domain:
            z_lo = (prob_lo[2] + dz * domain_lo[2]) * length_scale
            z_hi = (prob_lo[2] + dz * (domain_hi[2] + 1)) * length_scale
            if z_lo <= z_phys <= z_hi:
                x_lo = (prob_lo[0] + dx * domain_lo[0]) * length_scale
                x_hi = (prob_lo[0] + dx * (domain_hi[0] + 1)) * length_scale
                y_lo = (prob_lo[1] + dy * domain_lo[1]) * length_scale
                y_hi = (prob_lo[1] + dy * (domain_hi[1] + 1)) * length_scale
                add_box(ax, (x_lo, y_lo), (x_hi, y_hi), color=color, linestyle="--", linewidth=1.5)
                min_x = min(min_x, x_lo)
                max_x = max(max_x, x_hi)
                min_y = min(min_y, y_lo)
                max_y = max(max_y, y_hi)
                plotted = True
        if plotted:
            legend_handles.append(Line2D([0], [0], color=color, lw=2, label=f"Level {level_index}"))

    if plotted_count == 0:
        print(f"No grids intersect z = {z_phys}")
    else:
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)

    ax.set_aspect("equal", "box")
    ax.set_xlabel(f"x ({length_label})")
    ax.set_ylabel(f"y ({length_label})")
    ax.set_title(f"Grid projection at z = {z_phys}")
    if legend_handles:
        ax.legend(handles=legend_handles, loc="best", frameon=False)

    if output is None:
        plt.show()
    else:
        fig.savefig(output, dpi=200, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot x-y projection of AMReX grids at a given physical z.")
    parser.add_argument("grids_file", type=Path, nargs="?", default=Path("DiskGalaxyGrid.grids"))
    parser.add_argument("--inputs", type=Path, default=Path("inputs/DiskGalaxy.in"), help="inputs file with geometry.prob_lo/hi")
    parser.add_argument("--z", type=float, default=0.0, help="physical z coordinate to project onto (default: 0)")
    parser.add_argument(
        "--units",
        choices=("cgs", "kpc"),
        default="cgs",
        help="length units for z and axes (default: cgs)",
    )
    parser.add_argument("--show-domain", action="store_true", help="overlay level domain boxes as dashed outlines")
    parser.add_argument("--output", type=Path, help="save the plot to a file instead of showing it")
    args = parser.parse_args()

    levels = read_grids(args.grids_file)
    prob_lo, prob_hi = read_prob_bounds(args.inputs)
    if args.units == "kpc":
        length_scale = 1.0 / KPC_IN_CM
        length_label = "kpc"
    else:
        length_scale = 1.0
        length_label = "cgs"
    plot_grids(
        levels,
        prob_lo=prob_lo,
        prob_hi=prob_hi,
        z_phys=args.z,
        length_scale=length_scale,
        length_label=length_label,
        show_domain=args.show_domain,
        output=args.output,
    )


if __name__ == "__main__":
    main()
