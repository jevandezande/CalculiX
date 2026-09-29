#!/usr/bin/env python3.12
"""Timed linear-static block solve for the serial SPOOLES CalculiX build."""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

LENGTH = 10.0
WIDTH = 5.0
HEIGHT = 5.0
YOUNG = 210_000.0
POISSON = 0.3
FORCE = 1_000.0

# Aspect ratio 2:1:1. Baseline on an Apple M4 with 24 GB, serial SPOOLES,
# OpenBLAS pinned to one thread (two runs of `default` agreed within 1%):
#   check    4x2x2       123 equations    0.02 s
#   quick    32x16x16     28k equations    1.2 s    0.3 GiB
#   default  64x32x32    211k equations   44 s      3.8 GiB
# A sample of `default` put about 96% of the samples in spooles_factor and
# about 75% in DVdot33, the dense frontal update inside conda-forge SPOOLES
# (Utilities_DV.o, called from Chv_updateS). Assembly and stress recovery
# were about 2%. Use `default` to judge a change to that kernel.
PRESETS: dict[str, tuple[int, int, int]] = {
    "check": (4, 2, 2),
    "quick": (32, 16, 16),
    "default": (64, 32, 32),
}

# *NODE PRINT writes about seven significant figures, so the check is
# tighter than a wrong boundary condition and looser than the last printed digit.
RELATIVE_TOLERANCE = 1e-6
TIME_RE = re.compile(r"Total CalculiX Time:\s*([0-9.eE+-]+)")
EQUATION_RE = re.compile(
    r"number of equations\s+(\d+)\s+number of nonzero lower triangular matrix elements\s+(\d+)",
    re.MULTILINE,
)
DISPLACEMENT_RE = re.compile(
    r"^\s*(\d+)\s+([+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?)\s+"
    r"([+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?)\s*$",
    re.MULTILINE,
)
RSS_BYTES_RE = re.compile(r"(\d+)\s+maximum resident set size")
RSS_KB_RE = re.compile(r"maximum resident set size:\s*(\d+)")


class BenchmarkError(RuntimeError):
    """CalculiX run failed or the displacement check missed."""


@dataclass(frozen=True, slots=True)
class Mesh:
    """Regular C3D8 block with `nx` by `ny` by `nz` elements."""

    nx: int
    ny: int
    nz: int

    def __post_init__(self) -> None:
        if min(self.nx, self.ny, self.nz) < 1:
            raise ValueError("element counts must be positive")

    @property
    def elements(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def nodes(self) -> int:
        return (self.nx + 1) * (self.ny + 1) * (self.nz + 1)


@dataclass(frozen=True, slots=True)
class RunResult:
    """One finished CalculiX run."""

    mesh: Mesh
    equations: int
    nonzeros: int
    seconds: float
    max_rss_bytes: int
    ux: float
    uy: float
    uz: float


def expected_displacement() -> tuple[float, float, float]:
    """End-corner displacement for uniform traction and Poisson contraction.

    The fixed face restrains `ux` only. One corner pins `uy` and `uz`, and a
    second node on that face restrains `uz`, so the exact linear field is
    admissible. `C3D8` represents it, which makes a wrong answer a failed run
    rather than a looser mesh.

    Returns:
        `(ux, uy, uz)` at the loaded-face corner `(L, W, H)`.
    """
    ux = FORCE * LENGTH / (WIDTH * HEIGHT * YOUNG)
    strain_x = ux / LENGTH
    return ux, -POISSON * strain_x * WIDTH, -POISSON * strain_x * HEIGHT


def node_id(i: int, j: int, k: int, mesh: Mesh) -> int:
    """Return the 1-based id of grid node `(i, j, k)`."""
    return i * (mesh.ny + 1) * (mesh.nz + 1) + j * (mesh.nz + 1) + k + 1


def nodal_force(j: int, k: int, mesh: Mesh) -> float:
    """Return the consistent x-force at end-face node `(j, k)`."""
    share = 1.0
    if j in (0, mesh.ny):
        share *= 0.5
    if k in (0, mesh.nz):
        share *= 0.5
    return FORCE / mesh.elements * mesh.nx * share


def write_input(path: Path, mesh: Mesh) -> None:
    """Write a CalculiX deck whose only result request is one corner node."""
    ux, uy, uz = expected_displacement()
    with path.open("w", encoding="ascii") as deck:
        deck.write(
            "** Serial-SPOOLES speed benchmark.\n"
            "** Uniaxial tension on a regular C3D8 block. The end displacement\n"
            "** is the exact elasticity solution, so this deck is a patch test\n"
            f"** as well as a timing run. Expected tip corner U = {ux:.12e},"
            f" {uy:.12e}, {uz:.12e}.\n"
            f"** Mesh {mesh.nx} x {mesh.ny} x {mesh.nz}"
            f" ({mesh.elements} elements, {mesh.nodes} nodes).\n"
            "*NODE, NSET=NALL\n"
        )
        for i in range(mesh.nx + 1):
            x = LENGTH * i / mesh.nx
            for j in range(mesh.ny + 1):
                y = WIDTH * j / mesh.ny
                for k in range(mesh.nz + 1):
                    z = HEIGHT * k / mesh.nz
                    deck.write(f"{node_id(i, j, k, mesh)}, {x}, {y}, {z}\n")

        deck.write("*ELEMENT, TYPE=C3D8, ELSET=EALL\n")
        element = 1
        for i in range(mesh.nx):
            for j in range(mesh.ny):
                for k in range(mesh.nz):
                    nodes = (
                        node_id(i, j, k, mesh),
                        node_id(i + 1, j, k, mesh),
                        node_id(i + 1, j + 1, k, mesh),
                        node_id(i, j + 1, k, mesh),
                        node_id(i, j, k + 1, mesh),
                        node_id(i + 1, j, k + 1, mesh),
                        node_id(i + 1, j + 1, k + 1, mesh),
                        node_id(i, j + 1, k + 1, mesh),
                    )
                    deck.write(f"{element}, {', '.join(str(n) for n in nodes)}\n")
                    element += 1

        def write_set(name: str, nodes: list[int]) -> None:
            deck.write(f"*NSET, NSET={name}\n")
            for start in range(0, len(nodes), 8):
                deck.write(", ".join(str(n) for n in nodes[start : start + 8]) + "\n")

        fixed = [
            node_id(0, j, k, mesh)
            for j in range(mesh.ny + 1)
            for k in range(mesh.nz + 1)
        ]
        write_set("X0", fixed)
        write_set("PIN", [node_id(0, 0, 0, mesh)])
        write_set("NOTWIST", [node_id(0, mesh.ny, 0, mesh)])
        write_set("TIP", [node_id(mesh.nx, mesh.ny, mesh.nz, mesh)])

        deck.write(
            "*BOUNDARY\n"
            "X0, 1, 1, 0.\n"
            "PIN, 2, 3, 0.\n"
            "NOTWIST, 3, 3, 0.\n"
            "*MATERIAL, NAME=STEEL\n"
            "*ELASTIC\n"
            f"{YOUNG}, {POISSON}\n"
            "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n"
            "*STEP\n"
            "*STATIC\n"
            "*CLOAD\n"
        )
        for j in range(mesh.ny + 1):
            for k in range(mesh.nz + 1):
                deck.write(
                    f"{node_id(mesh.nx, j, k, mesh)}, 1, {nodal_force(j, k, mesh)}\n"
                )
        deck.write("*NODE PRINT, NSET=TIP\nU\n*END STEP\n")


def find_ccx(explicit: Path | None) -> Path:
    """Return the CalculiX executable, preferring an explicit path."""
    if explicit is not None:
        ccx = explicit.expanduser().resolve()
        if not ccx.is_file():
            raise BenchmarkError(f"CalculiX executable not found: {ccx}")
        return ccx
    candidate = Path(__file__).resolve().parents[1] / "build" / "pixi" / "CalculiX"
    if not candidate.is_file():
        raise BenchmarkError(
            f"CalculiX executable not found at {candidate}. Run `pixi run build` first."
        )
    return candidate


def _float_token(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def resident_bytes(time_stderr: str) -> int:
    """Return peak RSS reported by `/usr/bin/time`, or zero when it is absent."""
    mac = RSS_BYTES_RE.search(time_stderr)
    if mac is not None and sys.platform == "darwin":
        return int(mac.group(1))
    gnu = RSS_KB_RE.search(time_stderr)
    if gnu is not None:
        return int(gnu.group(1)) * 1024
    return 0


def launch_command(ccx: Path) -> list[str]:
    """Return the CalculiX command, wrapped in `/usr/bin/time` when that exists."""
    solver = [str(ccx), "-i", "block"]
    timer = Path("/usr/bin/time")
    if not timer.is_file():
        return solver
    if sys.platform == "darwin":
        return [str(timer), "-l", *solver]
    return [str(timer), "-f", "maximum resident set size: %M", *solver]


def parse_run(log_text: str, time_stderr: str, dat_text: str, mesh: Mesh) -> RunResult:
    """Collect timing, matrix size, and the tip displacement from CalculiX output."""
    if "Job finished" not in log_text:
        tail = "\n".join(log_text.splitlines()[-30:])
        raise BenchmarkError(f"CalculiX did not finish.\n{tail}")
    time_match = TIME_RE.search(log_text)
    equation_match = EQUATION_RE.search(log_text)
    displacements = DISPLACEMENT_RE.findall(dat_text)
    if time_match is None or equation_match is None or not displacements:
        raise BenchmarkError("could not parse CalculiX timing or displacement output")
    _node, ux, uy, uz = displacements[-1]
    return RunResult(
        mesh=mesh,
        equations=int(equation_match.group(1)),
        nonzeros=int(equation_match.group(2)),
        seconds=float(time_match.group(1)),
        max_rss_bytes=resident_bytes(time_stderr),
        ux=_float_token(ux),
        uy=_float_token(uy),
        uz=_float_token(uz),
    )


def check_displacement(result: RunResult) -> None:
    """Raise if the tip corner misses the exact elasticity solution."""
    expected = expected_displacement()
    got = (result.ux, result.uy, result.uz)
    scale = max(abs(expected[0]), 1e-30)
    misses = [
        f"{name}: got {value:.8e}, expected {target:.8e}"
        for name, value, target in zip(("ux", "uy", "uz"), got, expected, strict=True)
        if abs(value - target) > RELATIVE_TOLERANCE * scale
    ]
    if misses:
        raise BenchmarkError("displacement check failed\n" + "\n".join(misses))


def run_mesh(
    mesh: Mesh,
    ccx: Path,
    work: Path | None = None,
    keep: bool = False,
) -> RunResult:
    """Generate a deck, run CalculiX, and return the parsed timing result.

    Args:
        mesh: block discretization
        ccx: CalculiX executable
        work: directory for the deck and logs; a temporary directory when omitted
        keep: retain `work` after a successful run

    Returns:
        Parsed equations, wall time, memory, and tip displacement.

    Raises:
        BenchmarkError: solver did not finish or the patch test missed
    """
    owned = work is None
    if work is None:
        work = Path(tempfile.mkdtemp(prefix="ccx-bench-"))
    else:
        work.mkdir(parents=True, exist_ok=True)
    try:
        write_input(work / "block.inp", mesh)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["VECLIB_MAXIMUM_THREADS"] = "1"
        with (work / "block.log").open("w", encoding="utf-8", errors="replace") as log:
            completed = subprocess.run(
                launch_command(ccx),
                cwd=work,
                env=env,
                stdout=log,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        log_text = (work / "block.log").read_text(encoding="utf-8", errors="replace")
        dat_path = work / "block.dat"
        dat_text = (
            dat_path.read_text(encoding="utf-8", errors="replace") if dat_path.is_file() else ""
        )
        if completed.returncode != 0 and "Job finished" not in log_text:
            raise BenchmarkError(
                f"CalculiX exited {completed.returncode}\n{completed.stderr[-2000:]}"
            )
        result = parse_run(log_text, completed.stderr, dat_text, mesh)
        check_displacement(result)
        return result
    finally:
        if owned and not keep:
            shutil.rmtree(work, ignore_errors=True)


def format_result(result: RunResult) -> str:
    """Return one human-readable line for a finished run."""
    rss_gib = result.max_rss_bytes / (1024**3)
    return (
        f"{result.mesh.nx:>4} {result.mesh.ny:>4} {result.mesh.nz:>4}"
        f"  {result.mesh.elements:>8}  {result.equations:>8}  {result.nonzeros:>12}"
        f"  {result.seconds:>8.2f}s  {rss_gib:>5.2f} GiB"
    )


HEADER = (
    f"{'nx':>4} {'ny':>4} {'nz':>4}  {'elements':>8}  {'equations':>8}  "
    f"{'nonzeros':>12}  {'time':>9}  {'rss':>8}"
)


def median(values: list[float]) -> float:
    """Return the median of a non-empty series."""
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def run_repeated(mesh: Mesh, args: argparse.Namespace, ccx: Path) -> list[RunResult]:
    """Run one mesh `args.repeats` times and return every result."""
    keep = args.keep or args.work is not None
    results: list[RunResult] = []
    for _ in range(args.repeats):
        results.append(run_mesh(mesh, ccx, work=args.work, keep=keep))
        print(format_result(results[-1]))
    if len(results) > 1:
        print(f"median {median([item.seconds for item in results]):.2f}s")
    return results


def solve_preset(name: str, args: argparse.Namespace) -> list[RunResult]:
    """Run one named mesh size."""
    if name not in PRESETS:
        raise BenchmarkError(f"unknown preset {name}")
    if args.nx is None:
        mesh = Mesh(*PRESETS[name])
    else:
        mesh = Mesh(args.nx, args.ny, args.nz)
    ccx = find_ccx(args.ccx)
    print(HEADER)
    results = run_repeated(mesh, args, ccx)
    print("patch test passed")
    return results


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "preset",
        nargs="?",
        default="default",
        choices=tuple(PRESETS),
        help="mesh size to run (default: %(default)s)",
    )
    parser.add_argument("--nx", type=int)
    parser.add_argument("--ny", type=int)
    parser.add_argument("--nz", type=int)
    parser.add_argument("--ccx", type=Path)
    parser.add_argument("--work", type=Path, help="keep the deck and logs in this directory")
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument(
        "--sweep",
        action="store_true",
        help="run check, quick, and default in order",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the block benchmark."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.nx is not None and None in (args.ny, args.nz):
        parser.error("--nx, --ny, and --nz must be given together")
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.sweep:
        ccx = find_ccx(args.ccx)
        sizes = [Mesh(*PRESETS["check"]), Mesh(*PRESETS["quick"]), Mesh(*PRESETS["default"])]
        print(HEADER)
        for mesh in sizes:
            run_repeated(mesh, args, ccx)
        print("patch test passed")
        return
    solve_preset(args.preset, args)


if __name__ == "__main__":
    main()
