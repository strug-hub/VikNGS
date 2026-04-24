"""Tolerant text diff for VikNGS CLI output.

Output format (TSV): chr<TAB>pos<TAB>ref<TAB>alt<TAB>pvalue<TAB>test_description
Exact match on chr/pos/ref/alt/description; numeric tolerance on pvalue.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class DiffError:
    line_no: int
    message: str


def diff_pvalues(actual: str, expected: str, *, rtol: float = 1e-9, atol: float = 1e-12) -> list[DiffError]:
    errors: list[DiffError] = []
    a_lines = actual.strip().splitlines()
    e_lines = expected.strip().splitlines()

    if len(a_lines) != len(e_lines):
        errors.append(DiffError(0, f"line count: actual={len(a_lines)} expected={len(e_lines)}"))
        # still compare the overlap so we see what drifted
    for i, (a, e) in enumerate(zip(a_lines, e_lines), start=1):
        errors.extend(_diff_line(i, a, e, rtol=rtol, atol=atol))
    return errors


def _diff_line(line_no: int, actual: str, expected: str, *, rtol: float, atol: float) -> list[DiffError]:
    if actual == expected:
        return []
    af = actual.split("\t")
    ef = expected.split("\t")
    if len(af) != len(ef):
        return [DiffError(line_no, f"column count differs: {len(af)} vs {len(ef)} (a={actual!r}, e={expected!r})")]
    errors: list[DiffError] = []
    for col, (a, e) in enumerate(zip(af, ef)):
        if a == e:
            continue
        # try numeric compare
        try:
            av = float(a)
            ev = float(e)
        except ValueError:
            errors.append(DiffError(line_no, f"col {col} text mismatch: {a!r} vs {e!r}"))
            continue
        if not math.isclose(av, ev, rel_tol=rtol, abs_tol=atol):
            errors.append(DiffError(line_no, f"col {col} numeric mismatch: {av!r} vs {ev!r} (rtol={rtol})"))
    return errors


def format_errors(errors: list[DiffError], max_show: int = 20) -> str:
    if not errors:
        return ""
    head = errors[:max_show]
    msg = [f"{len(errors)} diff(s):"]
    for e in head:
        msg.append(f"  line {e.line_no}: {e.message}")
    if len(errors) > max_show:
        msg.append(f"  ... ({len(errors) - max_show} more)")
    return "\n".join(msg)
