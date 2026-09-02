#!/usr/bin/env python3
"""Time a command and sample its process-tree RSS/CPU for platform compares.

Usage:
  python docs/tmp_platform_run_metrics.py --out metrics.json --label aam -- \\
      python -u nps_active_space/scripts/generate_active_space.py ...
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

SCHEMA = "nps-activespace-platform-metrics/v1"
BYTES_PER_MB = 1024 * 1024
NAMED_PROCESS_FRAGMENTS = (
    "AAM_3.0.0",
    "Nord2000batch",
    "wine",
    "wineserver",
    "python",
)


def _utc_now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _git_head(repo: Path) -> str | None:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip() or None


def _host_info(repo: Path) -> dict[str, object]:
    total_ram_mb: float | None = None
    if sys.platform == "win32":
        raw = os.environ.get("PROCESSOR_IDENTIFIER")
        processor = raw.strip() if raw else None
        total_ram_mb = _windows_total_ram_mb()
    else:
        processor = platform.processor() or None
        try:
            pages = os.sysconf("SC_PHYS_PAGES")
            page_size = os.sysconf("SC_PAGE_SIZE")
            total_ram_mb = round(pages * page_size / BYTES_PER_MB, 1)
        except (ValueError, OSError):
            total_ram_mb = None

    return {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "hostname": platform.node(),
        "cpu_count": os.cpu_count(),
        "processor": processor,
        "total_ram_mb": total_ram_mb,
        "git": _git_head(repo),
        "cwd": str(repo),
        "env": {
            key: os.environ.get(key)
            for key in ("ACOUSTIC_MODEL", "AAM_PARALLEL_N", "AAM_CHUNK_SIZE", "AAM_NC")
            if os.environ.get(key)
        },
    }


def _windows_total_ram_mb() -> float | None:
    payload = _powershell_json(
        "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory",
        timeout_s=15,
    )
    if isinstance(payload, (int, float)):
        return round(float(payload) / BYTES_PER_MB, 1)
    return None


def _powershell_json(command: str, timeout_s: int = 20) -> object | None:
    try:
        proc = subprocess.run(
            [
                "powershell",
                "-NoProfile",
                "-NonInteractive",
                "-Command",
                command,
            ],
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    text = proc.stdout.strip()
    if proc.returncode != 0 or not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        try:
            return float(text)
        except ValueError:
            return text


def _sample_windows(root_pid: int) -> dict[str, object]:
    script = f"""
$ErrorActionPreference = 'SilentlyContinue'
function Get-TreeIds([int]$Id) {{
  $seen = New-Object 'System.Collections.Generic.HashSet[int]'
  $queue = New-Object System.Collections.Queue
  [void]$queue.Enqueue($Id)
  while ($queue.Count -gt 0) {{
    $current = [int]$queue.Dequeue()
    if (-not $seen.Add($current)) {{ continue }}
    Get-CimInstance Win32_Process -Filter "ParentProcessId=$current" |
      ForEach-Object {{ [void]$queue.Enqueue([int]$_.ProcessId) }}
  }}
  return $seen
}}
$ids = Get-TreeIds {root_pid}
$byName = @{{}}
$totalRss = [int64]0
$totalCpu = 0.0
$n = 0
foreach ($id in $ids) {{
  $proc = Get-Process -Id $id
  if (-not $proc) {{ continue }}
  $n++
  $totalRss += $proc.WorkingSet64
  if ($null -ne $proc.CPU) {{ $totalCpu += [double]$proc.CPU }}
  $name = $proc.ProcessName
  if (-not $byName.ContainsKey($name)) {{
    $byName[$name] = @{{ n = 0; rss_bytes = [int64]0 }}
  }}
  $byName[$name].n++
  $byName[$name].rss_bytes += $proc.WorkingSet64
}}
@{{
  rss_bytes = $totalRss
  cpu_s = $totalCpu
  n_procs = $n
  by_name = $byName
}} | ConvertTo-Json -Compress -Depth 6
"""
    payload = _powershell_json(script, timeout_s=25)
    if not isinstance(payload, dict):
        return {"rss_bytes": 0, "cpu_s": 0.0, "n_procs": 0, "by_name": {}}
    return payload


def _linux_children_by_ppid() -> dict[int, list[int]]:
    children: dict[int, list[int]] = {}
    proc_root = Path("/proc")
    if not proc_root.is_dir():
        return children
    for entry in proc_root.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            stat = (entry / "stat").read_text()
        except (OSError, UnicodeDecodeError):
            continue
        rparen = stat.rfind(")")
        if rparen < 0:
            continue
        fields = stat[rparen + 2 :].split()
        if len(fields) < 2:
            continue
        ppid = int(fields[1])
        children.setdefault(ppid, []).append(pid)
    return children


def _linux_tree_ids(root_pid: int) -> set[int]:
    children = _linux_children_by_ppid()
    ids = {root_pid}
    stack = [root_pid]
    while stack:
        current = stack.pop()
        for child in children.get(current, []):
            if child not in ids:
                ids.add(child)
                stack.append(child)
    return ids


def _linux_pid_rss_cpu(pid: int) -> tuple[int, float, str] | None:
    proc_dir = Path("/proc") / str(pid)
    try:
        stat = (proc_dir / "stat").read_text()
        statm = (proc_dir / "statm").read_text().split()
        comm = (proc_dir / "comm").read_text().strip()
    except (OSError, UnicodeDecodeError):
        return None
    rparen = stat.rfind(")")
    fields = stat[rparen + 2 :].split()
    if len(fields) < 13 or len(statm) < 2:
        return None
    ticks = os.sysconf("SC_CLK_TCK")
    page_size = os.sysconf("SC_PAGE_SIZE")
    cpu_s = (int(fields[11]) + int(fields[12])) / ticks
    rss_bytes = int(statm[1]) * page_size
    return rss_bytes, cpu_s, comm


def _sample_linux(root_pid: int) -> dict[str, object]:
    by_name: dict[str, dict[str, int | float]] = {}
    total_rss = 0
    total_cpu = 0.0
    n_procs = 0
    for pid in _linux_tree_ids(root_pid):
        info = _linux_pid_rss_cpu(pid)
        if info is None:
            continue
        rss_bytes, cpu_s, comm = info
        n_procs += 1
        total_rss += rss_bytes
        total_cpu += cpu_s
        bucket = by_name.setdefault(comm, {"n": 0, "rss_bytes": 0})
        bucket["n"] = int(bucket["n"]) + 1
        bucket["rss_bytes"] = int(bucket["rss_bytes"]) + rss_bytes
    return {
        "rss_bytes": total_rss,
        "cpu_s": total_cpu,
        "n_procs": n_procs,
        "by_name": by_name,
    }


def _sample_ps(root_pid: int) -> dict[str, object]:
    """Fallback when /proc is missing (macOS host). RSS from ps is in KiB."""
    try:
        proc = subprocess.run(
            ["ps", "-ax", "-o", "pid=,ppid=,rss=,comm="],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return {"rss_bytes": 0, "cpu_s": 0.0, "n_procs": 0, "by_name": {}}
    rows: dict[int, tuple[int, int, str]] = {}
    children: dict[int, list[int]] = {}
    for line in proc.stdout.splitlines():
        parts = line.split(maxsplit=3)
        if len(parts) < 4:
            continue
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
            rss_kib = int(parts[2])
        except ValueError:
            continue
        rows[pid] = (ppid, rss_kib, parts[3].strip())
        children.setdefault(ppid, []).append(pid)

    ids = {root_pid}
    stack = [root_pid]
    while stack:
        current = stack.pop()
        for child in children.get(current, []):
            if child not in ids:
                ids.add(child)
                stack.append(child)

    by_name: dict[str, dict[str, int]] = {}
    total_rss = 0
    n_procs = 0
    for pid in ids:
        if pid not in rows:
            continue
        _, rss_kib, comm = rows[pid]
        n_procs += 1
        rss_bytes = rss_kib * 1024
        total_rss += rss_bytes
        bucket = by_name.setdefault(comm, {"n": 0, "rss_bytes": 0})
        bucket["n"] += 1
        bucket["rss_bytes"] += rss_bytes
    return {
        "rss_bytes": total_rss,
        "cpu_s": 0.0,
        "n_procs": n_procs,
        "by_name": by_name,
    }


def sample_tree(root_pid: int) -> dict[str, object]:
    if sys.platform == "win32":
        return _sample_windows(root_pid)
    if Path("/proc").is_dir():
        return _sample_linux(root_pid)
    return _sample_ps(root_pid)


def _by_name_mb(by_name: object) -> dict[str, dict[str, float | int]]:
    if not isinstance(by_name, dict):
        return {}
    out: dict[str, dict[str, float | int]] = {}
    for name, raw in by_name.items():
        if not isinstance(raw, dict):
            continue
        rss = float(raw.get("rss_bytes") or 0)
        out[str(name)] = {
            "n": int(raw.get("n") or 0),
            "rss_mb": round(rss / BYTES_PER_MB, 1),
        }
    return out


def _named_rss_mb(by_name_mb: dict[str, dict[str, float | int]]) -> dict[str, float]:
    named: dict[str, float] = {}
    for name, stats in by_name_mb.items():
        if any(fragment.lower() in name.lower() for fragment in NAMED_PROCESS_FRAGMENTS):
            named[name] = float(stats["rss_mb"])
    return named


def _load_metrics(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _write_metrics(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _kill_tree(proc: subprocess.Popen[bytes]) -> None:
    if proc.poll() is not None:
        return
    if sys.platform == "win32":
        subprocess.run(
            ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
            capture_output=True,
            check=False,
        )
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except OSError:
        proc.terminate()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except OSError:
            proc.kill()


def run_timed(
    command: list[str],
    *,
    interval_s: float,
) -> dict[str, object]:
    if command and command[0].endswith(".py"):
        command = [sys.executable, "-u", *command]
    elif command and command[0] in {"python", "python3", "py"}:
        command = [sys.executable, *command[1:]]

    popen_kw: dict[str, object] = {}
    if sys.platform == "win32":
        popen_kw["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        popen_kw["start_new_session"] = True

    started = time.perf_counter()
    proc = subprocess.Popen(command, **popen_kw)
    samples: list[dict[str, object]] = []
    peak_rss_bytes = 0
    peak_n_procs = 0
    peak_named: dict[str, float] = {}
    peak_cpu_s = 0.0
    try:
        while True:
            snapshot = sample_tree(proc.pid)
            rss = int(snapshot.get("rss_bytes") or 0)
            n_procs = int(snapshot.get("n_procs") or 0)
            cpu_s = float(snapshot.get("cpu_s") or 0.0)
            by_name = _by_name_mb(snapshot.get("by_name"))
            named = _named_rss_mb(by_name)
            peak_rss_bytes = max(peak_rss_bytes, rss)
            peak_n_procs = max(peak_n_procs, n_procs)
            peak_cpu_s = max(peak_cpu_s, cpu_s)
            for name, rss_mb in named.items():
                peak_named[name] = max(peak_named.get(name, 0.0), rss_mb)
            samples.append({
                "t_s": round(time.perf_counter() - started, 1),
                "rss_mb": round(rss / BYTES_PER_MB, 1),
                "n_procs": n_procs,
                "cpu_s": round(cpu_s, 1),
                "named_rss_mb": named,
            })
            if proc.poll() is not None:
                break
            time.sleep(interval_s)
    except KeyboardInterrupt:
        _kill_tree(proc)
        raise
    wall_s = time.perf_counter() - started
    return {
        "command": command,
        "started_at": None,
        "wall_s": round(wall_s, 1),
        "exit_code": proc.returncode,
        "peak_rss_mb": round(peak_rss_bytes / BYTES_PER_MB, 1),
        "peak_cpu_s": round(peak_cpu_s, 1),
        "peak_n_procs": peak_n_procs,
        "peak_named_rss_mb": peak_named,
        "sample_interval_s": interval_s,
        "n_samples": len(samples),
        "samples": samples,
    }


def parse_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    if "--" in argv:
        sep = argv.index("--")
        wrapper_argv, command = argv[:sep], argv[sep + 1 :]
    else:
        wrapper_argv, command = argv, []

    parser = argparse.ArgumentParser(
        description="Run a command while sampling process-tree RSS and CPU.",
    )
    parser.add_argument("--out", type=Path, required=True, help="Metrics JSON path.")
    parser.add_argument("--label", required=True, help="Run label, e.g. aam or nmsim.")
    parser.add_argument(
        "--interval",
        type=float,
        default=5.0,
        help="Sampling interval in seconds (default 5).",
    )
    parser.add_argument(
        "--job-json",
        default="",
        help="JSON object of constant job fields stored once on the metrics file.",
    )
    parser.add_argument(
        "--job-file",
        type=Path,
        help="JSON file of constant job fields (same as --job-json).",
    )
    parser.add_argument(
        "--launcher",
        default="",
        help="Platform launcher label stored on the job, e.g. windows-native.",
    )
    args = parser.parse_args(wrapper_argv)
    if not command:
        parser.error("missing command after --")
    return args, command


def main(argv: list[str] | None = None) -> int:
    args, command = parse_args(sys.argv[1:] if argv is None else argv)
    repo = Path.cwd()
    out_path = args.out if args.out.is_absolute() else repo / args.out
    payload = _load_metrics(out_path)
    if payload.get("schema") != SCHEMA:
        payload = {
            "schema": SCHEMA,
            "host": _host_info(repo),
            "job": {},
            "runs": [],
        }
    job: dict[str, object] = {}
    sources = []
    if args.job_file:
        sources.append(args.job_file.read_text(encoding="utf-8"))
    if args.job_json:
        sources.append(args.job_json)
    for raw in sources:
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise SystemExit(f"job metadata is not valid JSON: {exc}") from exc
        if isinstance(parsed, dict):
            job.update(parsed)
    if args.launcher:
        job["launcher"] = args.launcher
    if job:
        existing = payload.get("job")
        merged = dict(existing) if isinstance(existing, dict) else {}
        merged.update(job)
        payload["job"] = merged

    print(f"[metrics] {args.label}: {' '.join(command)}", flush=True)
    started_at = _utc_now()
    result = run_timed(command, interval_s=max(1.0, args.interval))
    result["label"] = args.label
    result["started_at"] = started_at
    result["finished_at"] = _utc_now()

    runs = payload.get("runs")
    if not isinstance(runs, list):
        runs = []
    runs.append(result)
    payload["runs"] = runs
    payload["host"] = payload.get("host") or _host_info(repo)
    _write_metrics(out_path, payload)

    print(
        f"[metrics] {args.label}: wall={result['wall_s']}s "
        f"exit={result['exit_code']} peak_rss={result['peak_rss_mb']}MB "
        f"peak_cpu={result['peak_cpu_s']}s n_procs={result['peak_n_procs']} "
        f"-> {out_path}",
        flush=True,
    )
    return int(result["exit_code"] or 0)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
