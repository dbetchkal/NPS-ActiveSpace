# AAM site diagnostics - run via docs\tmp_aam_site_report.bat
param(
    [Parameter(Mandatory = $true)]
    [string] $Site
)

$ErrorActionPreference = "SilentlyContinue"

$runs = Join-Path $Site "Output_Data\aam\runs"
$log = Join-Path $Site "Output_Data\aam\active_space.log"
$pred = Join-Path $Site "Output_Data\aam\predictions"
$terrain = Join-Path $Site "Input_Data\aam\terrain"
$activ = Join-Path $Site "Output_Data\aam\ACTIVESPACES"

function Write-Section($title) {
    Write-Output ""
    Write-Output "=== $title ==="
}

Write-Output "AAM site report"
Write-Output "Site:     $Site"
Write-Output "Runs:     $runs"
Write-Output "Log:      $log"
Write-Output "Time:     $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"

Write-Section "Paths"
Write-Output "runs exists:        $(Test-Path $runs)"
Write-Output "active_space.log:   $(Test-Path $log)"
Write-Output "predictions dir:    $(Test-Path $pred)"
Write-Output "terrain dir:        $(Test-Path $terrain)"
Write-Output "ACTIVESPACES dir:   $(Test-Path $activ)"

Write-Section "Run directories (Output_Data\aam\runs)"
$runDirs = @(Get-ChildItem $runs -Directory -ErrorAction SilentlyContinue)
Write-Output "run subdir count:   $($runDirs.Count)"
$inpCount = @(Get-ChildItem $runs -Recurse -Filter "scenario.inp" -ErrorAction SilentlyContinue).Count
Write-Output "scenario.inp count: $inpCount"
if ($runDirs.Count -gt 0) {
    Write-Output "sample names (first 10):"
    $runDirs | Select-Object -First 10 -ExpandProperty Name | ForEach-Object { Write-Output "  $_" }
}

Write-Section "ACTIVESPACES layers (geojson count per altitude folder)"
if (Test-Path $activ) {
    Get-ChildItem $activ -Directory -ErrorAction SilentlyContinue | ForEach-Object {
        $n = @(Get-ChildItem $_.FullName -Filter "*.geojson" -ErrorAction SilentlyContinue).Count
        Write-Output "$($_.Name): $n geojson"
    }
} else {
    Write-Output "(ACTIVESPACES not found)"
}

Write-Section "Prediction cache CSV count by altitude prefix"
if (Test-Path $pred) {
    Get-ChildItem $pred -Filter "*.csv" -ErrorAction SilentlyContinue |
        ForEach-Object {
            if ($_.Name -match '^(\d+)m_') { $matches[1] }
            else { "other" }
        } |
        Group-Object | Sort-Object Name |
        ForEach-Object { Write-Output "$($_.Name)m: $($_.Count) csv" }
} else {
    Write-Output "(predictions dir not found)"
}

if (-not (Test-Path $log)) {
    Write-Section "Log"
    Write-Output "(active_space.log not found - skip log sections)"
    exit 0
}

Write-Section 'aam-run batches by altitude (from job name)'
Select-String -Path $log -Pattern '\[aam-run\]' |
    ForEach-Object {
        if ($_.Line -match '(\d+)m_(mesh|contour)') { "$($matches[1])" }
    } |
    Group-Object | Sort-Object { [int]$_.Name } |
    ForEach-Object { Write-Output "$($_.Name)m: $($_.Count) batches" }

Write-Section 'aam-run timing summary by altitude'
$batchRows = Select-String -Path $log -Pattern '\[aam-run\]' |
    ForEach-Object {
        $line = $_.Line
        $alt = $null
        $sec = $null
        if ($line -match '(\d+)m_(mesh|contour)') { $alt = $matches[1] }
        if ($line -match '(\d+\.\d+)s\s+inp=') { $sec = [double]$matches[1] }
        if ($alt -and $sec) {
            [pscustomobject]@{ Alt = $alt; Sec = $sec }
        }
    }
if ($batchRows) {
    $batchRows |
        Group-Object Alt |
        Sort-Object { [int]$_.Name } |
        ForEach-Object {
            $secs = $_.Group | ForEach-Object { $_.Sec }
            $sum = ($secs | Measure-Object -Sum).Sum
            $avg = ($secs | Measure-Object -Average).Average
            $max = ($secs | Measure-Object -Maximum).Maximum
            Write-Output (
                "$($_.Name)m: batches=$($_.Count) totalSec=$([math]::Round($sum, 1)) " +
                "avgSec=$([math]::Round($avg, 1)) maxSec=$([math]::Round($max, 1))"
            )
        }
} else {
    Write-Output '(no parseable aam-run lines)'
}

Write-Section 'ONE TRACK splits by altitude (hop below terrain)'
Select-String -Path $log -Pattern 'split ONE TRACK' |
    ForEach-Object {
        if ($_.Line -match '(\d+)m_') { "$($matches[1])" }
        else { "unknown" }
    } |
    Group-Object | Sort-Object { if ($_.Name -eq 'unknown') { 999999 } else { [int]$_.Name } } |
    ForEach-Object { Write-Output "$($_.Name)m: $($_.Count) split events" }

Write-Section 'last 15 aam-run lines'
Select-String -Path $log -Pattern '\[aam-run\]' |
    Select-Object -Last 15 |
    ForEach-Object { Write-Output $_.Line }

Write-Output ""
Write-Output "=== end ==="
