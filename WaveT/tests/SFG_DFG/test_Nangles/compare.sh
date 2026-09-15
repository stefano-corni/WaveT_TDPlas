#!/bin/bash
# Compare SFG/DFG outputs in out/ against frozen references in ref/.
# (make_sfg_dfg.x writes into out/, so references must not live only there.)
# Prints: >0 passed, 0 could not compare, <0 failed.

threshold=1.0e-8

ref_field="ref/field_mag.dat"
new_field="out/field_mag.dat"
ref_pol="ref/polarization_components.dat"
new_pol="out/polarization_components.dat"
ref_mu="ref/pol_0/run_opt_0_x_0/mu_t_pp.dat"
new_mu="out/pol_0/run_opt_0_x_0/mu_t_pp.dat"

if [ ! -f "$ref_field" ] || [ ! -f "$new_field" ]; then
  echo 0
  exit 0
fi

diff_field=$(awk '
  BEGIN { maxd = 0 }
  FNR == NR {
    if ($1 ~ /^#/ || NF < 3) next
    n++; a[n] = $2; b[n] = $3
    next
  }
  {
    if ($1 ~ /^#/ || NF < 3) next
    i++
    d1 = $2 - a[i]; if (d1 < 0) d1 = -d1
    d2 = $3 - b[i]; if (d2 < 0) d2 = -d2
    if (d1 > maxd) maxd = d1
    if (d2 > maxd) maxd = d2
  }
  END { printf "%.12e", maxd }
' "$ref_field" "$new_field")

diff_pol=0
if [ -f "$ref_pol" ] && [ -f "$new_pol" ]; then
  diff_pol=$(awk '
    BEGIN { maxd = 0 }
    FNR == NR {
      if ($1 ~ /^#/ || NF < 9) next
      n++; for (k = 1; k <= 9; k++) r[n, k] = $k
      next
    }
    {
      if ($1 ~ /^#/ || NF < 9) next
      i++
      for (k = 1; k <= 9; k++) {
        d = $k - r[i, k]; if (d < 0) d = -d
        if (d > maxd) maxd = d
      }
    }
    END { printf "%.12e", maxd }
  ' "$ref_pol" "$new_pol")
fi

diff_mu=0
if [ -f "$ref_mu" ] && [ -f "$new_mu" ]; then
  diff_mu=$(awk '
    BEGIN { maxd = 0 }
    FNR == NR {
      if ($1 ~ /^#/ || NF < 4) next
      n++; for (k = 2; k <= 4; k++) r[n, k] = $k
      next
    }
    {
      if ($1 ~ /^#/ || NF < 4) next
      i++
      for (k = 2; k <= 4; k++) {
        d = $k - r[i, k]; if (d < 0) d = -d
        if (d > maxd) maxd = d
      }
    }
    END { printf "%.12e", maxd }
  ' "$ref_mu" "$new_mu")
fi

awk 'BEGIN {
  thr = '"$threshold"'
  df = '"$diff_field"' + 0
  dp = '"$diff_pol"' + 0
  dm = '"$diff_mu"' + 0
  if (df < thr && dp < thr && dm < thr) print 1
  else print -1
}'
