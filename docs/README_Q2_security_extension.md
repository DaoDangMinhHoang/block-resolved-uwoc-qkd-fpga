# Q2 security extension rev D — DE1 / Cyclone II

> **Integration status.** This layer is merged into the project. The RTL lives in
> `verilog/`, the scripts in `python/`. The integration patch is applied to
> `verilog/top_module.v` and `verilog/uwoc_channel.v` with **four corrections** —
> read `docs/PATCH_top_module.md` § ERRATA before touching either file. The most
> consequential one: the rate presets are at **0x78–0x7F**, not 0x80–0x87, which
> is the qubit-command range. All five testbenches pass under ModelSim ASE,
> including the new `tb_index_align.v`.
>
> **Timing is met.** Quartus II 13.0 on EP2C20F484C7 at `create_clock -period
> 20.000`: worst-case setup slack **+3.689 ns**, Fmax **61.31 MHz**, hold
> +0.445 ns, TNS 0.000 on every corner, 8 248 / 18 752 LE (44 %), 8 / 52
> multipliers. The archive's −9.047 ns / 34.43 MHz is fixed. See
> `docs/PATCH_top_module.md` § ERRATA for how.
>
> <sub>Earlier revisions of this file quoted +3.800 ns / 61.73 MHz / 8 164 LE from
> the 08:57 build. The authority is `verilog/output_files/top_module.sta.summary`
> and `.fit.summary`, rebuilt 2026-09-08 18:00; the numbers above are that build.</sub>
>
> **Verified on hardware 2026-09-08.** The 100-block probe passed all 13 expected
> values, and the full campaign followed: 140 100 blocks × 65 536 = 9.18×10⁹
> attempts in ~15 minutes, `engine_overruns = engine_timeouts = dropped_blocks = 0`.
> Results in `README.md` §6.6 and `BAO_CAO_CAMPAIGN_Q2_2026-09-08.md`. Since then:
> the `--fast off` vs `--fast on` cross-check ran on 2026-09-09 (sifted gains agree
> within ±22 %), L1 was re-measured with 286 000 blocks, and since 2026-09-11 the
> `Y0 = 0` control exists at all five levels.

Additive layer on the RVIF design. `channel_monitor`, `adaptive_controller`,
`per_qubit_reporter` and `uart_reporter` are untouched.

```
rtl/
  security_block_monitor.v    block-resolved affine counts, basis-resolved,
                              window aligned to uwoc_channel's qub_index
  realtime_attempt_engine.v   paced attempt generator: exact f_rep, serial
  block_reporter_bin.v        96-byte packet per block, checksum, watchdog
  throughput_meter.v          on-chip attempts/s -- proves the real-time claim
  tb_security_chain.v         self-checking testbench for monitor + reporter
  PATCH_top_module.md         exact edits to top_module.v and uwoc_channel.v
python/
  send_security_commands.py   configure, rate preset, reset LAST
  capture_security_blocks.py  UART -> CSV with validation gates
  phase2_security_analysis.py debiased defect, BKM, radial-null, GLLP
```

All four modules lint clean under `iverilog -g2005`. The testbench passes. The
analysis was run end-to-end against synthetic captures with known ground truth.

## What changed in rev D

**The request queue is gone from the attempt engine.** The one-bit `pending`
register was worse than useless: in free running mode `request` was asserted
every clock, so every clock spent waiting for the channel incremented
`overrun_count` even though nothing was wrong, and a second paced tick arriving
while `pending` was already set was silently merged rather than queued. A missed
deadline is now counted and forgotten, with no catch-up burst. That is both
simpler and physically honest — if the hardware cannot sustain 10 MHz the right
report is "it ran slower", not a jittered stream that averages out. The block
grid is safe either way, because `qub_index` only advances on `sample_en`.

**Rate presets are decoded before the legacy branch.** `0x80`–`0x87` all have
bit 7 set and the stock decoder sends every such byte to the PC qubit-command
FIFO, so written in the obvious order `--rate` would never have reached the ROM.
`rate_sel` also occupies the three spare bits of `cfg_bus`, so a pacing change
inside a block raises `config_changed`.

**Cumulative diagnostics reset per campaign.** `engine_timeouts`,
`engine_overruns`, `dropped_blocks` and `watchdog_aborts` previously cleared only
on `rst_n`. One early overrun would have made `capture_valid` fail for every
later run until the board was power cycled.

**The throughput meter no longer retains a stale rate.** With presets changing
between campaigns, keeping the old value is wrong metadata: for the first second,
about 153 blocks at 10 Mattempt/s, every packet would have carried the previous
run's rate. It now clears to zero and reports `rate_valid`.

**The monitor re-arms at every boundary.** Staying armed made a configuration
change that landed between two blocks flag the next, perfectly homogeneous,
block as changed.

**GLLP uses the total detection probability, not the sifted gain.** `G_pool` is
`matched clicks / attempts`, smaller than `clicks_all / attempts` by the basis
factor, and the tagged fraction `Delta = p_multi/Q` was therefore overstated by
about 2x — enough to move the intensity optimum. `Q_all_pool` now matches Table V
(1.181e-3 against the published 1.19e-3) while `G_pool` sits at 5.90e-4. The
secret-rate branch also carries the sifting factor q = 1/2, so both returned
quantities are per emitted pulse.

**Confidence intervals carry the null's own uncertainty.** Each bootstrap
replicate is paired with a random null replicate instead of subtracting the null
mean as if it were known exactly. At L1 and L2 the null is the same size as the
signal, so its sampling error is not negligible.

**The per-trial cross-check was removed.** It was being debiased with the
per-attempt null, which is a different quantity. Rather than carry two nulls the
convention is simply fixed at per attempt.

**Basis sanity checks test both the total and each basis.** The old single test
accepted `trials` near N/4, which is exactly what a build with one dead basis
counter produces. Now the total is checked against N/2 and each basis against
N/4.

## What changed in rev C, and why

**`block_synced` is now a pulsed output, not internal state.** The previous
revision cleared the monitor's internal sync flag on the same edge that raised
`block_pulse`; the reporter samples one clock later, so every packet carried
zero and the host gate rejected the entire capture. `tb_security_chain.v` checks
flags bit 4 on every packet specifically to stop that regression returning.

**Mixture weights are attempt counts, not basis-matched trials.**
Proposition 1's `p_j` is the number of physical rounds, so the moment
denominator and the weight must both be `attempts`. Weighting by trials folds
basis-sampling noise into `p_j` for nothing. The per-trial convention is still
reported as `excess_pertrial_over_V` so the choice is visible rather than
buried.

**The GLLP sweep recomputes Q and e at each intensity with Y0 held fixed.** The
previous version scaled Q linearly with µ and held the QBER constant, which
contradicts the mechanism the whole paper rests on: as µ falls the signal falls
but the background does not, so the QBER rises. Now the measured pooled Q and e
are calibrated back to n̄ and e_pol and the sweep is recomputed properly. The
scan also reports `R_GLLP_secret`, the rate after error-correction leakage —
only that may be labelled a secure key rate.

**Unresolved levels are reported as one-sided 95% upper limits.** Proposition 1
gives Δ ≥ 0, so a negative point estimate is estimator fluctuation, not physics.
The figure draws them as limit markers.

**BKM gets the null and the bootstrap too.** It is built from the same noisy bin
means and carries the same bias; debiasing only the exact defect would make the
exact-versus-BKM figure invalid.

**Watchdog covers all three sending states, and UART ownership survives the
packet.** `security_uart_owner = security_fast_mode | sec_reporter_busy` in
top_module, otherwise a `0x70` mid-packet cuts the wire and the watchdog only
rescues the FSM, not the data.

## Real time

Two different claims, and only one of them was true before.

*Bounded control latency* was already true: no processor, no OS, no host in the
decision path.

*Wall-clock fidelity* was not. The published campaign ran at 4.21×10⁸ attempts
in 34 h = 3440 attempts/s against an emulated f_rep of 10 MHz, about **2900×
slower** than the link being modelled. One 5 ms coherence time took 19 s.

`realtime_attempt_engine` drives `uwoc_channel`'s `sample_en` directly and paces
it with a phase accumulator, so `rate_sel = 6` gives exactly one attempt every 5
clocks at 50 MHz. `attempts_per_sec` then reads 10 000 000 and a coherence block
takes 6.55 ms of wall clock. The two time bases in Section III-B collapse into
one, and the adaptive controller closes a 65536-attempt window every 6.55 ms —
a 153 Hz closed loop.

It also removes the throughput problem entirely:

| attempt rate | how | full L1–L5 campaign |
|---|---|---|
| 4 762/s | as shipped, TURBO_SLOT | 1458 h |
| 29 762/s | FAST_SLOT = 80, FSM path | 233 h |
| 10 000 000/s | real time, `--rate 10M` | **0.69 h** |

The engine is **strictly serial** by design: it never pulses `sample_en` for
attempt n+1 before the events of attempt n are out. `uwoc_channel` latches
`qub_index` at sample time (lines 413–417), so an overlapping engine would put
the monitor one index ahead of the `h` that was applied. `index_step_bad` is the
runtime guard. Do not pipeline the engine without exporting a per-attempt index
from the channel.

State plainly in the paper that the accelerated path bypasses the previously
validated OOK serialisation layer. Item 11 of the patch is the cross-check that
defends it.

## Validation

`tb_security_chain.v`, deliberately started mid-block:

```
captured 480 bytes, 5 complete packets
packet 0..4: magic a5 5a  ver 5  flags 75  attempts 64  chk ok
=== ALL CHECKS PASSED ===
```

Flags `0x75` has bit 4 set, i.e. `block_synced = 1` — the regression is fixed and
proven, not argued. The monitor correctly refused to report the partial first
block. `dropped_blocks`, `watchdog_aborts` and `protocol_errors` all zero.

`phase2_security_analysis.py` against synthetic captures generated from the
paper's own channel model:

| level | true Δ/V | recovered | 95% CI | resolved |
|---|---|---|---|---|
| L3 nominal, 20k blocks | 1.93×10⁻⁴ | 2.10×10⁻⁴ | [1.52, 2.85]×10⁻⁴ | yes |
| L4 nominal, 4k blocks | 1.35×10⁻³ | 1.03×10⁻³ | [0.74, 1.51]×10⁻³ | yes |
| L5 nominal, 2k blocks | 3.54×10⁻³ | 4.84×10⁻³ | [3.85, 5.96]×10⁻³ | yes |
| L3, L5 with Y0 = 0 | 0 | ≤4.5×10⁻⁵ | covers zero | **no** |

Every nominal interval covers the truth and all three radial-null controls come
out unresolved, which is exactly what the theorem requires. The BKM quadratic
tracks the exact defect at L1 and L3 and falls short by roughly 2.7× at L5 — the
weak-fluctuation-to-nonlinear transition, produced without tuning.

### Then on the real board, 2026-09-08

Same estimator at 25 m, `phase2_results/security_summary.csv` from the 2026-09-11
analysis (10 groups). L1 is the 286 000-block stage-B capture and its matched
control; every other row is 20 000 blocks. (The first 20 000-block L1 capture gave
8.40×10⁻⁷, below its own null — unresolved, and no longer used.)

| level | Δ/V measured | 95 % CI | resolved | model | BKM/exact |
|---|---|---|---|---|---|
| L1 nominal | 1.531×10⁻⁵ | [1.15, 2.13]×10⁻⁵ | yes | 1.34×10⁻⁵ | 0.99 |
| L2 nominal | 5.59×10⁻⁵ | [3.01, 10.2]×10⁻⁵ | yes | 5.33×10⁻⁵ | 0.97 |
| L3 nominal | 1.441×10⁻⁴ | [1.02, 2.11]×10⁻⁴ | yes | 1.93×10⁻⁴ | 0.83 |
| L4 nominal | 1.127×10⁻³ | [9.68×10⁻⁴, 1.29×10⁻³] | yes | 1.35×10⁻³ | 0.49 |
| L5 nominal | 3.476×10⁻³ | [3.18, 3.81]×10⁻³ | yes | 3.54×10⁻³ | 0.45 |
| L1 `Y0 = 0` | n.r. | [−7.3×10⁻⁷, 1.49×10⁻⁶] | no, UL 1.34×10⁻⁶ | 0 | — |
| L2 `Y0 = 0` | n.r. | [−8.9×10⁻⁶, 2.70×10⁻⁵] | no, UL 2.50×10⁻⁵ | 0 | — |
| L3 `Y0 = 0` | n.r. | [−7.8×10⁻⁶, 2.91×10⁻⁵] | no, UL 2.57×10⁻⁵ | 0 | — |
| L4 `Y0 = 0` | 1.29×10⁻⁵ | [1.9×10⁻⁶, 4.76×10⁻⁵] | marginal, see below | 0 | — |
| L5 `Y0 = 0` | n.r. | [−5.4×10⁻⁶, 2.68×10⁻⁵] | no, UL 2.45×10⁻⁵ | 0 | — |

The model lands inside the interval at four levels; at L4 it sits 3–4.5 % above
the upper bound, depending on the bootstrap run. The `SIGMA2_HO_20M` explanation
offered earlier does not hold: the comment above `SIGMA2_HO_20M` in
`phase2_security_analysis.py` shows that at 25 m the bitstream also draws from the
20 m ROM row (level offset 0), so model and hardware use the same σ². A remaining
candidate is that the model uses the continuous fading law while the ROM holds a
256-entry quantised inverse CDF. Not yet checked; the paper says the offset is
unattributed.

The measured BKM shortfall at L5 is **2.22×**, against 2.7× on synthetic data. Both
show the same transition; quote the measured one.

The L4 control's interval clears zero by a hair (1.29×10⁻⁵, 1.1 % of the nominal
L4 value), which reads as a residual defect with the background switched off.
Across 4/6/8/12/16/24 bins it reads `4.1 → 2.2 → 12.7 → 9.5 → 10.3 → 6.0` (×10⁻⁶,
2026-09-11 run): positive throughout, but without the monotonic climb and
saturation of a real signal. The L1, L2, L3 and L5 controls all cover zero, and
among five controls one marginal exclusion at 95 % happens about 23 % of the time.
Report the controls as bounding any non-`Y0` residual (≤ 4.2×10⁻⁵, i.e. ≤ 3.7 % of
the L4 signal) — **never as zero**, and never as "85× / 3 117×", which divided by a
control point estimate sitting on zero.

## Run order

1. ~~`uwoc_channel` timing fix. Positive slack at 50 MHz before anything else.~~ ✅
2. ~~`evt_bmatch` / `evt_bz` on lost attempts.~~ ✅
3. ~~Testbench must print `ALL CHECKS PASSED`.~~ ✅ 6/6
4. ~~100-block probe. See the table in PATCH item 11.~~ ✅ 2026-09-08, 13/13 values
5. ~~Cross-check one point through the stock OOK path against fast mode.~~ ✅
   2026-09-09 — pooled over L1–L5 the sifted gains agree (OOK/fast 0.985, CI
   [0.786, 1.210]) and QBER differs by 0.30 pp (z = 0.59); resolution ±22 %. Use
   `fpga_collect.py --sections B --scale 0.05`, not `bb84_uwoc_measure.py` (it
   `readline()`s per qubit, 86 q/s), and back up `data/` first: `save_clicks()`
   overwrites.
6. ~~Campaign, then `phase2_security_analysis.py --bin-scan`.~~ ✅ 2026-09-08

## Figures produced

- `defect_vs_turbulence_*.png` — Δ_cg/V per level, nominal against the Y0 = 0
  control, bootstrap error bars, upper-limit markers for unresolved levels.
  The headline figure.
- `GT_turb*.png` — the (G,T) trajectory with the pooled ray drawn. A straight
  line through the origin is radial and carries zero defect; the curvature away
  from it is the entire mechanism.
- `QBER_h_turb*.png` — QBER against irradiance, measurement over model.
- `exact_vs_bkm_*.png` — both debiased, both with intervals.
- `CONTROL_GT_turb*.png` — the Y0 = 0 and Y0 nominal trajectories on the SAME
  axes. With the background off the trajectory is a straight ray through the
  origin and the defect is exactly zero; switching it on bends the curve away
  from the ray. One picture, the whole result.
- `basis_symmetry_*.png` — e_Z against e_X. The emulator applies e_pol
  independently of basis, so these must agree; the figure justifies the reduced
  2-D analysis with data instead of assertion, and shows the basis-resolved
  counts are a validation artefact rather than four independent constraints.
- `gllp_intensity_turb*.png` — click-rate optimum against certified GLLP
  optimum, one figure per operating point, never averaged across configurations.
- `bin_scan.png` — convergence in the number of irradiance bins.

## Pick one repetition rate

The manuscript says f_rep = 10 MHz and also that τ_coh = 5 ms spans "roughly
5×10⁴ pulses, which the emulator rounds to 65536". Both cannot hold: 65536
pulses at 10 MHz is 6.55 ms, and 5 ms at 65536 pulses implies 13.1072 MHz. Keep
**10 MHz**, because it already enters eq. (8) and every SKR in Table III, and
state that the coherence block corresponds to 6.55 ms. `--rate 13M1` exists for
the other convention but rescales every rate in the paper by 1.31.

## Two things not to put in the paper as written

`blocks_planning` in `--model-only` uses an empirical factor of 20 calibrated
against Monte Carlo of this pipeline. It is for planning. If a sample-size
statement is needed, quote measured bootstrap power instead.

> Borne out by the campaign: the table asked for 285 623 / 71 975 / 19 851 / 2 856 /
> 1 089 blocks at L1…L5, and 20 000 blocks resolved L3, L4 and L5 as predicted and
> left L1 unresolved as predicted — but **also resolved L2**, where it had asked for
> 3.6× more. Four of five rows right, one conservative. Do not quote it as a power
> calculation.

`SIGMA2_HO_20M` holds the Table I values, which are quoted at the 20 m
reference while the sweep runs at 25 m. Fine for planning; regenerate from the
Nikishov integral at the actual range before a final figure. The model overlay
is likewise the analytical channel model, not a replay of the RTL's Q24 integer
arithmetic — call it an analytical overlay until you replay the ROM contents.

> **This is no longer hypothetical.** At L4 the model gives 1.346×10⁻³ against a
> measured interval of [9.77×10⁻⁴, 1.307×10⁻³] — 3 % outside. Every other level
> lands inside. A single level missing high, in the one place a 20 m-versus-25 m
> mismatch would bite hardest, is the expected signature. Regenerating σ²_ho at
> 25 m is the first thing to try before treating it as a hardware discrepancy.
