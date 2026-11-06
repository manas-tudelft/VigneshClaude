# Code guide: `ControlInterface4.py`

A reference for every class, function and method in the script. Line numbers refer to the current version of
the file and will drift as it is edited; the names stay valid. To run the program, see
[HOW_TO_RUN.md](HOW_TO_RUN.md). For the equations behind the functions, see
[MECHANICS_AND_CONTROL.md](MECHANICS_AND_CONTROL.md).

## Contents

1. [Big picture](#1-big-picture)
2. [Conventions](#2-conventions)
3. [Configuration constants](#3-configuration-constants)
4. [Plant dynamics and controller design](#4-plant-dynamics-and-controller-design)
5. [Material analysis functions](#5-material-analysis-functions)
6. [Material model (simulator ground truth and overlay)](#6-material-model)
7. [SimulatedArduino](#7-simulatedarduino)
8. [DataStore and small helpers](#8-datastore-and-small-helpers)
9. [StageGUI](#9-stagegui)
10. [DemoScript](#10-demoscript)
11. [Self-test](#11-self-test)
12. [Entry points](#12-entry-points)
13. [Call flows](#13-call-flows)
14. [The `results` dictionary](#14-the-results-dictionary)
15. [Where to change what](#15-where-to-change-what)

---

## 1. Big picture

Everything is in one file, in this order:

```
Configuration            constants: plant, weights, load cell, modes           (line 46)
Plant dynamics           A,B,C,D matrices, reachability, LQR design            (line 86)
Material analysis        E, k, yield, UTS, hardening, relaxation, regime       (line 306)
Material model           elastoplastic + Maxwell sample (simulator + overlay)  (line 650)
SimulatedArduino         stage + sample + firmware in software                 (line 765)
DataStore + helpers      growable column buffer, small utilities               (line 1071)
StageGUI                 the PyQt5 window: controls, tabs, serial, analysis    (line 1121)
DemoScript               scripted test for --demo                              (line 2189)
Self-test                15 unittest checks for --selftest                     (line 2305)
open_link / main         command line, link to the stage, event loop           (line 2489)
```

How the pieces talk to each other:

```
              serial (omega,I,d,F,u per line)                        every 20 ms
 Arduino ───────────────────────────────────►  read_serial ► handle_line ► DataStore
 (or SimulatedArduino)  ◄────────────────────  send(...)                       │
        gains / RAMP / JOG / PRELOAD / ...                                      │ every 500 ms
                                                                                ▼
        apply_design ◄── maybe_reschedule ◄── update_model ◄── run_analysis (material analysis functions)
            │                                                   │
   design_lqr (dynamics)                                        ▼
                                                  results dict ► plots, Results label, Dynamics tab
```

- **Fast loop (20 ms timer, `on_tick`)**: read serial lines, advance the test program, redraw the signal plots.
- **Slow loop (500 ms timer, `run_analysis`)**: analyse the whole data record, update the measured plant model,
  reschedule the gains, redraw the material and dynamics views.
- **The firmware (or the simulator) runs the control loop at 1 kHz.** The GUI only designs gains and sends
  references; it never closes the loop itself.

## 2. Conventions

| Item | Convention |
|---|---|
| State vector | `x = [ω, I, d, F]` = shaft speed (rad/s), motor current (A), elongation (m), force (N). |
| Input | `u = V_dc · D_cycle` (V), saturated at ±24 V. |
| Augmented state | `[ω, I, d, F, z]` where `z` is the integral of the tracking error. |
| Gain vector | `K = [k_omega, k_I, k_d, k_F, k_int]`; the control law is `u = −(K · [ω, I, e_d, e_F, z])`. |
| Modes | `MODE_DISP = 0` tracks `d`; `MODE_FORCE = 1` tracks `F`. `-1` means "off" in logged data. |
| Units | SI inside the code (m, N, Pa, s). The GUI shows mm, MPa and %. |
| Strain origin | `d0` is the displacement at the end of pre-tension; strain is `(d − d0)/L0`. |
| Sample inputs | Only `A_cs` and `L0`. E, k and everything else are calculated from `d` and `F`. |

## 3. Configuration constants

| Name | Value | Meaning |
|---|---|---|
| `SERIAL_PORT`, `BAUD_RATE` | `'COM9'`, 115200 | Default serial link. |
| `b`, `J`, `k_t`, `L`, `R`, `p`, `V_dc` | see file | Damping, inertia, torque constant, inductance, resistance, lead-screw parameter, supply voltage. |
| `G_SCREW` | `p/π` | Sample elongation per radian of the output shaft (m/rad). |
| `Q_OMEGA`, `Q_CURRENT`, `Q_DISP`, `Q_INT`, `R_INPUT` | 1e-3, 1e-3, 5e8, 3e8, 1e-3 | LQR weights (the existing tuning). |
| `SECONDARY_WEIGHT` | 1e-3 | Weight of the non-tracked output relative to the tracked one. |
| `TAU_PLACEHOLDER` | 1.0 s | Force-state time constant until a relaxation hold has been fitted. |
| `FORCE_LIMIT_N` | 5.0 | LSB200 maximum load; default for the force limit. |
| `PRELOAD_N` | 0.01 | Default preload force. |
| `HX711_RATE_HZ` | 80 | HX711 sample rate (one data line per sample). |
| `MODE_DISP`, `MODE_FORCE`, `MODE_NAMES` | 0, 1, dict | Controller modes and their display names. |
| `C_TRACK` | dict of 1×4 arrays | Selects the tracked output (`d` or `F`) per mode. |
| `C_MEASURED` | 2×4 array | The two physically measured outputs, `d` (encoder) and `F` (load cell). |
| `HOLD, LOAD, RELOAD, UNLOAD` | 0, 1, 2, 3 | Branch labels used by `segment_branches`. |
| `QT_AVAILABLE`, `_MainWindowBase` | – | `True` and `QMainWindow` when PyQt5/pyqtgraph import; otherwise `False` and `object`, so `--selftest` runs without Qt. |

## 4. Plant dynamics and controller design

Lines 86–303. Pure functions of numbers; no Qt, no serial.

### `PlantModel` (dataclass, line 90)

The sample part of the plant, i.e. the force row of `A`. Fields: `k_i` (instantaneous stiffness, N/m), `k_r`
(relaxed or tangent stiffness, N/m), `tau` (relaxation time, s), `regime` (text), `measured` (False until the
stiffness has been measured), `softening` (True when the tangent stiffness is ≤ 0) and `k_e` (measured elastic
stiffness, used to scale the displacement-mode F weight). `PlantModel()` is the *stiffness-free* model
(`k_i = k_r = 0`).

### `plant_matrices(k_i, k_r, tau)` (line 102)

Returns `A, B, C, D` of the 4-state model.

- Rows of `A`: shaft dynamics `J·ω̇`, electrical `L·İ`, kinematics `ḋ = g·ω`, force `Ḟ = k_i·g·ω + (k_r/τ)·d − F/τ`.
- `B = [0, 1/L, 0, 0]ᵀ`, `C = I₄`, `D = 0`.
- Used by: `reachability`, `design_lqr`, the self-test.

### `regime_parameters(regime, k_e=None, k_tan=None, k_v=0.0, tau_v=None)` (line 122)

Maps the *measured* material state to a `PlantModel`.

| Regime | `k_i` | `k_r` |
|---|---|---|
| elastic, unloading, reloading | `k_e + k_v` | `k_e` |
| plastic | `k_tan + k_v` | `k_tan` (floored at `0.01·k_e`) |
| softening | as plastic, with `softening=True` | |

Returns the stiffness-free model if `k_e` is missing or not positive. `tau` falls back to `TAU_PLACEHOLDER`.
Used by `StageGUI.update_model`.

### `_state_scales(k_i)` (line 140)

Typical magnitudes `[1, 0.1, 1e-4, max(k_i,1)·1e-4]` of ω, I, d, F. Scaling the state is needed because the raw
reachability matrix has condition number ≈ 1e20. Used by `reachability` and `observability`.

### `reachability(A, B, model, tol=1e-12)` (line 146)

- Computes the analytic `det W_r = g²·k_t³·(k_i − k_r) / (J³·L⁴·τ)`.
- Runs a scaled PBH test at each eigenvalue.
- Returns `dict(det, rank, uncontrollable, stabilizable, modes)`.
  `rank` is 4 minus the number of uncontrollable modes; `stabilizable` is True when every uncontrollable mode
  is stable.
- For a purely elastic sample (`k_i = k_r`) the determinant is zero, the rank is 3, and the uncontrollable mode
  is `F − k·d` at `−1/τ`.

### `observability(A, C, model, tol=1e-12)` (line 161)

Scaled PBH test for the outputs `C`. Returns `dict(rank, unobservable, detectable)`. With `C_MEASURED`
(encoder and load cell) the rank is 4.

### `lqr_weights(mode, k_i, k_e=0.0)` (line 176)

Builds the 5×5 `Q` for `[ω, I, d, F, z]`.

- **Displacement mode:** `diag(Q_OMEGA, Q_CURRENT, Q_DISP, SECONDARY_WEIGHT·Q_DISP/k_ref², Q_INT)` with
  `k_ref = k_e` (or `k_i` if `k_e` is unknown; the F weight is 0 if both are zero).
- **Force mode:** `diag(Q_OMEGA, Q_CURRENT, SECONDARY_WEIGHT·Q_DISP, Q_DISP/k_i², Q_INT/k_i²)`. Scaling with
  `1/k_i²` keeps the force-loop bandwidth at about 12 rad/s in every regime.
- Raises `ValueError` in force mode when `k_i ≤ 0` (stiffness not measured yet).

### `ControllerDesign` (dataclass, line 192)

The result of a design: `mode`, `model`, `K`, `A`, `B`, `C`, `D`, `poles_open`, `poles_closed`, `reach`, `obs`.

### `design_lqr(mode, model)` (line 206)

LQR with integral action on the tracked output. Builds the augmented system
`[[A, 0], [−C_track, 0]]`, calls `control.lqr` with `lqr_weights` and `R_INPUT`, and returns a
`ControllerDesign` including open- and closed-loop poles, reachability and observability.
Used by `apply_design`, `update_model`, `stiffness_free_design`, the simulator's initial gains and the tests.

### `stiffness_free_design()` (line 222)

The displacement-mode design for `PlantModel()` (nothing known about the sample), cached in `_K0_CACHE`. Its
gains are the `K0` compiled into the sketch; they are stable for stiffness from 10 to 1e6 N/m, so the
controller can run before anything has been measured.

### `format_gain_command(r, K, mode)` (line 231)

Formats the protocol-v2 gain message `r,k_omega,k_I,k_d,k_F,k_int,mode\n` (8 significant digits).

### `dominant_poles(poles, limit=1e3)` (line 236)

The poles with `|Re| < limit` rad/s, sorted slowest first. The electrical and mechanical poles (≈ 1e5 rad/s)
are filtered out so that the slow, interesting poles remain.

### `motor_quantities(omega, current, voltage)` (line 243)

Input power `V·I`, mechanical power `V·I − R·I²` and torque `k_t·I`. If `voltage` is NaN it is reconstructed as
`R·I + k_t·ω`. Used for the optional drive-train force estimate.

### `sample_force(omega, tau_m)` (line 255)

`(π/p)·(τ_m − b·ω)`: the force implied by the motor torque. A cross-check only; gearbox and lead-screw
friction are not modelled. Shown in grey on the force plot when "Show drive-train force estimate" is ticked.

### `describe_design(design, A_cs, L0, K_sent=None)` (line 261)

Builds the text of the Dynamics tab: mode, regime, `k_i`, `k_r`, `τ`, the moduli `E_i`, `E_r`, matrices
(helper `mat`), gains (designed and, if given, sent), open- and closed-loop poles (helper `poles`), dominant
poles, `det W_r`, reachable rank, any uncontrollable mode, and the observable rank. The inner helper `wrap`
wraps long pole lists to 74 characters.

## 5. Material analysis functions

Lines 306–647. All of them work on plain NumPy arrays from `d`, `F`, `A_cs`, `L0`. Nothing is assumed about
the material.

### Preparation

| Function | Does |
|---|---|
| `strain_stress(d, F, A_cs, L0, d0=0)` (313) | Returns `eps, sig, eps_true, sig_true`: engineering strain `(d−d0)/L0`, stress `F/A_cs`, true strain `ln(1+ε)`, true stress `σ(1+ε)`. |
| `windowed_rate(t, y, window=0.3)` (321) | `dy/dt` from the difference over a centred time window; robust against encoder steps. |
| `displacement_quantum(d)` (337) | Smallest displacement step in the data (5th percentile of non-zero steps), i.e. about one encoder count. Used to scale the hold-detection tolerance. |
| `segment_branches(t, d, rate_tol, d_tol, window=0.3)` (344) | Labels every sample `LOAD` (straining on the envelope), `RELOAD` (straining below the previous maximum), `UNLOAD` (d decreasing) or `HOLD` (rate within tolerance). Returns `(labels, rate)`. |

### Elastic properties

| Function | Does |
|---|---|
| `_linfit(x, y)` (358) | Straight-line fit; returns `slope, intercept, R², residuals`. |
| `fit_elastic(eps, sig, d, F, load_mask, window_eps=None)` (366) | **Young's modulus E** (slope of σ–ε) and **stiffness k** (slope of F–d) on the virgin loading branch. Default window: σ between 10 % and 40 % of σ_max on the first rise; or a user window in strain. Returns `dict(E, sigma0, k, F0, r2, n, method, eps_lo, eps_hi, resid_std)` or `None`. E and k are linked by `E = k·L0/A_cs`. |
| `_smooth(y, n=5)` (396) | Moving average that keeps the end points. |
| `proportional_limit(eps, sig, load_mask, fit)` (422) | First sustained deviation (5 points) from the elastic line beyond the fit window. |

### Plastic properties

| Function | Does |
|---|---|
| `offset_yield(eps, sig, load_mask, fit, offset=0.002)` (405) | Offset-method yield: first crossing of the smoothed loading curve with the elastic line shifted by `offset`. Returns `dict(eps, sig, offset, index)`. |
| `ultimate(eps, sig)` (441) | UTS and the strain at UTS; detects fracture when stress falls below 10 % of UTS after the maximum. Returns `dict(uts, eps_u, index, broken[, eps_break, sig_break])`. |
| `tangent_modulus(eps, sig, load_mask, n_eval=150, rel_width=0.03)` (458) | Sliding-window regression slope dσ/dε along the loading envelope (uses cumulative sums, so it is fast). Returns `(strain_centres, slopes)`. |
| `plastic_strain(eps, sig, fit)` (487) | `ε_p = ε − (σ − σ0)/E`. |
| `hardening(eps, sig, load_mask, fit, yld, ult)` (492) | Linear hardening modulus `H = dσ/dε_p` and mean tangent modulus `E_t` between yield and UTS. |
| `hollomon_fit(eps, sig, load_mask, yld, ult)` (509) | Hollomon law `σ_true = K·ε_true^n` by a log–log fit between yield and UTS. |
| `energies(eps, sig, load_mask, fit, yld, ult)` (525) | Modulus of resilience `σ_y²/(2E)` and toughness (area under the loading curve up to fracture). |

### Holds: relaxation and creep

| Function | Does |
|---|---|
| `detect_holds(t, rate, rate_tol, min_duration=2, mask=None)` (542) | Index ranges `(start, end)` where `|rate| ≤ rate_tol` for at least `min_duration` seconds. |
| `_exp_model(tt, y_inf, delta, tau)` (559) | `y_inf + delta·exp(−tt/τ)`. |
| `fit_hold(t, y)` (563) | Fits that exponential with `curve_fit`, from initial guesses taken from the data. For force at constant displacement this is **stress relaxation**; for displacement at constant force it is **creep**. Returns `dict(y0, y0_measured, y_inf, delta, tau, r2, noise, significant, t_start, duration)` or `None`. `significant` requires `|delta| > 3·noise` and `τ < 5·duration`. |
| `relaxation_results(fit, d_hold, d0, A_cs, L0, v_load=None, t_load=None, t_ref=None)` (597) | Turns a force-hold fit into `F0`, `F_inf`, the **relaxation force `dF = F0 − F_inf`**, its percentage, `τ`, the secant moduli `E0`, `E_inf`, and, if the loading rate and duration are known, the Maxwell-arm stiffness `k_v = dF / (v·τ·(1 − e^(−t_r/τ)))` and `E_v = k_v·L0/A_cs`. If `t_ref` is given, `F0` is extrapolated back to the true start of the hold. |

### Regime

`classify_regime(t, d, F, results, eps_now, slack_force, rate_tol, d_tol, window=0.5)` (622) returns the live
regime as text: `fractured`, `unmeasured`, `slack` (mean force below `slack_force`), `hold`, `unloading`,
`elastic`, `plastic` (on the envelope beyond the yield strain) or `softening` (live tangent stiffness < 0).

## 6. Material model

Lines 650–762. Used as the simulator's *ground truth* and, built from measured values only, as the "identified
model" overlay.

### `MaterialParams` (dataclass, line 654)

Sample parameters: `E`, `sigma_y0`, `Q_voce`, `beta_voce`, `H`, `E_v`, `tau_v`, `eps_neck`, `eps_fracture`,
`A_cs`, `L0`, `slack`. The defaults describe a soft polymer film and are **only** the simulator's ground truth;
the GUI never uses them for a real sample.

### `MaterialModel` (class, line 671)

A 1-D elasto-visco-plastic sample: a rate-independent elastoplastic backbone (Voce plus linear isotropic
hardening, radial return) in parallel with a Maxwell arm (`E_v`, `tau_v`). Necking is a damage factor that
reaches fracture; the sample carries no compression.

| Method | Does |
|---|---|
| `__init__(P)` (677) | Stores the parameters and resets the state. |
| `reset()` (681) | Clears plastic strain, hardening variable, Maxwell stress, previous strain, `broken` and stress. |
| `yield_stress(a)` (689) | `σ_y0 + Q·(1 − e^(−β·a)) + H·a` for hardening variable `a`. |
| `step(eps, dt)` (693) | Advances the state to strain `eps`: trial stress, Newton iteration for the plastic increment when over yield, exact update of the Maxwell arm, damage after `eps_neck`, fracture when damage reaches 1. Returns the stress. |
| `force(d, dt)` (733) | Force for absolute grip displacement `d`: `A_cs · step((d − slack)/L0, dt)`. |

### `simulate_material_force(t, d, P)` (line 738)

Runs a fresh `MaterialModel` over a measured displacement history and returns the predicted force.

### `identified_material(res, A_cs, L0, d0)` (line 749)

Builds a `MaterialParams` from the *measured* results only (E from the fit, yield and hardening from the
analysis, `E_v` and `τ` from a significant relaxation fit). Returns `None` if there is no elastic fit. Drives
the magenta dashed "identified model" curve.

## 7. SimulatedArduino

Lines 765–1068. The stage, the sample and the firmware in software. It speaks protocol v2 exactly like the
sketch and offers the parts of `serial.Serial` the GUI uses: `in_waiting`, `readline`, `write`, `close`.

**Class constants:** `TS = 0.002` s (simulated control period), `ENCODER_RES = 1e-7` m, `FORCE_NOISE = 1e-4` N
rms, `HX711_OFFSET = 0.015` N (zero offset before TARE), `FRICTION_CURRENT = 0.005` A (Coulomb friction),
`PRELOAD_SPEED = 1e-5` m/s, `D_LIMIT = 0.02` m, `OMEGA_FILTER = 0.005` s.

**State machine** (`self.state`): `idle`, `run`, `jog`, `preload`, `manual`. A running ramp is
`self.ramp = [target, rate, kind]` with kind `ramp`, `jog` or `jog_idle`.

### Construction and serial interface

| Member | Does |
|---|---|
| `__init__(params=None, speed=1.0, realtime=True, seed=1)` (781) | Creates the material, random generator, state, and starts from the stiffness-free gains `K0`. `realtime=False` is used by the tests (time moves only with `advance`). Prints the boot banner. |
| `in_waiting` (property, 812) | In real-time mode first advances the simulation to match the wall clock, then returns the number of queued lines. |
| `readline()` (817) | Pops the next queued line (bytes) and stores its simulated time in `last_line_time`. |
| `write(data)` (823) | Collects text until newlines and passes each complete line to `_command`. |
| `close()` (831) | No-op. |

### Time

| Member | Does |
|---|---|
| `now()` (835) | Simulated time in seconds. |
| `advance(duration)` (838) | Runs `_step` for `duration` seconds. |
| `_advance_realtime(max_chunk=1.0)` (842) | Catches up with `speed × wall time`; if the GUI fell more than `max_chunk` behind, it slows the simulation down instead of jumping. |
| `drain()` (850) | Returns and clears all queued output lines as text (tests). |
| `_print(msg)` (856) | Queues a text line. |

### Measurement and control law

| Member | Does |
|---|---|
| `_d_meas()` (860) | Encoder displacement, quantised to `ENCODER_RES`. |
| `_x_meas()` (864) | Measured state `(ω_filtered, I_est, d, F)` with `I_est = (u − k_t·ω)/R`. |
| `_law(K, mode, r, z, x)` (static, 869) | The shared control law. Returns `(u, y)` where `u = −(K0·ω + K1·I + K2·e_d + K3·e_F + K4·z)`, with `e_d = d − r, e_F = F` in displacement mode and `e_d = d, e_F = F − r` in force mode, and `y` the tracked output. |
| `_bumpless(K, mode, r, u_now)` (877) | Chooses the integrator `z` so the new law outputs the present `u_now`; no voltage jump at start or gain change. |
| `_stop()` (884) | Back to `idle`; clears ramp and manual voltage. |
| `_control()` (889) | One control step: advances the preload or ramp reference, ends a jog once it has settled (0.2 s within 2 counts, or 5 s timeout), evaluates the law, integrates `z` with **conditional integration** (anti-windup) and saturates `u` to ±24 V. Jog and preload always use the displacement gains `K_disp`. |
| `_hx711_sample()` (929) | One 80 Hz load-cell sample: noisy force, TARE/CAL averaging over 16 samples, force limit, preload completion, displacement limit, and queues the telemetry line `omega,I,d,F,u`. |
| `_step()` (958) | Advances 2 ms: HX711 sample if due, control, then the **quasi-static motor** with a Coulomb friction deadband (`ω = 0` when the drive is below friction), the encoder angle, the sample force from `MaterialModel`, the speed filter and the time. |

### Commands (same set as the sketch)

| Member | Does |
|---|---|
| `_command(line)` (975) | Dispatches `START`, `STOP`, `RESET`, `TARE`, `CAL:`, `FLIM:`, `RAMP:`, `R:`, `JOG:`, `PRELOAD:`, `MANUAL:`; anything else is parsed as a gain message. Bad input prints `ERR bad command: …`. |
| `_gains(line)` (1023) | Parses `r,k_omega,k_I,k_d,k_F,k_int,mode`. If the mode is unchanged and a ramp is running, the ramp keeps its own reference (gain scheduling during a ramp); a mode change cancels the ramp. Switches bumplessly when running. Replies `GAINS OK mode=…`. |
| `_jog(direction, distance, speed)` (1041) | Distance-based move. When running (displacement mode) it starts a `jog` ramp; when idle it enters the `jog` state at the present position. Distance 0 means continuous (target ±`D_LIMIT`). Refused in force mode. |
| `_preload(F_target)` (1060) | Enters `preload`: the reference creeps outwards at `PRELOAD_SPEED` until the measured force reaches `F_target`. |

## 8. DataStore and small helpers

### `DataStore` (class, line 1074)

Growable column buffer with fields `t, omega, I, d, F, u, r, mode`.

| Member | Does |
|---|---|
| `__init__(capacity=8192)` / `clear()` | Allocate (or reallocate) a NaN-filled buffer. |
| `__len__()` | Number of stored samples. |
| `append(row)` | Adds a row, doubling the buffer when full. |
| `arrays()` | Dict of NumPy **views** of each column (no copy). |
| `last(field)` | Newest value of a field, or NaN when empty. |

### Helper functions

| Function | Does |
|---|---|
| `_decimate(n, max_points)` (1104) | A slice that thins `n` points to at most `max_points`, for plotting. |
| `_rel_change(new, old)` (1108) | `|new − old| / max(|old|, 1e-12)`; used to decide whether the model changed by more than 25 %. |
| `_spin(lo, hi, value, step, decimals)` (1112) | Creates a configured `QDoubleSpinBox`. |

## 9. StageGUI

Lines 1121–2186. The window: controls in a scrollable panel on the left, three tabs on the right.

**Class constants:** `RATE_TOL = 5e-7` m/s (hold threshold; raised to 4 encoder counts per window if that is
larger), `RATE_WINDOW = 0.5` s, `D_TOL = 1e-6` m (on-envelope tolerance), `HOLD_MIN_S = 2.0` s,
`CREEP_RATE_TOL = 2e-3` N/s, `PLOT_POINTS = 4000`.

**Main attributes:** `link` (serial port or simulator), `sim` (the simulator or `None`), `data`
(`DataStore`), `events` (time-stamped log), `running`, `jogging`, `ctrl_mode`, `protocol` (`None`, `"v2"` or
`"legacy"`), `d0` and `t_origin` (strain origin), `model` (`PlantModel`), `design` (`ControllerDesign`),
`sent_K`, `r_gui` and `ramp_gui` (the GUI's copy of the reference), `results`, `permanent_set`, `test`
(test-program state), `model_hist`, `pole_hist`, `overlay`.

### Construction and layout

| Method | Does |
|---|---|
| `__init__(link, sim=None, link_name="")` (1133) | Initialises state, builds the UI, starts the 20 ms timer (`on_tick`) and the 500 ms timer (`run_analysis`), draws the first status and Dynamics text, and sends the default force limit. |
| `_group(title)` (1176) | A titled group box with a form layout. |
| `_build_ui()` (1180) | The left panel and all button connections. Group boxes: status labels; **Sample** (A_cs, L0); **Control**; **Test program**; **Load cell**; **Manual control (distance-based)**; **Analysis settings**; **Results**; **Data**; **Arduino log**. |
| `_build_signal_tab()` (1355) | *Signals* tab: d and reference, F (load cell, model, drive estimate), ω, u and estimated I. |
| `_build_material_tab()` (1383) | *Material* tab: stress–strain, force–displacement, tangent modulus, last hold. |
| `_build_dynamics_tab()` (1414) | *Dynamics* tab: text view, k_i/k_r history plot, pole map. |

### Helpers

| Method | Does |
|---|---|
| `clock()` (1442) | Seconds: simulator time in `--sim`, wall time since start on the real stage. |
| `log(text)` (1446) | Appends to the Arduino log box and prints to the console. |
| `send(cmd)` (1450) | Writes a command (adds `\n`) to the link and logs it as `>> …`. |
| `add_event(name)` (1458) | Adds `(clock(), name)` to the event list that goes into the CSV. |
| `sample_geometry()` (1461) | `(A_cs [m²], L0 [m])` from the two sample inputs. |
| `last_d()`, `last_F()` (1465, 1469) | Newest displacement and force (0 if none). |
| `check_protocol()` (1473) | True once a v2 data line has been seen; otherwise shows a warning dialog. |
| `_update_hint_labels(*_)` (1482) | Shows the strain, force and strain-rate equivalents next to the inputs. |
| `_mode_changed(idx)` (1489) | Enables the displacement or stress input according to the mode. |
| `reference_from_inputs(mode)` (1493) | Reference to send: `d0 + d_input` in displacement mode, `σ·A_cs` in force mode. |
| `current_reference(t)` (1499) | The reference the controller should have at time `t`, following the GUI's copy of a running ramp. Used to log `r` with every sample. |
| `update_status_display()` (1507) | Status label (ACTIVE / JOG / OFF), start-button text, and the regime/model line. |

### Controller

| Method | Does |
|---|---|
| `apply_design(mode, r, reason, keep_ramp=False)` (1522) | Designs the LQR for the current measured model with `design_lqr`, sends the gain message, records the poles and an event, refreshes the Dynamics tab. Refuses force mode when the sample is softening, and shows a message if force mode is not yet possible. Returns True on success. |
| `send_parameters()` (1547) | Button *Send Parameters*: `apply_design` for the selected mode and reference. |
| `toggle_controller()` (1553) | Start/Stop button. |
| `start_controller()` (1559) | If no gains were sent yet, first sends the stiffness-free displacement hold at the present position; then `START`. |
| `stop_controller()` (1569) | `STOP`; clears running, jogging, test and ramp state. |
| `update_reference()` (1577) | Sends `R:<r>` for the current inputs (asks first if the mode was changed). |
| `maybe_reschedule()` (1587) | Gain scheduling: if the model changed (`schedule_pending`), auto scheduling is on, the controller is running and at least 1 s has passed, resends the gains designed on the latest model (bumpless in the firmware). If the sample started softening while in force mode, it switches to displacement hold. |

### Test program

| Method | Does |
|---|---|
| `_ensure_displacement_control(reason)` (1605) | Makes sure displacement control is running, starting bumplessly at the present position. |
| `start_ramp_hold()` (1618) | Button *Ramp → Hold*: sends `RAMP:<target>,<rate>` and starts the test state machine. |
| `start_unload()` (1632) | Button *Unload*: ramps the reference back to `d0`. |
| `on_ramp_done()` (1645) | Called on the firmware's `RAMP done`; starts the hold timer, or finishes an unload. |
| `finish_unload(why)` (1657) | Holds the present position and records the **permanent set** `(d − d0)/L0`. |
| `update_test_program()` (1667) | Every tick: ends the hold when its time is over, and ends an unload when the force falls to the preload. |

### Load cell and manual control

| Method | Does |
|---|---|
| `tare()` (1679) | Sends `TARE`. |
| `calibrate()` (1683) | Converts the entered mass to newtons (`m·g`), asks for confirmation, sends `CAL:<N>`. |
| `pre_tension()` (1694) | Sends `PRELOAD:<N>`; the strain origin is set when the firmware reports completion. |
| `set_force_limit()` (1703) | Sends `FLIM:<N>`. |
| `set_strain_origin()` (1707) | Stores the present `d` and time as the strain origin `d0`, `t_origin`; sets the displacement input to 0. |
| `jog(direction)` (1716) | Sends `JOG:<dir>,<distance m>,<speed m/s>`; if the controller is running, mirrors the move in the GUI's reference. Refused in force mode. |
| `stop_manual_control()` (1734) | Same as `stop_controller()`. |

### Data

| Method | Does |
|---|---|
| `reset_data()` (1738) | Clears data, events, results and the model, locks force mode again, sends `RESET`, empties all curves. |
| `save_data_to_file(path=None)` (1765) | Writes the CSV (header with A_cs, L0, `d0`, calculated properties, events; then the 13 data columns). Without a path it asks with a file dialog. |

### Serial and message handling

| Method | Does |
|---|---|
| `on_tick()` (1803) | 20 ms: `read_serial`, `update_test_program`, `update_signal_plots`. |
| `read_serial()` (1808) | **Simulator:** reads complete lines with their own simulated time. **Real port:** non-blocking read into a buffer, splits at newlines, keeps a partial last line, and spreads the time stamps of the data lines at 1/80 s so they stay monotonic. |
| `handle_line(line, t)` (1843) | Parses a line with 5 numbers (v2) or 2 numbers (legacy `omega,d`) into the `DataStore`, tagged with the current reference and mode; anything else goes to `handle_message`. |
| `handle_message(line)` (1868) | Logs the text and reacts: *Pre-tensioning complete* → sets the strain origin; `RAMP done` → `on_ramp_done`; `JOG done`; `FORCE LIMIT` / `DISPLACEMENT LIMIT` → aborts the test; banner with "protocol v2" → marks v2. |

### Analysis (every 500 ms)

| Method | Does |
|---|---|
| `fit_window()` (1894) | The user's E-fit window as strain limits, or `None` for automatic. |
| `live_tangent_stiffness(t, d, F, load)` (1901) | Slope of F against d over the last 3 s of loading; used to detect plastic or softening behaviour. |
| `_loaded_before(t, d, i0, quantum)` (1907) | True if the displacement rose over the 2 s before index `i0` (a hold after loading, not after unloading). |
| `analyse_holds(t, d, F, mode, labels, rate, rate_tol, quantum, A_cs, L0)` (1912) | Finds holds: a displacement hold that follows loading is a **stress relaxation** hold; a force hold in force mode is a **creep** hold. Fits the latest one with `fit_hold`, finds the loading ramp before it for `v_load` and `t_load`, and returns the result dictionary (see section 14). |
| `run_analysis()` (1949) | The main analysis: selects data since the strain origin, computes strain/stress, labels branches, runs every analysis function of section 5, classifies the regime, stores the `results` dictionary, then calls `update_model`, `maybe_reschedule`, the results text, the plot updates, and handles fracture (stops the motor once). Skipped until at least 20 samples and a v2 stream exist. |
| `update_model(res)` (1999) | Builds the measured plant model with `regime_parameters` from the fitted k, the live tangent stiffness and the relaxation fit. If it changed by more than 25 % (or measured/softening flags changed) it stores it, appends to `model_hist`, unlocks force mode, redesigns the controller for display and sets `schedule_pending`. |
| `results_text_lines(res=None)` (2029) | The result lines shown in the Results box and written into the CSV header. Quantities not yet available say "not measured yet". |

### Plots

| Method | Does |
|---|---|
| `update_signal_plots()` (2085) | Redraws the Signals tab from decimated data; updates the live force label. |
| `update_material_plots(d, F, eps, sig, eps_t, sig_t, res)` (2110) | Redraws σ–ε (with E fit, offset line, yield/UTS/fracture markers, optional true curve), F–d (with the stiffness line), the tangent-modulus curve and the last hold with its exponential fit. |
| `update_overlay(a, res, A_cs, L0)` (2156) | Every 2 s of wall time, simulates the identified material over the measured displacement to obtain the "identified model" force curve. |
| `update_dynamics_view()` (2173) | Refreshes the Dynamics text (only when changed), the k_i/k_r history and the pole map. |

## 10. DemoScript

Lines 2192–2302. Started by `--demo`. A list of step methods is polled every 100 ms; each returns True when
done. It acts only through the real button handlers.

| Step | Does |
|---|---|
| `__init__`, `start`, `_elapsed`, `_tick` (2196–2215) | Set up the step list and a timer; run the current step and advance when it returns True. |
| `_jog`, `_wait_jog` | Jog 0.005 mm outwards at 0.05 mm/s and wait for it to finish. |
| `_tare` | Tare. |
| `_preload`, `_wait_preload` | Pre-tension and wait for the strain origin. |
| `_ramp_hold`, `_wait_test` | Ramp to 0.5 mm at 0.01 mm/s, hold 30 s. |
| `_unload` | Unload at 0.02 mm/s. |
| `_fracture`, `_wait_fracture` | Ramp to 2.0 mm at 0.02 mm/s until the sample breaks. |
| `_report` (2276) | Runs a final analysis; prints the calculated properties and the simulator's ground truth with the expected relaxation force; if a folder was given, saves three screenshots and the CSV; then quits the application. |

## 11. Self-test

Lines 2305–2486. Run with `--selftest`. Uses `unittest`; no window, no serial port.

**Helpers:** `_elastic_params(**kw)` (a purely elastic ground truth for exact checks), `_data_lines(lines)`
(parses simulator output into an array).

**`TestDynamics`**

| Test | Checks |
|---|---|
| `test_reachability_determinant` | Analytic `det W_r` equals the numeric determinant. |
| `test_elastic_mode_uncontrollable_but_stable` | Elastic case has rank 3 with the uncontrollable mode at `−1/τ`; viscoelastic case has rank 4. |
| `test_lqr_stable_in_all_regimes` | Closed loop is stable in both modes for elastic, viscoelastic and plastic models; force-mode dominant pole near −12 rad/s. |
| `test_stiffness_free_design_is_robust` | `K0` stabilises plants from 10 to 1e6 N/m. |
| `test_force_mode_needs_measured_stiffness` | `design_lqr(MODE_FORCE, PlantModel())` raises `ValueError`. |
| `test_observable_from_encoder_and_load_cell` | Observable rank is 4 from `[d, F]`. |

**`TestMaterialAnalysis`**

| Test | Checks |
|---|---|
| `test_elastic_yield_hardening_stiffness` | From synthetic elastoplastic data: E, k, yield stress, H and E_t within 2 %; `E = k·L0/A_cs`. |
| `test_relaxation_fit` | A noisy synthetic hold returns τ and ΔF within 2 %. |
| `test_regime_parameters` | No stiffness → not measured; plastic regime gives `k_i = k_tan + k_v`, `k_r = k_tan`. |

**`TestClosedLoop`** (simulated stage in virtual time). Helpers: `_sim(params)` creates a non-real-time
`SimulatedArduino` and discards its banner; `_cmd(sim, text, duration)` sends one command, advances the
simulation by `duration` seconds and returns the output lines.

| Test | Checks |
|---|---|
| `test_displacement_step` | A 0.1 mm step is reached within 0.3 µm. |
| `test_force_step_with_measured_stiffness` | A 0.2 N force step is reached within 2 mN. |
| `test_preload_stops_at_force` | Preload stops near 0.01 N, controller idle. |
| `test_force_limit_cuts_motor` | `FLIM` stops the motor and reports `FORCE LIMIT`. |
| `test_jog_distance` | A 0.1 mm jog moves 0.1 mm within 0.3 µm and ends idle. |
| `test_ramp_rate` | A ramp reaches its target at the commanded rate within 2 %. |

`run_selftest()` (2480) loads the three test cases, runs them verbosely and returns exit status 0 or 1.

## 12. Entry points

### `open_link(args)` (line 2492)

Returns `(link, sim, name)`. For `--sim` or `--demo` it creates a `SimulatedArduino`. Otherwise it imports
`pyserial`, opens the port, waits 2 s for the Uno reset and returns it. If that fails, it asks whether to
start the simulator instead; answering No returns `(None, None, "")` and the program exits.

### `main(argv=None)` (line 2511)

Parses the command line (`--port`, `--baud`, `--sim`, `--sim-speed`, `--demo`, `--screenshot`, `--selftest`),
runs the self-test and returns if asked, checks that Qt is installed, sets the demo speed to 20× if left at 1,
creates the `QApplication`, opens the link, shows `StageGUI`, starts the `DemoScript` for `--demo`, runs the
event loop, and closes the link on exit.

## 13. Call flows

### One data sample

```
HX711 sample ready (80 Hz) ► firmware prints "omega,I,d,F,u"
► on_tick ► read_serial ► handle_line ► DataStore.append((t, omega, I, d, F, u, r, mode))
```

### Analysis and gain scheduling (every 500 ms)

```
run_analysis
 ├─ strain_stress, displacement_quantum, segment_branches
 ├─ fit_elastic ► offset_yield, proportional_limit, ultimate, tangent_modulus,
 │               hardening, hollomon_fit, energies
 ├─ analyse_holds ► detect_holds, fit_hold, relaxation_results
 ├─ classify_regime
 ├─ update_model ► regime_parameters ► (model changed > 25 %) ► design_lqr, schedule_pending = True
 ├─ maybe_reschedule ► apply_design ► design_lqr ► send(format_gain_command(...))
 └─ results_text_lines, update_material_plots, update_overlay, update_dynamics_view
```

### Ramp → Hold test

```
Ramp -> Hold button ► start_ramp_hold ► _ensure_displacement_control ► send("RAMP:target,rate")
firmware ramps the reference at 1 kHz ► prints "RAMP done" ► handle_message ► on_ramp_done
► test stage "hold" ► update_test_program waits hold seconds ► "HOLD end"
(meanwhile run_analysis detects the hold and fits the relaxation)
```

### Pre-tension

```
pre_tension ► send("PRELOAD:0.01") ► firmware creeps outwards until F >= 0.01 N
► prints "Pre-tensioning complete, controller stopped" ► handle_message ► set_strain_origin
```

### Force mode

```
fit_elastic gives k ► update_model ► model.measured = True ► force mode unlocked in the combo box
Send Parameters ► apply_design(MODE_FORCE, sigma*A_cs) ► design_lqr (Q scaled with 1/k_i^2)
► gains sent ► Start ► firmware tracks F
(plastic: the tangent stiffness falls ► rescheduling redesigns K ► same ~12 rad/s force bandwidth)
```

## 14. The `results` dictionary

Stored in `StageGUI.results` by `run_analysis`; every entry may be `None` until enough data exists.

| Key | Content |
|---|---|
| `fit` | `E` (Pa), `sigma0`, `k` (N/m), `F0`, `r2`, `n`, `method`, `eps_lo`, `eps_hi`, `resid_std` |
| `yield` | `eps`, `sig`, `offset`, `index` |
| `prop` | `eps`, `sig` (proportional limit) |
| `ult` | `uts`, `eps_u`, `index`, `broken`, and if broken `eps_break`, `sig_break` |
| `tangent` | `(strain_centres, E_t_values)` |
| `hard` | `H`, `E_t`, `r2`, `n` |
| `hollomon` | `K`, `n`, `r2` |
| `energy` | `resilience`, `toughness` |
| `relax` | Relaxation: `kind="relaxation"`, `F0`, `F_inf`, `dF`, `pct`, `tau`, `r2`, `significant`, `eps_hold`, `E0`, `E_inf`, `k_v`, `E_v`, `duration`, `t_start`, plus `fit`, `t_seg`, `y_seg`. Creep: `kind="creep"`, `d_start`, `d_inf`, `dd`, `tau`, `F_hold`, … |
| `k_tan_now` | Live tangent stiffness (N/m) over the last 3 s of loading |
| `fractured` | True once fracture has been detected |
| `regime` | Text from `classify_regime` |
| `eps_now`, `sig_now`, `F_now` | Latest strain, stress, force |

## 15. Where to change what

| I want to… | Change |
|---|---|
| Use another serial port or baud rate | `SERIAL_PORT`, `BAUD_RATE` (or `--port`, `--baud`) |
| Retune the controller | `Q_*`, `R_INPUT`, `SECONDARY_WEIGHT`; or the weights in `lqr_weights` |
| Change the default force limit or preload | `FORCE_LIMIT_N`, `PRELOAD_N` (and `force_limit` in the sketch) |
| Change the stage or motor constants | The plant constants at the top, **and** the same constants in the sketch |
| Change the elastic fit window | The two *E fit window* inputs, or the 10–40 % default in `fit_elastic` |
| Change the yield offset | The *Yield offset* input (default 0.2 %) |
| Make hold detection more or less sensitive | `StageGUI.RATE_TOL`, `RATE_WINDOW`, `HOLD_MIN_S` |
| Add a new calculated property | Write a function in section 5; call it in `run_analysis`; add a key to `res`; print it in `results_text_lines` |
| Add a new plot | Create it in a `_build_*_tab` method, store the curve as `self.c_…`, fill it in the matching `update_*` method, and clear it in `reset_data` |
| Add a new serial command | Add it to `handle()` in the sketch **and** to `SimulatedArduino._command`; send it from a `StageGUI` method with `send(...)` |
| Test with a different virtual sample | Edit `MaterialParams` defaults (simulator ground truth only) |
