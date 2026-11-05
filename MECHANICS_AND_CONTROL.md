# Mechanics and control of the straining stage with a force state

This document derives the model and the controller implemented in `ControlInterface4.py`. It also
describes how the GUI calculates the material properties. The wiring and firmware are in
`LOADCELL_WIRING.md`.

The only sample inputs are the cross-section $A_{cs}$ and the gauge length $L_0$. The displacement $d$
comes from the encoder and the force $F$ from the FUTEK LSB200 load cell. Everything else, including
$E$ and the stiffness $k$, is calculated from these four quantities.

## 1. State-space model with the force as a state

The previous model was $x = [\omega, I, d]^T$. In it, the sample was a spring, $F = \frac{EA}{L_0} d$,
so the force was not an independent quantity. Adding the load cell makes $F$ a measured state:

$$x = \begin{bmatrix} \omega & I & d & F \end{bmatrix}^T, \qquad u = V_{dc} D_{cycle}, \qquad g = \frac{p}{\pi}$$

$$
\begin{aligned}
J\dot\omega &= -b\,\omega + k_t I - g F\\
L\dot I &= -k_t\,\omega - R I + u\\
\dot d &= g\,\omega\\
\dot F &= k_i\, g\,\omega + \frac{k_r}{\tau}\, d - \frac{1}{\tau} F \quad(+\ \text{constant, removed by the integrator})
\end{aligned}
$$

$$
A=\begin{bmatrix}
-\frac{b}{J} & \frac{k_t}{J} & 0 & -\frac{g}{J}\\
-\frac{k_t}{L} & -\frac{R}{L} & 0 & 0\\
g & 0 & 0 & 0\\
k_i\,g & 0 & \frac{k_r}{\tau} & -\frac{1}{\tau}
\end{bmatrix},\quad
B=\begin{bmatrix}0\\ \frac{1}{L}\\ 0\\ 0\end{bmatrix},\quad
C=I_4,\quad D=0
$$

- $\omega$ and $d$ come from the encoder and $F$ from the load cell.
- $I$ is estimated as $(u - k_t\omega)/R$, which is valid because $L/R \approx 1.8\ \mu s$.
- Substituting $F = k d$ into the first row gives the report's term $-\frac{p}{\pi}\frac{EA}{L_0 J}$.
  The 3-state model is therefore the special case of a purely elastic sample.

The material enters through three **measured** parameters:

- $k_i$: instantaneous stiffness, the force response to a fast step in $d$;
- $k_r$: relaxed or tangent stiffness, the force change per unit $d$ once the material has relaxed;
- $\tau$: relaxation time.

## 2. Material mechanics: elastic and inelastic regimes

The sample is modelled as an elastoplastic backbone $F_{ep}(d)$ in parallel with a Maxwell arm
$F_v$:

$$F = F_{ep}(d) + F_v,\qquad \dot F_v = k_v\,\dot d - \frac{F_v}{\tau}$$

**Elastic.** $\sigma = E\varepsilon$, $F = k d$ and $k = \dfrac{E A_{cs}}{L_0}$. Unloading and reloading
below the previous maximum are also elastic, with the same $E$.

**Viscoelastic (Maxwell arm, standard linear solid).**
- The instantaneous modulus is $E_0 = E + E_v$ and the relaxed modulus is $E_\infty = E$.
- At constant displacement the force relaxes as

  $$F(t) = F_\infty + \Delta F\, e^{-(t-t_h)/\tau},\qquad \Delta F = F_0 - F_\infty\ \text{(the relaxation force)}$$

- After a ramp at rate $v$ lasting $t_r$, $\Delta F = k_v\, v\, \tau\,(1 - e^{-t_r/\tau})$, which gives
  $k_v$ and $E_v = k_v L_0 / A_{cs}$.
- At constant force the displacement creeps instead, with the retardation time
  $\tau_{creep} = \tau\,k_i/k_r$.

**Plastic.**
- Yield occurs when $\sigma > \sigma_Y(\bar\varepsilon_p)$, with isotropic hardening.
- The simulator's ground truth uses $\sigma_Y = \sigma_{y0} + Q(1-e^{-\beta\bar\varepsilon_p}) + H\bar\varepsilon_p$
  with a radial return.
- On the loading envelope the tangent modulus is
  $E_t = \dfrac{d\sigma}{d\varepsilon} = \dfrac{E\,H}{E+H}$ for linear hardening.
- The total strain splits as $\varepsilon = \sigma/E + \varepsilon_p$. The plastic part stays after
  unloading and is the permanent set.

**Necking and fracture.** Beyond the UTS the tangent stiffness becomes negative. At fracture the force
drops to zero.

**Linearised force row.** Around the operating point $F_{ep} \approx k_{tan} d + c$, so

$$\dot F = (k_{tan} + k_v)\,\dot d - \frac{F - F_{ep}(d)}{\tau}
\;\Rightarrow\; k_i = k_{tan} + k_v,\quad k_r = k_{tan}$$

| Regime | $k_i$ | $k_r$ | $\tau$ |
|---|---|---|---|
| Elastic loading, unloading or reloading | $k_e + k_v$ | $k_e$ | $\tau_v$ |
| Plastic loading (on the envelope) | $k_{tan} + k_v$ | $k_{tan} = E_t A_{cs}/L_0$ | $\tau_v$ |
| Hold (relaxation or creep) | as the preceding loading regime; $k_v$ and $\tau$ from the hold fit | | |
| Softening or necking ($k_{tan} \le 0$) | force control is ill-posed; the GUI keeps displacement control | | |

All entries are measured:
- $k_e$ is the slope of the $F$–$d$ fit;
- $k_{tan}$ is the slope of the last 3 s of the loading envelope, floored at $0.01\,k_e$;
- $k_v$ and $\tau$ come from the last significant relaxation fit.

Until a hold has been fitted, $k_v = 0$ and $\tau = 1$ s, where $\tau$ is only a numerical placeholder.

## 3. Reachability, stabilizability and observability

$W_r = [B,\ AB,\ A^2B,\ A^3B]$. $B$ and $AB$ are zero in the $d$ and $F$ rows, so $W_r$ is block
upper-triangular and

$$\det W_r = \det\begin{bmatrix}0 & \frac{k_t}{JL}\\ \frac{1}{L} & -\frac{R}{L^2}\end{bmatrix}\cdot
\det\begin{bmatrix}\frac{g k_t}{JL} & g\,a_3\\ \frac{k_i g k_t}{JL} & k_i g a_3 + \frac{g k_t}{JL\tau}(k_r-k_i)\end{bmatrix}
= \frac{g^2 k_t^3\,(k_i - k_r)}{J^3 L^4 \tau}$$

where $a_3 = -\frac{b k_t}{J^2 L} - \frac{k_t R}{J L^2}$. The first three rows of $W_r$ are the
3×3 matrix of the report.

- **Inelastic sample** ($k_i \neq k_r$): $W_r$ has rank 4, so the system is fully reachable.
- **Purely elastic sample** ($k_i = k_r = k$): $W_r$ has rank 3. The mode $F - k d$ obeys
  $\frac{d}{dt}(F - k d) = -(F - k d)/\tau$ whatever $u$ is, so it is uncontrollable, with
  $\lambda = -1/\tau$.

  Physically, the force and the displacement of a spring cannot be set independently. The mode is
  stable, so the system is **stabilizable** and LQR remains valid. The GUI's Dynamics tab shows this
  live: rank 3 before a relaxation has been measured, rank 4 once $k_v > 0$.
- **Numerical check.** A raw rank test is unreliable because $\operatorname{cond}(W_r) \sim 10^{20}$.
  The code uses the analytic determinant plus a PBH test on a state-scaled system, with scales
  $[1\ \text{rad/s},\ 0.1\ \text{A},\ 10^{-4}\ \text{m},\ 10^{-4} k_i\ \text{N}]$.
- **Observability** from the physical sensors $[d, F]$ alone has rank 4 in every regime. $d$ gives
  $\omega = \dot d / g$, and $\dot\omega$ then gives $I$. The firmware's derived $\omega$ and estimated
  $I$ are therefore justified.

## 4. Controller

**LQR with integral action on the tracked output** $y = C_{tr}x$: $y = d$ in displacement mode and
$y = F$ in force mode.

$$A_{aug} = \begin{bmatrix}A & 0\\ -C_{tr} & 0\end{bmatrix},\quad B_{aug} = \begin{bmatrix}B\\0\end{bmatrix},\quad
K = \operatorname{lqr}(A_{aug}, B_{aug}, Q, R_u),\quad R_u = 10^{-3}$$

**Weights,** for the state order $[\omega, I, d, F, z]$:

| Mode | Weights | Notes |
|---|---|---|
| Displacement | $\operatorname{diag}(10^{-3},\,10^{-3},\,5\cdot10^8,\,10^{-3}\cdot5\cdot10^8/k_e^2,\,3\cdot10^8)$ | The existing tuning, plus a small force weight scaled by the measured elastic stiffness. |
| Force | $\operatorname{diag}(10^{-3},\,10^{-3},\,5\cdot10^5,\,5\cdot10^8/k_i^2,\,3\cdot10^8/k_i^2)$ | Weights scaled by $1/k_i^2$, because the plant gain from $d$ to $F$ is $k_i$. |

The scaling holds the force-loop poles at about $-12$ rad/s in every regime:

| $k_i$ (N/m) | $k_r$ (N/m) | Dominant closed-loop poles, force mode (rad/s) |
|---|---|---|
| 1000 | 1000 | −0.5, −0.78, −12.0 |
| 1000 | 200 | −0.1, −0.78, −12.0 |
| 1000 | 50 | −0.025, −0.78, −12.0 |

The slowest pole is the physical creep mode, $-k_r/(\tau k_i)$. With fixed weights the force loop would
become 20–50× slower after yield. That is why the gains are scheduled.

**Time-scale separation.** The electrical and mechanical poles are at about $-2.8\cdot10^5 \pm
4.5\cdot10^5 j$ rad/s. The closed-loop dynamics that matter are at most 12 rad/s. A continuous-time
design implemented at 1 kHz is therefore valid, and so is the 80 samples/s load cell (500 rad/s).

**Stiffness-free start.** Nothing is known about the sample at first, so the first design uses
$k_i = k_r = 0$. This gives

$$K_0 = [0.0385,\ 0.156,\ 7.51\cdot10^5,\ \approx 0,\ -5.48\cdot10^5]$$

Applied to true plants from $k = 10$ to $10^6$ N/m, and to elastic, viscoelastic and plastic samples,
$K_0$ keeps the same dominant poles ($-0.78$ and $-12$ rad/s). The reason is that the drive is almost
a velocity source: $\partial\omega/\partial F = -Rg/k_t^2 \approx -7\cdot10^{-5}$ rad/s per N.

Displacement control is therefore safe from the first second. $K_0$ is compiled into the firmware, so
jog and preload work before any gains are sent. **Force mode stays locked until $k$ has been
measured**, because its gains scale with $k_i$.

**Control law** in the firmware and the simulator:

$$u = -\left(k_\omega\omega + k_I I + k_d\,\tilde d + k_F\,\tilde F + k_{int} z\right),\qquad \dot z = r - y$$

Here $\tilde d = d - r$ and $\tilde F = F$ in displacement mode, and $\tilde d = d$ and $\tilde F = F - r$
in force mode. Feeding the reference to the proportional term as well gives fast ramp tracking. The
integrator still removes any offset, for example from gearbox friction.
- Anti-windup is by conditional integration, and $u$ is saturated at $\pm V_{dc}$.
- When the gains, mode or model change, the integrator is reset to
  $z = -(u_{now} + K_x\tilde x)/k_{int}$, so $u$ does not jump (bumpless transfer).

**Gain scheduling.** The GUI recomputes the LQR when:
- $k_i$, $k_r$ or $\tau$ changes by more than 25 %,
- the model becomes measured, or
- softening starts.

It resends the gains at most once per second, and only while the controller runs. If softening
starts in force mode, it switches to a displacement hold.

## 5. How the GUI calculates the properties

All of these are calculated automatically, every 0.5 s, from $d$, $F$, $A_{cs}$ and $L_0$.

**Strain and stress.**
- $\varepsilon = (d - d_0)/L_0$, where $d_0$ is set when the preload completes; $\sigma = F/A_{cs}$.
- True values: $\varepsilon_t = \ln(1+\varepsilon)$ and $\sigma_t = \sigma(1+\varepsilon)$.

**Branches.** Each sample is labelled loading (on the envelope), reloading, unloading or hold. The label
comes from the displacement rate over 0.5 s, with a threshold of at least 4 encoder counts per window.

**Young's modulus and stiffness.**
- A linear fit of $\sigma$–$\varepsilon$ and of $F$–$d$ on the virgin loading branch, so $E = k L_0/A_{cs}$.
- The default window is $\sigma$ between 10 % and 40 % of $\sigma_{max}$, on the first rise. A window
  can be set in mm instead.
- On a viscoelastic sample a ramp measures an apparent modulus between $E_\infty$ and $E_0$, which
  depends on the rate.

**Yield and the elastic limit.**
- Yield uses the 0.2 % offset method (configurable): the first crossing with
  $\sigma = E(\varepsilon - 0.002) + \sigma_0$.
- The proportional limit is the first sustained deviation from the elastic line.

**Ultimate strength and fracture.**
- UTS is the maximum $\sigma$; the strain at UTS is the uniform elongation.
- Fracture is detected when, after the UTS, $\sigma$ falls below 10 % of the UTS at a larger strain.

**Plastic regime.**
- The tangent modulus along the envelope comes from a sliding regression.
- $H = d\sigma/d\varepsilon_p$, with $\varepsilon_p = \varepsilon - (\sigma-\sigma_0)/E$, and the mean $E_t$
  are taken between yield and UTS.
- Hollomon $\sigma_t = K\varepsilon_t^n$ is fitted in log–log.
- The permanent set is the strain after *Unload*.

**Energies.**
- Resilience is $\sigma_y^2/(2E)$.
- Toughness is the area under the loading curve.

**Relaxation.**
- It is fitted on the last hold of at least 2 s that follows loading. A hold after unloading shows
  recovery, not relaxation.
- The fit is $F(t) = F_\infty + \Delta F e^{-(t-t_h)/\tau}$, by least squares, extrapolated back to the
  start of the hold.
- Reported: $F_0$, $F_\infty$, the **relaxation force $\Delta F$** and its percentage, $\tau$, the secant
  moduli $E_0$ and $E_\infty$ at the hold strain, and $k_v$ and $E_v$.
- A fit is used for the model only when $|\Delta F| > 3\times$ the residual noise.

**Creep.** In a force-mode hold, $d(t)$ is fitted with the same exponential.

**Live regime.** Slack, elastic, plastic, unloading, hold, softening or fractured.

**Model prediction.** A material model built only from the identified values (E, yield, H, E_v and τ)
is driven by the measured $d(t)$. Its force is overlaid on the measured force in the *Signals* tab.

## 6. Verification

`python ControlInterface4.py --selftest` runs 15 checks without hardware, and all pass:

- **Model and controller:**
  - the analytic $\det W_r$ equals the numeric one;
  - the elastic case has exactly one uncontrollable, stable mode at $-1/\tau$;
  - LQR is stable in all regimes and both modes, and force mode keeps a pole near −12 rad/s;
  - $K_0$ is stable from $k = 10$ to $10^6$ N/m;
  - force mode is refused without a measured $k$;
  - the system is observable from $[d, F]$.
- **Analysis:** on synthetic elastoplastic data, E, k, the 0.2 % yield, H and $E_t$ are recovered
  within 2 %. On a noisy relaxation curve, τ and ΔF are recovered within 2 %.
- **Closed loop (simulated stage):**
  - displacement and force steps settle on the reference;
  - preload stops at the set force;
  - `FLIM` cuts the motor;
  - `JOG` moves the requested distance;
  - `RAMP` runs at the commanded rate within 2 %.

`python ControlInterface4.py --sim --demo` runs this sequence through the real GUI buttons:
jog → tare → preload → ramp to 0.5 mm at 0.01 mm/s → 30 s hold → unload → ramp to fracture.
The simulated sample is a viscoelastic-plastic polymer film (0.01 mm², 10 mm).

| Property (calculated from d, F, A_cs, L0) | Demo result | Simulator ground truth |
|---|---|---|
| E (apparent, at 0.1 %/s) | 1.08 GPa | backbone 1.0 GPa, instantaneous 1.3 GPa |
| Stiffness k | 1080 N/m | backbone 1000 N/m |
| Yield (0.2 % offset) | 32.9 MPa at 3.10 % | σ_y0 = 30 MPa, plus hardening and viscous overstress |
| UTS | 52.5 MPa at 12.0 % | necking starts at 12 % |
| Fracture strain | 14.8 % | 15 % |
| Relaxation force ΔF (30 s hold at 5 %) | 0.0246 N | 0.0240 N expected |
| Relaxation time τ | 7.6 s | 8 s |
| Maxwell modulus E_v | 0.32 GPa | 0.30 GPa |
| Permanent set after unloading | 1.9 % | (plastic strain at 5 %) |

These are simulation results. They show that the analysis and the control chain work end to end,
not how accurate the real stage is. Real data will add gearbox friction, alignment effects and the
compliance of the rig.
