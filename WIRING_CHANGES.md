# Wiring and firmware changes for motor current sensing

`ControlInterface4.py` now uses the 3-state model x = [ω, I, d]. It estimates
Young's modulus from the motor current, so the motor current has to be measured
and sent to the PC.

## 1. Current sensor (recommended: INA226 breakout, 0.1 Ω shunt)

The force on the sample is F = (π/p)·(k_t·I − b·ω), which is about **58.6 kN per A**.
A thin sample therefore needs very fine current resolution:

| Sensor                     | Current resolution on Uno | Force resolution |
|----------------------------|---------------------------|------------------|
| ACS712-05B (analog, A0)    | ~26 mA                    | ~1.5 kN          |
| INA226 + 0.1 Ω (I²C, 16 bit) | ~25 µA                  | ~1.5 N           |
| INA226 + 1 Ω (I²C)         | ~2.5 µA (max 82 mA)       | ~0.15 N          |

Use the INA226. A Hall sensor such as the ACS712 is only usable with very stiff or thick samples.
The INA226 is bidirectional, so it handles both directions of the H-bridge.
Its common-mode range is 0–36 V, so it can sit on the 24 V motor leg.
If you use a 1 Ω shunt, add it to `R` in the code (R = 55.041 + 1).

## 2. Changes to the wiring diagram

All existing connections stay unchanged: encoder → optocoupler → D2/D3,
driver IN3/IN4 → D9/D10, buck converter, and bench supply.
Only the following is added:

| From                          | To                        | Purpose                     |
|-------------------------------|---------------------------|-----------------------------|
| Motor driver **OUT3**         | INA226 **IN+** (VIN+)     | Shunt in series with motor  |
| INA226 **IN−** (VIN−)         | Motor terminal (was OUT3) | Shunt in series with motor  |
| INA226 **VCC**                | Arduino **5V**            | Sensor logic supply         |
| INA226 **GND**                | Arduino **GND** (common ground rail) | Common reference |
| INA226 **SDA**                | Arduino **A4**            | I²C data                    |
| INA226 **SCL**                | Arduino **A5**            | I²C clock                   |
| INA226 **A0, A1**             | GND                       | I²C address 0x40            |

```
                 Motor Driver                    INA226                 Motor
  OUT3 ───────────────────────────────► IN+ [shunt] IN− ───────────────► M (+)
  OUT4 ─────────────────────────────────────────────────────────────────► M (−)

  Arduino 5V  ─── VCC      Arduino A4 ─── SDA
  Arduino GND ─── GND      Arduino A5 ─── SCL
```

In other words, cut the blue wire between OUT3 and the motor and insert the shunt
(IN+ on the driver side). If the measured force comes out negative while
the sample is being stretched ("Outwards"), swap IN+ and IN−, or negate the
current in the firmware.

## 3. Arduino firmware changes (the sketch is not in this repo)

1. **Read the current** each control cycle, for example with the `INA226_WE` or
   `INA226` (Rob Tillaart) library. Set the averaging to 16 or more samples so the PWM ripple
   is averaged out.
2. **Serial output** must become `omega,I,d` or `omega,I,d,V`
   (I in A, V = the applied voltage u = V_dc·D_cycle in V). The GUI accepts both formats.
   If V is omitted, it is reconstructed as V = R·I + k_t·ω.
3. **Gain message** now has four gains: `r,k_omega,k_I,k_d,k_int`
   (previously `r,k_omega,k_d,k_int`). The control law becomes
   `u = -(k_omega*omega + k_I*I + k_d*d + k_int*z)`, using the same integrator
   `z` and the same sign convention as before. Only the `k_I*I` term is new.

Note: the electrical time constant L/R ≈ 1.8 µs is much shorter than the Arduino
loop period, so the current behaves as a quasi-static state. The LQR gains for ω, d and the
integrator are practically identical to the old 2-state design, and k_I is small.

## 4. How E is calculated in the GUI

1. Input power: P_in = V·I. Mechanical power: P_mech = P_in − I²R (= τ·ω).
2. Motor torque: τ_m = k_t·I. This is used for the force rather than P_mech/ω, because
   that ratio is singular at standstill, where a quasi-static test spends most of its time.
3. Sample force: F = (π/p)·(τ_m − b·ω). The inertia term J·ω̇ is neglected.
4. σ = F/A_cs and ε = d/L0. E is the slope of a linear fit of σ against ε.
   The strain window for the fit can be set in the GUI.

Unmodelled Coulomb friction in the gearbox and lead screw shows up as a force offset
and as hysteresis. The fit intercept absorbs the offset. For best results, fit
only on the loading branch.

_____ Claude output text

I've updated [ControlInterface4.py](ControlInterface4.py) to calculate Young's modulus from the measured motor current, and to plot stress–strain and force–displacement curves in the GUI. The wiring and firmware changes are in this file (sections 1–4 above). The control and E-calculation math works in a test with made-up data, which recovered the E it was given exactly. I couldn't open the GUI because PyQt5 wasn't installed on the machine where the code was written, so the window itself is untested.

**Main issue: the current is tiny for thin samples.** Because of the gearbox and lead screw, 1 A of motor current equals about 58,600 N on the sample. For the GUI's default sample (0.01 mm², 10 mm long, about 1–2 GPa), 1 % strain needs only about 3.4 µA. A Hall sensor like the ACS712 can't measure that on an Arduino, and even an INA226 can only resolve about 1.5 N. Gearbox friction will also affect the result. This method works much better with stiffer or thicker samples. A larger 1 Ω shunt with the INA226 improves resolution about ten times but limits it to about 82 mA.

**Code changes**
- **System model:** uses the 3×3 A, B, C, D matrices from the report, with state [ω, I, d]. The controller (LQR) now produces four gains: ω, I, d and the integral term.
- **Reachability check:** I check the determinant of the reachability matrix from the report instead of a numerical rank. The numerical rank wrongly reports rank 2 because the matrix is so badly scaled.
- **New gains:** the gains for ω, d and the integral term are almost the same as before, plus a small current gain. That's expected, because the current settles in about 1.8 µs, far faster than the Arduino loop.
- **Calculating E**, in this order:
  1. Input power is voltage × current, and mechanical power is that minus the I²R loss.
  2. Motor torque is k_t × current.
  3. The force on the sample comes from the torque, minus the bearing friction (b·ω).
  4. Stress is force ÷ area and strain is displacement ÷ length.
  5. E is the slope of a straight-line fit of stress against strain.
- **Why torque and not power:** I get the force from torque rather than from power ÷ speed, because that division blows up when the motor is nearly stopped, which is most of a slow tensile test. The two give the same answer when the motor is moving. Power is still calculated, plotted and saved.
- **GUI:** there are six plots: displacement, ω, current, power, stress–strain (with the fitted line), and force–displacement. A label shows the measured E and how well the line fits.
- **New controls:** you can set a strain limit for the fit. A checkbox, "Use measured E for controller gains", is on by default. The old E input is now an "initial guess", used only until E has been measured.
- **Saved file:** the CSV now also includes current, voltage, power, torque, force, stress and strain.

**Wiring changes** (everything already in the diagram stays the same):
- Cut the blue wire between motor driver OUT3 and the motor, and put the INA226 shunt in that gap: OUT3 to IN+, IN− to the motor.
- INA226 VCC to Arduino 5V, GND to the common ground, SDA to A4, SCL to A5. A4 and A5 are free on the current board.

**You'll need to change the Arduino sketch** (it isn't in this repo):
- Read the current every loop.
- Send `omega,I,d`, or `omega,I,d,V` with V the applied voltage.
- Read four gains in the form `r,k_omega,k_I,k_d,k_int`, and add the `k_I*I` term to the control law.

The old sketch's two-value serial lines won't be read as data; the GUI will just print them as Arduino messages.
