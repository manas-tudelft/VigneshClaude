# How to run the straining-stage files

This guide covers every file in the repository: what it is, how to install and run it, and in what order to
bring up the real stage. For the *why* of the model, see [MECHANICS_AND_CONTROL.md](MECHANICS_AND_CONTROL.md).
For the wiring and the electronics, see [LOADCELL_WIRING.md](LOADCELL_WIRING.md). For every function in the
Python script, see [CODE_GUIDE.md](CODE_GUIDE.md).

> **No YAML files.** The repository contains no `.yaml`/`.yml` files and none are needed. The only
> configuration files are `requirements.txt` (Python packages) and the constants at the top of
> `ControlInterface4.py` and `StrainingStage.ino`.

## 1. What is in the repository

| File | What it is | Run / open it with |
|---|---|---|
| `ControlInterface4.py` | The whole PC side in one file: GUI, plant model, LQR design, material analysis, simulated stage, demo and self-test. | `python ControlInterface4.py` |
| `firmware/StrainingStage/StrainingStage.ino` | Arduino Uno sketch (serial protocol v2): 1 kHz control loop, encoder, motor PWM, HX711 load cell. | Arduino IDE or `arduino-cli` |
| `requirements.txt` | Python packages the script needs. | `pip install -r requirements.txt` |
| `LOADCELL_WIRING.md` | Parts list, wiring tables, grounding, calibration, protocol tables. | Read |
| `docs/wiring_diagram_loadcell.png` / `.svg` | Your original wiring diagram with the LSB200 and HX711 added in colour. | Any image viewer |
| `MECHANICS_AND_CONTROL.md` | Model, regimes, reachability, LQR, property formulas, verification. | Read (GitHub renders the maths) |
| `CODE_GUIDE.md` | Reference of every class, function and method in `ControlInterface4.py`. | Read |
| `WIRING_CHANGES.md` | **Superseded** (INA226 current-sensing plan). Kept for reference only. | Do not follow |
| `Screenshot ….png`, `Prompt.txt` | Original input material. No code uses them. | – |

## 2. Install the PC side

1. Install Python 3 (64-bit) and, optionally, create a virtual environment:

   ```
   python -m venv .venv
   .venv\Scripts\activate          # Windows
   ```

2. Install the packages:

   ```
   pip install -r requirements.txt
   ```

   This installs `numpy`, `scipy`, `control`, `PyQt5`, `pyqtgraph` and `pyserial`.
   - `--selftest` needs only `numpy`, `scipy` and `control`. Qt is imported inside a `try`, so the self-test
     runs without it.
   - `pyserial` is imported only when a real serial port is opened, so `--sim` works without it.

## 3. Run the GUI

All commands are run from the repository folder.

| Command | What happens |
|---|---|
| `python ControlInterface4.py` | Opens the GUI on the real stage on `COM9` at 115200 baud. |
| `python ControlInterface4.py --port COM5` | Same, on another serial port. |
| `python ControlInterface4.py --port COM5 --baud 115200` | Also sets the baud rate (must match `Serial.begin` in the sketch). |
| `python ControlInterface4.py --sim` | Simulated stage and sample. No hardware, no `pyserial`. Real-time speed. |
| `python ControlInterface4.py --sim --sim-speed 10` | Simulation running 10× faster than real time. |
| `python ControlInterface4.py --sim --demo` | Scripted test through the real buttons (defaults to 20× speed), prints the results next to the simulator's ground truth, then closes. |
| `python ControlInterface4.py --sim --demo --screenshot out` | The same, and also saves `demo_signals.png`, `demo_material.png`, `demo_dynamics.png` and `demo_data.csv` into the folder `out`. |
| `python ControlInterface4.py --selftest` | Runs the 15 built-in checks without a window or a port and exits with status 0 (all passed) or 1. |

Notes:
- To change the default port permanently, edit `SERIAL_PORT = 'COM9'` near the top of the script.
- If the serial port cannot be opened, the GUI asks whether to start the simulated stage instead.
- After the port opens the script waits 2 s, because the Uno resets when the port is opened.
- **Only one program can hold the port.** Close the Arduino IDE Serial Monitor before starting the GUI.

### Try it first in simulation

```
python ControlInterface4.py --selftest
python ControlInterface4.py --sim --demo
```

The demo takes about 20 s of wall time (it runs the simulator at 20× speed). It runs: jog → tare → pre-tension → ramp to 0.5 mm with a
30 s hold → unload → ramp to fracture. The simulated sample has known properties (E ≈ 1 GPa, relaxation
time 8 s, fracture at 15 % strain), so you can see how close the calculated values come.

## 4. Upload the Arduino sketch

### 4.1 Prepare

1. Install the **Arduino IDE** (or `arduino-cli`).
2. Install the library **"HX711 Arduino Library" by Bogdan Necula** (Library Manager). `EEPROM` ships with the
   core.
3. Open `firmware/StrainingStage/StrainingStage.ino`. The folder name must stay equal to the sketch name.
4. Edit the constants marked `TODO` near the top of the sketch:

   | Constant | What to put there |
   |---|---|
   | `COUNTS_PER_REV` | Encoder counts per **output-shaft** revolution after ×4 decoding. Copy it from your old sketch. |
   | `ENC_SIGN` | `+1` or `-1`, so that moving *outwards* gives positive d. |
   | `MOTOR_SIGN` | `+1` or `-1`, so that a positive voltage u moves *outwards*. |
   | `LC_CAPACITY_N` | Rated capacity of your LSB200 in newtons (from its calibration certificate). |
   | `LC_RATED_OUTPUT_MV_V` | Rated output in mV/V (from the certificate). |
   | `OMEGA_TAU` and the `pinMode` lines | Keep these in line with your old sketch if it had a better speed filter or a different encoder input mode. |

   The plant constants (`P_LEAD`, `K_T`, `R_MOTOR`, `V_DC`) are the same as in the GUI. Do not change one
   side without the other.

### 4.2 Upload

**Arduino IDE:** Tools → Board → *Arduino Uno*, Tools → Port → your port, then Upload.

**arduino-cli:**

```
arduino-cli lib install "HX711 Arduino Library"
arduino-cli compile --fqbn arduino:avr:uno firmware/StrainingStage
arduino-cli upload  --fqbn arduino:avr:uno -p COM5 firmware/StrainingStage
```

> **Status of the sketch:** it was written to mirror the simulator inside `ControlInterface4.py`, which is
> tested, but the sketch itself has **not been compiled or run on hardware**. Expect to fix small compile
> errors or constants on the first build.

### 4.3 Bench test with the Serial Monitor (before the GUI)

Open the Serial Monitor at **115200 baud** with line ending **"Newline"**.

1. With the **24 V motor supply switched off**, you should see the banner
   `Straining stage FW v2 (load cell, protocol v2)`, followed by data lines
   `omega,I,d,F,u` about 80 times per second. Data lines are sent once per HX711 sample, so if you see the
   banner but **no data lines**, the HX711 is not wired (or not powered): check DT on D4, SCK on D5, VCC and GND.
2. Press on the load cell by hand. `F` (the fourth number) must change. If tension reads negative, swap the
   green and white signal wires (see [LOADCELL_WIRING.md](LOADCELL_WIRING.md)).
3. Turn the motor shaft by hand. `d` (the third number) must change, positive when moving outwards. If not,
   flip `ENC_SIGN`.
4. Send `TARE` (cell unloaded). The sketch replies `TARE done, offset = …`, and `F` goes to about 0.
5. Only then switch on the motor supply. Send `FLIM:5` and, with the sample out, try
   `JOG:1,0.0001,0.00005` (outwards 0.1 mm at 0.05 mm/s). The motor must move outwards and the sketch
   replies `JOG done`. If it moves inwards, flip `MOTOR_SIGN`.

Close the Serial Monitor before starting the GUI.

### 4.4 Serial commands (PC → Arduino, one per line)

| Command | Meaning |
|---|---|
| `r,k_omega,k_I,k_d,k_F,k_int,mode` | Gains and reference (mode 0: track d in m, mode 1: track F in N). The GUI sends these for you. |
| `R:<r>` | New reference. |
| `RAMP:<target>,<rate>` | Ramp the reference (m and m/s, or N and N/s); replies `RAMP done`. Controller must be running. |
| `JOG:<dir>,<distance m>,<speed m/s>` | Distance-based move; replies `JOG done`. Distance 0 = continuous until `STOP`. |
| `PRELOAD:<N>` | Move outwards until F ≥ N; replies `Pre-tensioning complete, controller stopped`. |
| `START`, `STOP`, `RESET` | Start the controller, stop the motor, zero the encoder and integrator. |
| `TARE`, `CAL:<known N>`, `FLIM:<N>` | Load-cell zero, scale factor (both saved in EEPROM), and force limit. |
| `MANUAL:<dir>,<magnitude %>` | Legacy open-loop duty. The GUI no longer uses it. |

Arduino → PC: one data line `omega,I,d,F,u` per HX711 sample (rad/s, A estimated, m, N, V), and text lines
such as `RAMP done`, `FORCE LIMIT: …`, `DISPLACEMENT LIMIT, motor stopped`, `GAINS OK mode=…`.

## 5. Wire it

Follow [LOADCELL_WIRING.md](LOADCELL_WIRING.md) and `docs/wiring_diagram_loadcell.png`. In short:

| From | To |
|---|---|
| LSB200 red (+Excitation) / black (−Excitation) | HX711 E+ / E− |
| LSB200 green (+Signal) / white (−Signal) | HX711 A+ / A− (swap if tension reads negative) |
| HX711 VCC / GND | Arduino 5V / GND |
| HX711 DT / SCK | Arduino **D4** / **D5** |
| HX711 RATE | VCC, for 80 samples per second |

Check the lead colours against the FUTEK certificate of your own cell. Everything else in the original
diagram stays as it is. Fit a **mechanical overload stop**: the 5 N software cut-off reacts within 12.5 ms,
but the drive can push much harder.

## 6. A complete test on the real stage

1. Power up, start the GUI, and check the label *Firmware protocol: v2*.
2. With the sample slack, press **Tare**. Check the force limit (default 5 N) and press **Set force limit**.
3. Enter the two sample inputs: **A_cs** (mm²) and **L0** (mm). Nothing else about the material is entered.
4. *(Once per set-up)* Hang a known mass on the cell, enter it in grams and press **Calibrate with hanging
   mass**.
5. **Manual control (distance-based):** enter a distance (mm) and a speed (mm/s), then press *Outwards* or
   *Inwards* to take up the slack. A distance of 0 moves continuously until *Stop*.
6. **Pre-tension (preload, sets strain zero):** the stage moves outwards until the preload force (default
   0.01 N) is reached. That position becomes d = 0 and ε = 0.
7. **Ramp → Hold:** enter the target d (mm from the preload zero), the rate (mm/s) and the hold time (s).
   - The stress–strain and force–displacement curves build up live in the *Material* tab.
   - E, k, yield, UTS and hardening appear in *Results* as soon as there is enough data.
   - During the hold, the relaxation force ΔF and time constant τ are fitted.
8. **Unload:** ramps back down to the preload force and reports the permanent set.
9. **Force (stress) mode** unlocks only after k has been measured. Choose *Force*, enter a stress in MPa,
   press **Send Parameters to Arduino**, then **Start Controller** (or **Update Reference**).
10. **Save Data to File** writes a CSV. **Reset Measurement** clears everything for the next sample.

The *Auto gain scheduling* checkbox (on by default) recomputes and resends the gains while the controller is
running, whenever the measured stiffness changes by more than 25 % (at most once per second).

### What the three tabs show

- **Signals:** d and its reference, F (load cell, the identified-model prediction and an optional drive-train
  estimate), ω, and the drive voltage with the estimated current.
- **Material:** stress–strain with the E fit, the 0.2 % offset line, and yield/UTS/fracture markers; the
  force–displacement curve with the stiffness line; the tangent modulus; and the last hold with its
  exponential fit.
- **Dynamics:** the current A, B, C, D, gains K, poles, reachability and observability; the measured k_i and
  k_r over time; and the pole map.

### The CSV file

Comment lines (starting with `#`) hold A_cs, L0, the strain origin, the calculated properties and the event
log. The data columns are:

`Time(s), Displacement(m), Elongation(m), Strain(%), Force(N), Stress(Pa), TrueStrain(-), TrueStress(Pa),
Omega(rad/s), Current_est(A), Voltage(V), Reference(m or N), Mode(0=disp,1=force,-1=off)`

To load it in Python (skips the comment, blank and header lines), or just open it in Excel:

```python
import io
import numpy as np

with open("measurement.csv") as f:
    rows = [ln for ln in f if ln.strip() and not ln.startswith(("#", "Time"))]
data = np.loadtxt(io.StringIO("".join(rows)), delimiter=",")   # shape (samples, 13)
t, d, strain_pct, F, stress = data[:, 0], data[:, 1], data[:, 3], data[:, 4], data[:, 5]
```

## 7. Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| "Could not open COM9" | Wrong port or another program holds it (Arduino Serial Monitor). Use `--port`, close the other program, or accept the simulation prompt. |
| Warning "No protocol-v2 data line has been received yet" | The old sketch is still on the Uno, or the baud rate differs. Upload `StrainingStage.ino`; the baud must be 115200. |
| Label says *LEGACY (omega,d)* | Old firmware: force features are disabled. Upload the new sketch. |
| Force reads negative when stretching | Swap the HX711 A+ and A− wires (green and white). |
| *Outwards* gives negative d | Flip `ENC_SIGN` in the sketch. |
| Motor runs the wrong way for positive voltage | Flip `MOTOR_SIGN` in the sketch. |
| Force mode greyed out | By design: it unlocks once k has been measured. Load the sample in displacement control first. |
| "FORCE LIMIT … motor stopped" | The load exceeded the limit (default 5 N). Reduce the speed or target; do not raise the limit above the cell's rated capacity. |
| "DISPLACEMENT LIMIT" | Travel exceeded 20 mm (`D_LIMIT`). Press *Reset* and re-zero. |
| Readings jump when the motor runs | Grounding or cabling: see section 6 of [LOADCELL_WIRING.md](LOADCELL_WIRING.md). |
| Force-mode warning about softening | The sample is necking: force control is ill-posed. The GUI switches to displacement control. |
| `ModuleNotFoundError` | `pip install -r requirements.txt` in the active environment. |
