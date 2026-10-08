# Authors: Maxim Beckers & Twan ten Have

import sys
import serial
import numpy as np
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import time
import control as ct

# -----------------------
# System Parameters
# -----------------------
b = 1e-5
J = 6.20646e-6
k_t = 243.8 * 53.601e-3
L = 1e-4
R = 55.041
p = 0.0007
V_dc = 24.0
k_spring = 0.0  # Will be updated when parameters are set
A_cs_stored = 0.0    # Will be updated when parameters are set
E_stored = 0.0     # E used for the controller design (measured if available, else initial guess)
L0_stored = 0.0    # Will be updated when parameters are set
E_measured = None  # Young's modulus estimated from motor torque/power (Pa)
E_fit_offset = 0.0 # Stress intercept of the fit (Pa), absorbs sensor offset/pretension
E_fit_r2 = 0.0     # Coefficient of determination of the fit
pretension_value = 0.0
running_state = False



# -----------------------
# Serial setup
# -----------------------
ser = serial.Serial('COM9', 115200, timeout=1)
time.sleep(2)

def send_parameters(mode, ref_strain, ref_stress, A_cs, L0, E_guess):
    """Calculate controller gains from system dynamics and send to Arduino.
    E_guess is only used until E has been measured from the motor torque."""

    E = E_measured if (use_measured_checkbox.isChecked() and E_measured) else E_guess

    # Calculate r based on mode
    if mode == "Strain": # Strain mode
        r = ref_strain * L0
    else:  # Stress mode
        r = ref_stress * L0 / E
    
    A, B, C, D = calculate_system_dynamics(A_cs, L0, E)
    
    global A_cs_stored
    global L0_stored
    global E_stored
    
    A_cs_stored = A_cs  # Store for stress calculations
    L0_stored = L0    # Store for stress calculations
    E_stored = E      # Store for stress calculations
    
    A_aug = np.block([
        [A, np.zeros((3, 1))],
        [-C, np.zeros((1, 1))]
    ])

    B_aug = np.vstack((B, [[0]]))

    # Calculate controller gains, state order: [omega, I, d, integral of error]
    Q = np.diag([1e-3, 1e-3, 5e8, 3e8])
    Q_u = np.array([[1e-3]])

    K = ct.lqr(A_aug, B_aug, Q, Q_u)[0]
    
    # Format gains for Arduino as comma-separated values
    gains_str = ",".join([f"{gain:.6f}" for gain in K.flatten()])
    r_str = f"{r:.9f}"
    
    # Send to Arduino: "r,k_omega,k_I,k_d,k_int"
    command = f"{r_str},{gains_str}\n"
    ser.write(command.encode())
    print(f"Sent: r={r}, K={K.flatten()}, E used={E/1e9:.4f} GPa")

def calculate_system_dynamics(A_cs, L0, E):
    """
    Calculate system matrices from physical parameters.
    """

    global k_spring
    k_spring = E * A_cs / L0


    # -----------------------
    # Plant, state x = [omega, I, d], input u = V_dc * D_cycle
    # -----------------------
    A = np.array([
        [-b/J,       k_t/J,  -(p/np.pi) * (k_spring/J)],
        [-k_t/L,     -R/L,   0],
        [p/np.pi,    0,      0]
    ])

    B = np.array([
        [0],
        [1/L],
        [0]
    ])

    C = np.array([
        [0, 0, 1]
    ])

    D = np.array([[0]])

    # Reachability: W_r = [B, AB, A^2 B]. A numerical rank test is unreliable
    # here (cond(W_r) ~ 1e20 because 1/L dominates), so use the analytic
    # determinant det(W_r) = -p*k_t^2 / (pi*J^2*L^3).
    det_Wr = -p * k_t**2 / (np.pi * J**2 * L**3)
    if det_Wr == 0:
        raise ValueError("System is not reachable, cannot design controller")

    return A, B, C, D

def update_reference(mode, ref_strain, ref_stress):
    """Send new reference value to Arduino without recalculating gains"""
    # Calculate r based on mode
    if mode == "Strain": # Strain mode
        new_r = ref_strain * L0_stored
    else:  # Stress mode
        new_r = ref_stress * L0_stored / E_stored
    
    command = f"R:{new_r:.9f}\n"
    ser.write(command.encode())
    print(f"Reference updated to: {new_r}")

def motor_quantities(omega, current, voltage):
    """Motor input power, mechanical power and produced torque.
    voltage may be NaN when the Arduino does not report it; it is then
    reconstructed from the electrical equation V = R*I + k_t*omega (L*dI/dt ~ 0)."""
    omega = np.asarray(omega)
    current = np.asarray(current)
    voltage = np.where(np.isnan(voltage), R * current + k_t * omega, voltage)
    P_in = voltage * current          # electrical input power (W)
    P_mech = P_in - R * current**2    # power converted to mechanical (W) = tau*omega
    tau_m = k_t * current             # torque produced by the motor (Nm)
    return P_in, P_mech, tau_m

def sample_force(omega, tau_m):
    """Axial force on the sample from J*domega = tau_m - b*omega - (p/pi)*F.
    Quasi-static test: J*domega is neglected (J ~ 6e-6 kg m^2)."""
    return (np.pi / p) * (np.asarray(tau_m) - b * np.asarray(omega))

def estimate_youngs_modulus(strain, stress):
    """Fit stress = E*strain + offset over the selected strain window (strain as fraction)."""
    global E_measured, E_fit_offset, E_fit_r2
    max_strain = fit_limit_input.value() * 1e-2
    mask = strain <= max_strain if max_strain > 0 else np.ones_like(strain, dtype=bool)
    eps, sig = strain[mask], stress[mask]
    if len(eps) < 20 or np.ptp(eps) <= 0:
        return
    slope, offset = np.polyfit(eps, sig, 1)
    residual = sig - (slope * eps + offset)
    ss_tot = np.sum((sig - sig.mean())**2)
    E_measured = slope if slope > 0 else None
    E_fit_offset = offset
    E_fit_r2 = 1 - np.sum(residual**2) / ss_tot if ss_tot > 0 else 0.0

def send_manual_control(direction, magnitude):
    """Send manual control command to Arduino
    direction: 1 for clockwise, -1 for counterclockwise
    magnitude: 0-100 (percentage of max control input)
    """
    command = f"MANUAL:{direction},{magnitude:.2f}\n"
    ser.write(command.encode())
    print(f"Manual control: direction={direction}, magnitude={magnitude:.2f}%")

def manual_control_clockwise():
    """Send clockwise control command"""
    magnitude = manual_control_input.value()
    send_manual_control(1, magnitude)

def manual_control_counterclockwise():
    """Send counterclockwise control command"""
    magnitude = manual_control_input.value()
    send_manual_control(-1, magnitude)

def stop_manual_control():
    """Stop manual control (send 0 magnitude)"""
    send_manual_control(0, 0)
    manual_control_input.setValue(0)

def start_controller():
    """Send command to Arduino to start controller"""
    global running_state
    command = "START\n"
    ser.write(command.encode())
    print("Controller started")
    running_state = True
    events_data.append((time.time() - t0, "START"))  # Log start event
    start_button.setText("Stop Controller")
    start_button.clicked.disconnect()
    start_button.clicked.connect(stop_controller)
    update_status_display()

def stop_controller():
    """Send command to Arduino to stop controller"""
    global running_state
    command = "STOP\n"
    ser.write(command.encode())
    print("Controller stopped")
    running_state = False
    events_data.append((time.time() - t0, "STOP"))  # Log stop event
    start_button.setText("Start Controller")
    start_button.clicked.disconnect()
    start_button.clicked.connect(start_controller)
    update_status_display()

def pre_tension():
    """Send command to Arduino to perform pre-tensioning"""
    global running_state
    pretension_reference = pretension_input.value()*L0_stored/100  # Convert % strain to m displacement
    command = f"PRETENSION:{pretension_reference:.9f}\n"
    ser.write(command.encode())
    print(f"Pre-tensioning with reference: {pretension_reference:.9f} m")
    events_data.append((time.time() - t0, "PRETENSION"))  # Log pre-tension event
    start_button.setText("Stop Controller")
    start_button.clicked.disconnect()
    start_button.clicked.connect(stop_controller)
    running_state = True
    update_status_display()

def save_data_to_file():
    """Save measurement data to CSV file"""
    file_dialog = QtWidgets.QFileDialog()
    file_path, _ = file_dialog.getSaveFileName(
        None, 
        "Save Data", 
        "", 
        "CSV Files (*.csv);;Text Files (*.txt)"
    )
    
    if file_path:
        try:
            with open(file_path, 'w') as f:
                # Write sample parameters
                f.write(f"# Sample Parameters:\n")
                f.write(f"# Cross-sectional Area (A_cs): {A_cs_stored:.6e} m^2\n")
                f.write(f"# Length (L0): {L0_stored:.6e} m\n")
                f.write(f"# Young's Modulus used for controller: {E_stored:.6e} Pa\n")
                if E_measured:
                    f.write(f"# Young's Modulus measured (motor torque): {E_measured:.6e} Pa (R^2 = {E_fit_r2:.4f})\n")
                
                # Write events
                f.write(f"\n# Measurement Events:\n")
                for event_time, event_type in events_data:
                    f.write(f"# {event_type} at t={event_time:.6f}s\n")
                
                # Write header
                f.write("\nTime(s),Displacement(m),Omega(rad/s),Current(A),Voltage(V),P_in(W),P_mech(W),Torque(Nm),Force(N),Stress(Pa),Strain(%)\n")

                # Write data
                w = np.array(w_data)
                i_m = np.array(i_data)
                v = np.array(v_data)
                d = np.array(d_data)
                P_in, P_mech, tau_m = motor_quantities(w, i_m, v)
                force = sample_force(w, tau_m)
                stress = force / A_cs_stored if A_cs_stored > 0 else np.full_like(force, np.nan)
                strain = d / L0_stored * 100 if L0_stored > 0 else np.full_like(force, np.nan)
                for k in range(len(t_data)):
                    f.write(f"{t_data[k]:.6f},{d[k]:.9f},{w[k]:.6f},{i_m[k]:.6e},{v[k]:.6f},"
                            f"{P_in[k]:.6e},{P_mech[k]:.6e},{tau_m[k]:.6e},{force[k]:.6e},{stress[k]:.6e},{strain[k]:.6f}\n")
            
            print(f"Data saved to: {file_path}")
            QtWidgets.QMessageBox.information(None, "Success", f"Data saved to:\n{file_path}")
        except Exception as e:
            print(f"Error saving file: {e}")
            QtWidgets.QMessageBox.critical(None, "Error", f"Failed to save data:\n{str(e)}")

def reset_data():
    """Clear all measurement data and reset plots"""
    t_data.clear()
    d_data.clear()
    global E_measured
    w_data.clear()
    i_data.clear()
    v_data.clear()
    events_data.clear()
    E_measured = None
    for curve in (curve_d, curve_w, curve_i, curve_p_in, curve_p_mech, curve_s, curve_s_fit, curve_f):
        curve.setData([], [])
    E_label.setText("Measured E: -")
    command = "RESET\n"
    ser.write(command.encode())

def update_status_display():
    """Update the status label to show if controller is active or off"""
    if running_state:
        status_label.setText("Status: ACTIVE")
        status_label.setStyleSheet("QLabel { color: green; font-weight: bold; }")
    else:
        status_label.setText("Status: OFF")
        status_label.setStyleSheet("QLabel { color: red; font-weight: bold; }")    

# -----------------------
# App + Window
# -----------------------
app = QtWidgets.QApplication(sys.argv)
win = pg.GraphicsLayoutWidget(title="Motor Control Live")
win.resize(1200, 800)
win.show()


# -----------------------
# User input dialog
# -----------------------
# Create input widget
input_widget = QtWidgets.QWidget()
input_layout = QtWidgets.QVBoxLayout()

# Status display label
status_label = QtWidgets.QLabel("Status: OFF")
status_label.setStyleSheet("QLabel { color: red; font-weight: bold; }")
status_font = status_label.font()
status_font.setPointSize(12)
status_label.setFont(status_font)
input_layout.addWidget(status_label)
input_layout.addSpacing(10)  # Add some space below status

# Reference value input
ref_label = QtWidgets.QLabel("Reference Value (Desired Strain) in %:")
ref_input = QtWidgets.QDoubleSpinBox()
ref_input.setRange(0, 100)
ref_input.setValue(0)
ref_input.setSingleStep(0.1)
ref_input.setDecimals(6)

ref2_label = QtWidgets.QLabel("Reference Value (Desired Stress) in MPa:")
ref2_input = QtWidgets.QDoubleSpinBox()
ref2_input.setRange(0, 1000)
ref2_input.setValue(0)
ref2_input.setSingleStep(1)
ref2_input.setDecimals(6)

# Mode selector
mode_label = QtWidgets.QLabel("Control Mode:")
mode_combo = QtWidgets.QComboBox()
mode_combo.addItems(["Strain", "Stress"])

# Cross-sectional area input
area_label = QtWidgets.QLabel("Cross-sectional Area (A_cs) in mm^2:")
area_input = QtWidgets.QDoubleSpinBox()
area_input.setRange(0.001, 10)
area_input.setValue(0.01)
area_input.setSingleStep(0.01)
area_input.setDecimals(5)

# Length input
length_label = QtWidgets.QLabel("Length (L0) in mm:")
length_input = QtWidgets.QDoubleSpinBox()
length_input.setRange(0.001, 100)
length_input.setValue(10)
length_input.setSingleStep(0.01)
length_input.setDecimals(6)

# Initial Young's modulus guess, only used for the controller until E is measured
modulus_label = QtWidgets.QLabel("Initial guess Young's Modulus (E) in GPa:")
modulus_input = QtWidgets.QDoubleSpinBox()
modulus_input.setRange(0.001, 20)
modulus_input.setValue(1)
modulus_input.setSingleStep(0.1)
modulus_input.setDecimals(6)

# Young's modulus estimation from motor torque
use_measured_checkbox = QtWidgets.QCheckBox("Use measured E for controller gains")
use_measured_checkbox.setChecked(True)

fit_limit_label = QtWidgets.QLabel("Fit E up to strain (%) (0 = all data):")
fit_limit_input = QtWidgets.QDoubleSpinBox()
fit_limit_input.setRange(0, 100)
fit_limit_input.setValue(0)
fit_limit_input.setSingleStep(0.01)
fit_limit_input.setDecimals(4)

E_label = QtWidgets.QLabel("Measured E: -")
E_label.setStyleSheet("QLabel { font-weight: bold; }")

# Send button
send_button = QtWidgets.QPushButton("Send Parameters to Arduino")
send_button.clicked.connect(lambda: send_parameters(
    mode_combo.currentText(),  # Control mode
    ref_input.value()*1e-2,  # Convert % to decimal
    ref2_input.value()*1e6,  # Convert MPa to Pa
    area_input.value()*1e-6,  # Convert mm^2 to m^2
    length_input.value()*1e-3,  # Convert mm to m
    modulus_input.value()*1e9 # Convert GPa to Pa
))

# Start/Stop controller button
start_button = QtWidgets.QPushButton("Start Controller")
start_button.clicked.connect(start_controller)

# Update reference button
update_ref_button = QtWidgets.QPushButton("Update Reference")
update_ref_button.clicked.connect(lambda: update_reference(
    mode_combo.currentText(),  # Control mode
    ref_input.value()*1e-2,  # Convert % to decimal
    ref2_input.value()*1e6,  # Convert MPa to Pa
))

# Reset data button
reset_data_button = QtWidgets.QPushButton("Reset Measurement (Clear Data)")
reset_data_button.clicked.connect(reset_data)


# Manual control slider
manual_control_label = QtWidgets.QLabel("Manual Control:")
manual_control_input = QtWidgets.QDoubleSpinBox()
manual_control_input.setRange(0, 100)
manual_control_input.setValue(0)

# Manual control buttons
manual_control_cw_button = QtWidgets.QPushButton("Outwards")
manual_control_ccw_button = QtWidgets.QPushButton("Inwards")
manual_control_stop_button = QtWidgets.QPushButton("Stop Manual Control")
manual_control_stop_button.clicked.connect(stop_manual_control)
manual_control_cw_button.clicked.connect(manual_control_clockwise)
manual_control_ccw_button.clicked.connect(manual_control_counterclockwise)

# Pre-tension strain input
pretension_label = QtWidgets.QLabel("Pre-tension strain (%):")
pretension_input = QtWidgets.QDoubleSpinBox()
pretension_input.setRange(0, 10)
pretension_input.setValue(0.01)
pretension_input.setSingleStep(0.001)
pretension_input.setDecimals(5)

# Pre-tension button
pretension_button = QtWidgets.QPushButton("Pre-tension")
pretension_button.clicked.connect(pre_tension)

# Save data button
save_button = QtWidgets.QPushButton("Save Data to File")
save_button.clicked.connect(save_data_to_file)

# Add to layout
input_layout.addWidget(ref_label)
input_layout.addWidget(ref_input)
input_layout.addWidget(ref2_label)
input_layout.addWidget(ref2_input)
input_layout.addWidget(mode_label)
input_layout.addWidget(mode_combo)
input_layout.addWidget(area_label)
input_layout.addWidget(area_input)
input_layout.addWidget(length_label)
input_layout.addWidget(length_input)
input_layout.addWidget(modulus_label)
input_layout.addWidget(modulus_input)
input_layout.addWidget(use_measured_checkbox)
input_layout.addWidget(fit_limit_label)
input_layout.addWidget(fit_limit_input)
input_layout.addWidget(E_label)
input_layout.addWidget(send_button)
input_layout.addWidget(start_button)
input_layout.addWidget(update_ref_button)
input_layout.addWidget(manual_control_label)
input_layout.addWidget(manual_control_input)
input_layout.addWidget(manual_control_cw_button)
input_layout.addWidget(manual_control_ccw_button)
input_layout.addWidget(manual_control_stop_button)
input_layout.addWidget(pretension_label)
input_layout.addWidget(pretension_input)
input_layout.addWidget(pretension_button)
input_layout.addWidget(reset_data_button)
input_layout.addWidget(save_button)
input_layout.addStretch()

input_widget.setLayout(input_layout)

# Create main container window
main_window = QtWidgets.QMainWindow()
main_window.setWindowTitle("Motor Control Live")
main_window.resize(1200, 800)

# Create central widget with layout
central_widget = QtWidgets.QWidget()
main_layout = QtWidgets.QHBoxLayout()  # Side-by-side layout
main_layout.addWidget(input_widget, 1)  # Left: inputs
main_layout.addWidget(win, 3)           # Right: plots (wider)
central_widget.setLayout(main_layout)

main_window.setCentralWidget(central_widget)
main_window.show()


# -----------------------
# Plots
# -----------------------
plot_d = win.addPlot(title="Displacement")
curve_d = plot_d.plot(pen='y')
plot_d.setLabel('left', 'Displacement (m)')

plot_w = win.addPlot(title="Omega")
curve_w = plot_w.plot(pen='g')
plot_w.setLabel('left', 'Omega (rad/s)')

win.nextRow()

plot_i = win.addPlot(title="Motor Current (sensor)")
curve_i = plot_i.plot(pen='c')
plot_i.setLabel('left', 'Current (A)')

plot_p = win.addPlot(title="Motor Power")
plot_p.addLegend()
curve_p_in = plot_p.plot(pen='m', name="P_in = V*I")
curve_p_mech = plot_p.plot(pen='w', name="P_mech = P_in - I^2 R")
plot_p.setLabel('left', 'Power (W)')

win.nextRow()

plot_s = win.addPlot(title="Stress-Strain")
plot_s.addLegend()
curve_s = plot_s.plot(pen=None, symbol='o', symbolSize=3, symbolPen=None, symbolBrush='r', name="Measured")
curve_s_fit = plot_s.plot(pen=pg.mkPen('w', width=2, style=QtCore.Qt.DashLine), name="Linear fit (E)")
plot_s.setLabel('left', 'Stress (Pa)')
plot_s.setLabel('bottom', 'Strain (%)')

plot_f = win.addPlot(title="Force-Displacement")
curve_f = plot_f.plot(pen=None, symbol='o', symbolSize=3, symbolPen=None, symbolBrush='y')
plot_f.setLabel('left', 'Force (N)')
plot_f.setLabel('bottom', 'Displacement (m)')

# -----------------------
# Data buffers
# -----------------------
max_points = 10000

t_data = []
d_data = []
w_data = []
i_data = []   # measured motor current (A)
v_data = []   # applied motor voltage (V), NaN if not reported by the Arduino
events_data = []

t0 = time.time()

# -----------------------
# Timer update function
# -----------------------
def update():
    global t_data, d_data, w_data, i_data, v_data

    while ser.in_waiting:
        line = ser.readline().decode().strip()

        if line:
            try:
                # Arduino sends "omega,I,d" or "omega,I,d,V"
                values = list(map(float, line.split(",")))
                if len(values) == 3:
                    omega, current, d = values
                    voltage = np.nan
                elif len(values) == 4:
                    omega, current, d, voltage = values
                else:
                    raise ValueError

                t = time.time() - t0

                t_data.append(t)
                d_data.append(d)
                w_data.append(omega)
                i_data.append(current)
                v_data.append(voltage)

                # limit buffer size
                t_data = t_data[-max_points:]
                d_data = d_data[-max_points:]
                w_data = w_data[-max_points:]
                i_data = i_data[-max_points:]
                v_data = v_data[-max_points:]


            except ValueError:
                print(f"[Arduino Info]: {line}")
                # Update controller status with pretension
                if "Pre-tensioning complete, controller stopped" in line:
                    global running_state
                    running_state = False
                    update_status_display()
                    start_button.setText("Start Controller")
                    start_button.clicked.disconnect()
                    start_button.clicked.connect(start_controller)
            except Exception as e:
                # Catches any other unexpected processing errors
                print(f"[Python Error processing line]: {line} | Error: {e}")

    # update plots
    curve_d.setData(t_data, d_data)
    curve_w.setData(t_data, w_data)
    curve_i.setData(t_data, i_data)
    if not t_data:
        return

    d = np.array(d_data)
    w = np.array(w_data)
    P_in, P_mech, tau_m = motor_quantities(w, np.array(i_data), np.array(v_data))
    force = sample_force(w, tau_m)
    curve_p_in.setData(t_data, P_in)
    curve_p_mech.setData(t_data, P_mech)
    curve_f.setData(d, force)

    if L0_stored > 0 and A_cs_stored > 0:
        strain = d / L0_stored         # fraction
        stress = force / A_cs_stored   # Pa
        curve_s.setData(strain * 100, stress)
        estimate_youngs_modulus(strain, stress)
        if E_measured:
            eps_line = np.array([strain.min(), strain.max()])
            curve_s_fit.setData(eps_line * 100, E_measured * eps_line + E_fit_offset)
            E_label.setText(f"Measured E: {E_measured/1e9:.4f} GPa (R² = {E_fit_r2:.3f})")

# -----------------------
# Timer (refresh rate)
# -----------------------
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(20)  # 50 Hz GUI refresh


# -----------------------
# Start app
# -----------------------
sys.exit(app.exec_())