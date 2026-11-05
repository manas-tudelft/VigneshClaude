// StrainingStage.ino - reference firmware for the straining stage with load-cell force feedback
// Serial protocol v2, used by ControlInterface4.py. Its logic mirrors SimulatedArduino in that file.
//
// Hardware (see LOADCELL_WIRING.md):
//   Arduino Uno R3
//   Encoder A/B via optocoupler OUT1/OUT2 -> D2/D3 (interrupts, x4 decoding)
//   Motor driver IN3/IN4 -> D9/D10 (Timer1 PWM, ENB jumpered)
//   FUTEK LSB200 load cell -> HX711 (gain 128, RATE tied to VCC = 80 SPS): DOUT -> D4, SCK -> D5
// Library: "HX711 Arduino Library" by Bogdan Necula (Arduino Library Manager)
//
// Not compiled in the session that wrote it: check the TODO constants against your old sketch.
//
// PC -> Arduino (one command per line):
//   r,k_omega,k_I,k_d,k_F,k_int,mode   gains + reference (mode 0: track d [m], mode 1: track F [N])
//   R:<r>                              new reference
//   RAMP:<target>,<rate>               ramp the reference (m, m/s or N, N/s), replies "RAMP done"
//   JOG:<dir>,<distance m>,<speed m/s> distance-based manual move, replies "JOG done" (distance 0 = continuous)
//   PRELOAD:<N>                        move outwards until F >= N, replies "Pre-tensioning complete, controller stopped"
//   START, STOP, RESET, TARE, CAL:<known N>, FLIM:<N>, MANUAL:<dir>,<magnitude %> (legacy open loop)
// Arduino -> PC, one line per HX711 sample (80 Hz):  omega,I,d,F,u   (rad/s, A estimated, m, N, V)

#include <HX711.h>
#include <EEPROM.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ---------- pins ----------
const uint8_t PIN_ENC_A = 2;       // optocoupler OUT1
const uint8_t PIN_ENC_B = 3;       // optocoupler OUT2
const uint8_t PIN_IN3 = 9;         // motor driver IN3 (OC1A)
const uint8_t PIN_IN4 = 10;        // motor driver IN4 (OC1B)
const uint8_t PIN_HX_DOUT = 4;     // HX711 DT
const uint8_t PIN_HX_SCK = 5;      // HX711 SCK

// ---------- stage constants ----------
const float COUNTS_PER_REV = 2000.0;   // TODO: encoder counts per OUTPUT-shaft revolution after x4 decoding
const int8_t ENC_SIGN = 1;             // TODO: +1 or -1 so that moving outwards gives positive d
const int8_t MOTOR_SIGN = 1;           // TODO: +1 or -1 so that positive u moves outwards
const float P_LEAD = 0.0007;           // m, lead-screw parameter p (same as the GUI)
const float G_SCREW = P_LEAD / PI;     // m of elongation per rad of the output shaft
const float RAD_PER_COUNT = 2.0 * PI / COUNTS_PER_REV;
const float D_PER_COUNT = G_SCREW * RAD_PER_COUNT;
const float V_DC = 24.0;
const float K_T = 243.8 * 53.601e-3;   // Nm/A at the output shaft
const float R_MOTOR = 55.041;          // ohm
const float TS = 0.001;                // s, control period (1 kHz)
const unsigned long TS_US = 1000;
const float OMEGA_TAU = 0.005;         // s, speed low-pass (TODO: or use the omega estimate of your sketch)
const float D_LIMIT = 0.02;            // m, software travel limit
const float PRELOAD_SPEED = 1e-5;      // m/s
const uint16_t PWM_TOP = 1023;         // Timer1 10-bit phase-correct PWM, about 7.8 kHz

// ---------- load cell: FUTEK LSB200 (values from the calibration certificate) ----------
const float LC_CAPACITY_N = 4.903;          // TODO: rated capacity of your LSB200 in N
const float LC_RATED_OUTPUT_MV_V = 2.0;     // TODO: rated output in mV/V from the certificate
// HX711 at gain 128: full scale 2^23 counts = 0.5*E/128 (ratiometric), so counts/N = RO*1e-3*2^31/F_rated
const float SCALE_FROM_CERT = LC_RATED_OUTPUT_MV_V * 1e-3 * 2147483648.0 / LC_CAPACITY_N;
float force_limit = 5.0;                    // N, FLIM (maximum load of the LSB200 set-up)

// Stiffness-free displacement gains K0 = [k_omega, k_I, k_d, k_F, k_int], designed by
// ControlInterface4.py without any sample knowledge. Jog and preload work with these
// before the GUI has sent gains.
const float K0[5] = {0.038540367, 0.15629437, 751289.64, -0.00090061443, -547722.56};

// ---------- state ----------
enum State : uint8_t { IDLE, RUN, JOG, PRELOAD, MANUAL };
enum RampKind : uint8_t { RK_NONE, RK_RAMP, RK_JOG_RUN, RK_JOG_IDLE };

HX711 hx;
State state = IDLE;
float K[5], K_disp[5];
uint8_t mode = 0;
float r = 0.0, z = 0.0, u_out = 0.0, manual_u = 0.0;
RampKind ramp_kind = RK_NONE;
float ramp_r0 = 0.0, ramp_target = 0.0, ramp_rate = 0.0;
unsigned long ramp_steps = 0;
float preload_force = 0.0, preload_r0 = 0.0;
unsigned long preload_steps = 0;
float jog_settle = 0.0, jog_timer = 0.0;

volatile long enc_count = 0;
volatile uint8_t enc_state = 0;
long enc_zero = 0;
long count_hist[10];
uint8_t hist_i = 0;
float omega_f = 0.0, d_meas = 0.0, F_meas = 0.0;

long hx_offset = 0;
float hx_scale = SCALE_FROM_CERT;          // counts per N
uint8_t collect_kind = 0;                  // 1 = TARE, 2 = CAL
float collect_value = 0.0;
long collect_sum = 0;
uint8_t collect_n = 0;

char buf[112];
uint8_t blen = 0;
unsigned long t_next = 0;

struct Calib { uint16_t magic; long offset; float scale; };
const uint16_t EE_MAGIC = 0x4C43;

// ---------- encoder (x4 quadrature decoding) ----------
const int8_t ENC_TABLE[16] = {0, -1, 1, 0, 1, 0, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0};

void encoderISR() {
  uint8_t pins = PIND;                                   // D2 = PD2 (A), D3 = PD3 (B)
  uint8_t s = (((pins >> 2) & 1) << 1) | ((pins >> 3) & 1);
  enc_state = ((enc_state << 2) | s) & 0x0F;
  enc_count += ENC_TABLE[enc_state];
}

// ---------- motor ----------
void pwm_init() {
  pinMode(PIN_IN3, OUTPUT);
  pinMode(PIN_IN4, OUTPUT);
  TCCR1A = _BV(COM1A1) | _BV(COM1B1) | _BV(WGM11);       // phase-correct PWM, TOP = ICR1 (mode 10)
  TCCR1B = _BV(WGM13) | _BV(CS10);                       // no prescaler: 16 MHz / (2 * 1023) = 7.8 kHz
  ICR1 = PWM_TOP;
  OCR1A = 0;
  OCR1B = 0;
}

void motor(float u) {
  u = constrain(u, -V_DC, V_DC) * MOTOR_SIGN;
  uint16_t duty = (uint16_t)(fabs(u) / V_DC * PWM_TOP + 0.5);
  if (u >= 0) { OCR1B = 0; OCR1A = duty; }               // IN3 PWM, IN4 low: outwards
  else        { OCR1A = 0; OCR1B = duty; }               // IN4 PWM, IN3 low: inwards
}

// ---------- calibration storage ----------
void load_calib() {
  Calib c;
  EEPROM.get(0, c);
  if (c.magic == EE_MAGIC && isfinite(c.scale) && c.scale != 0.0) { hx_offset = c.offset; hx_scale = c.scale; }
}

void save_calib() {
  Calib c = {EE_MAGIC, hx_offset, hx_scale};
  EEPROM.put(0, c);
}

// ---------- control law (identical to SimulatedArduino._law) ----------
float law(const float *Kc, uint8_t m, float rr, float zz, float w, float I, float *y) {
  float e_d, e_F;
  if (m == 0) { e_d = d_meas - rr; e_F = F_meas; *y = d_meas; }
  else        { e_d = d_meas; e_F = F_meas - rr; *y = F_meas; }
  return -(Kc[0] * w + Kc[1] * I + Kc[2] * e_d + Kc[3] * e_F + Kc[4] * zz);
}

float current_estimate() {                               // quasi-static: L*dI/dt ~ 0 (L/R ~ 2 us)
  return (u_out - K_T * omega_f) / R_MOTOR;
}

void bumpless(const float *Kn, uint8_t m, float rr, float u_now) {
  // choose the integrator so that the new law gives the present output u_now
  float w = omega_f, I = current_estimate();
  float e_d = (m == 0) ? d_meas - rr : d_meas;
  float e_F = (m == 0) ? F_meas : F_meas - rr;
  float rest = Kn[0] * w + Kn[1] * I + Kn[2] * e_d + Kn[3] * e_F;
  z = (Kn[4] != 0.0) ? -(u_now + rest) / Kn[4] : 0.0;
}

void stop_all() {
  state = IDLE;
  ramp_kind = RK_NONE;
  manual_u = 0.0;
  u_out = 0.0;
  motor(0.0);
}

void start_ramp(float target, float rate, RampKind kind) {
  ramp_r0 = r;
  ramp_target = target;
  ramp_rate = fabs(rate);
  ramp_steps = 0;
  ramp_kind = kind;
}

float control_output() {
  if (state == IDLE) return 0.0;
  if (state == MANUAL) return manual_u;
  if (state == PRELOAD) {
    preload_steps++;
    r = preload_r0 + PRELOAD_SPEED * preload_steps * TS;   // step counter: no float accumulation
  } else if (ramp_kind != RK_NONE) {
    ramp_steps++;
    float travel = ramp_rate * ramp_steps * TS;
    float span = ramp_target - ramp_r0;
    if (travel >= fabs(span)) {
      r = ramp_target;
      RampKind k = ramp_kind;
      ramp_kind = RK_NONE;
      if (k == RK_RAMP) Serial.println(F("RAMP done"));
      else if (k == RK_JOG_RUN) Serial.println(F("JOG done"));
      else { jog_settle = 0.0; jog_timer = 0.0; }
    } else {
      r = ramp_r0 + (span > 0 ? travel : -travel);
    }
  }
  if (state == JOG && ramp_kind == RK_NONE) {               // settle at the jog target, then motor off
    jog_timer += TS;
    jog_settle = (fabs(d_meas - r) <= 2.0 * D_PER_COUNT) ? jog_settle + TS : 0.0;
    if (jog_settle >= 0.2 || jog_timer >= 5.0) {
      stop_all();
      Serial.println(F("JOG done"));
      return 0.0;
    }
  }
  bool disp_only = (state == JOG || state == PRELOAD);
  const float *Kc = disp_only ? K_disp : K;
  uint8_t m = disp_only ? 0 : mode;
  float w = omega_f, I = current_estimate(), y;
  float u = law(Kc, m, r, z, w, I, &y);
  float z_new = z + TS * (r - y);
  float u_new = law(Kc, m, r, z_new, w, I, &y);
  if (fabs(u_new) <= V_DC || fabs(u_new) < fabs(u)) { z = z_new; u = u_new; }   // anti-windup
  return constrain(u, -V_DC, V_DC);
}

void control_step() {
  noInterrupts();
  long c = enc_count;
  interrupts();
  c = ENC_SIGN * (c - enc_zero);
  d_meas = c * D_PER_COUNT;
  long c_old = count_hist[hist_i];                          // speed over a 10 ms window, then low-pass
  count_hist[hist_i] = c;
  hist_i = (hist_i + 1) % 10;
  float w_raw = (c - c_old) * RAD_PER_COUNT / (10.0 * TS);
  omega_f += (w_raw - omega_f) * (TS / OMEGA_TAU);
  float u = control_output();
  u_out = (state == IDLE) ? 0.0 : u;
  motor(u_out);
}

// ---------- load cell ----------
void send_telemetry() {
  Serial.print(omega_f, 6); Serial.print(',');
  Serial.print(current_estimate(), 6); Serial.print(',');
  Serial.print(d_meas, 9); Serial.print(',');
  Serial.print(F_meas, 6); Serial.print(',');
  Serial.println(u_out, 4);
}

void on_force_sample(long raw) {
  F_meas = (raw - hx_offset) / hx_scale;
  if (collect_kind) {
    collect_sum += raw;
    if (++collect_n >= 16) {
      long mean = collect_sum / 16;
      if (collect_kind == 1) {
        hx_offset = mean;
        save_calib();
        Serial.print(F("TARE done, offset = ")); Serial.println(mean);
      } else if (labs(mean - hx_offset) > 10 && collect_value != 0.0) {
        hx_scale = (mean - hx_offset) / collect_value;
        save_calib();
        Serial.print(F("CAL done, scale = ")); Serial.print(hx_scale, 3); Serial.println(F(" counts/N"));
      } else {
        Serial.println(F("ERR CAL: no load on the load cell"));
      }
      collect_kind = 0;
    }
  }
  if (state != IDLE && fabs(F_meas) > force_limit) {
    stop_all();
    Serial.print(F("FORCE LIMIT: F = ")); Serial.print(F_meas, 4);
    Serial.print(F(" N > ")); Serial.print(force_limit, 4); Serial.println(F(" N, motor stopped"));
  }
  if (state == PRELOAD && F_meas >= preload_force) {
    stop_all();
    Serial.println(F("Pre-tensioning complete, controller stopped"));
  }
  if (state != IDLE && fabs(d_meas) > D_LIMIT) {
    stop_all();
    Serial.println(F("DISPLACEMENT LIMIT, motor stopped"));
  }
  send_telemetry();
}

// ---------- commands ----------
bool parse_floats(const char *s, float *out, uint8_t n) {
  char *end;
  for (uint8_t i = 0; i < n; i++) {
    out[i] = strtod(s, &end);
    if (end == s) return false;
    s = end;
    if (i < n - 1) {
      if (*s != ',') return false;
      s++;
    }
  }
  return *s == 0;
}

void cmd_gains(char *cmd) {
  float v[7];
  if (!parse_floats(cmd, v, 7)) { Serial.print(F("ERR bad command: ")); Serial.println(cmd); return; }
  uint8_t m = (uint8_t)lround(v[6]);
  if (m > 1) { Serial.println(F("ERR bad mode")); return; }
  float Kn[5] = {v[1], v[2], v[3], v[4], v[5]};
  float rr = v[0];
  if (m != mode) ramp_kind = RK_NONE;
  else if (ramp_kind != RK_NONE) rr = r;                    // gain schedule during a ramp: keep the ramp
  if (state == RUN) bumpless(Kn, m, rr, u_out);
  memcpy(K, Kn, sizeof K);
  mode = m;
  r = rr;
  if (ramp_kind != RK_NONE) { ramp_r0 = r; ramp_steps = 0; }
  if (m == 0) memcpy(K_disp, Kn, sizeof K_disp);
  Serial.print(F("GAINS OK mode=")); Serial.println(m);
}

void cmd_jog(float dir, float dist, float speed) {
  if (speed <= 0) { Serial.println(F("ERR JOG: speed must be > 0")); return; }
  float sgn = dir > 0 ? 1.0 : -1.0;
  if (state == RUN) {
    if (mode != 0) { Serial.println(F("ERR JOG: stop force control first")); return; }
    start_ramp(dist > 0 ? r + sgn * dist : sgn * D_LIMIT, speed, RK_JOG_RUN);
    return;
  }
  stop_all();
  state = JOG;
  r = d_meas;
  bumpless(K_disp, 0, r, 0.0);
  start_ramp(dist > 0 ? r + sgn * dist : sgn * D_LIMIT, speed, RK_JOG_IDLE);
  jog_settle = jog_timer = 0.0;
}

void cmd_preload(float F_target) {
  if (F_meas >= F_target) { Serial.println(F("Pre-tensioning complete, controller stopped")); return; }
  stop_all();
  state = PRELOAD;
  preload_force = F_target;
  r = preload_r0 = d_meas;
  preload_steps = 0;
  bumpless(K_disp, 0, r, 0.0);
}

void handle(char *cmd) {
  float v[3];
  if (!strcmp(cmd, "START")) {
    if (state != RUN) {
      stop_all();
      state = RUN;
      bumpless(K, mode, r, 0.0);
      Serial.println(F("Controller started"));
    }
  } else if (!strcmp(cmd, "STOP")) {
    stop_all();
    Serial.println(F("Controller stopped"));
  } else if (!strcmp(cmd, "RESET")) {
    stop_all();
    noInterrupts();
    enc_zero = enc_count;
    interrupts();
    for (uint8_t i = 0; i < 10; i++) count_hist[i] = 0;
    z = r = 0.0;
    Serial.println(F("RESET done"));
  } else if (!strcmp(cmd, "TARE")) {
    collect_kind = 1; collect_sum = 0; collect_n = 0;
  } else if (!strncmp(cmd, "CAL:", 4)) {
    collect_value = atof(cmd + 4); collect_kind = 2; collect_sum = 0; collect_n = 0;
  } else if (!strncmp(cmd, "FLIM:", 5)) {
    force_limit = fabs(atof(cmd + 5));
    Serial.print(F("FLIM = ")); Serial.print(force_limit, 4); Serial.println(F(" N"));
  } else if (!strncmp(cmd, "RAMP:", 5)) {
    if (!parse_floats(cmd + 5, v, 2)) Serial.println(F("ERR bad RAMP"));
    else if (state != RUN) Serial.println(F("ERR RAMP: controller not running"));
    else start_ramp(v[0], v[1], RK_RAMP);
  } else if (!strncmp(cmd, "R:", 2)) {
    r = atof(cmd + 2);
    ramp_kind = RK_NONE;
  } else if (!strncmp(cmd, "JOG:", 4)) {
    if (parse_floats(cmd + 4, v, 3)) cmd_jog(v[0], fabs(v[1]), fabs(v[2]));
    else Serial.println(F("ERR bad JOG"));
  } else if (!strncmp(cmd, "PRELOAD:", 8)) {
    cmd_preload(atof(cmd + 8));
  } else if (!strncmp(cmd, "MANUAL:", 7)) {            // legacy open-loop duty, not used by the GUI
    if (!parse_floats(cmd + 7, v, 2)) { Serial.println(F("ERR bad MANUAL")); return; }
    stop_all();
    if (v[0] != 0 && v[1] > 0) {
      state = MANUAL;
      manual_u = (v[0] > 0 ? 1.0 : -1.0) * min(v[1], 100.0) / 100.0 * V_DC;
    }
  } else {
    cmd_gains(cmd);
  }
}

void read_serial() {
  while (Serial.available()) {
    char c = Serial.read();
    if (c == '\n' || c == '\r') {
      if (blen) { buf[blen] = 0; handle(buf); blen = 0; }
    } else if (blen < sizeof(buf) - 1) {
      buf[blen++] = c;
    }
  }
}

// ---------- main ----------
void setup() {
  Serial.begin(115200);
  pinMode(PIN_ENC_A, INPUT_PULLUP);                    // TODO: keep the input mode of your old sketch
  pinMode(PIN_ENC_B, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(PIN_ENC_A), encoderISR, CHANGE);
  attachInterrupt(digitalPinToInterrupt(PIN_ENC_B), encoderISR, CHANGE);
  pwm_init();
  motor(0.0);
  hx.begin(PIN_HX_DOUT, PIN_HX_SCK);                   // channel A, gain 128
  load_calib();
  memcpy(K, K0, sizeof K);
  memcpy(K_disp, K0, sizeof K_disp);
  Serial.println(F("Straining stage FW v2 (load cell, protocol v2)"));
  t_next = micros() + TS_US;
}

void loop() {
  read_serial();
  unsigned long now = micros();
  if ((long)(now - t_next) >= 0) {
    t_next += TS_US;
    if ((long)(now - t_next) > 10000L) t_next = now + TS_US;   // resynchronise after a long stall
    control_step();
  }
  if (hx.is_ready()) on_force_sample(hx.read());       // about 80 times per second
}
