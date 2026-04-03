#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();

#define SERVO_FREQ   50
#define LINSERVOMID  288
#define NUM_SERVOS   9

const int SERVO_CHANNELS[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const int ENCODER_PINS[]   = {9, 10, 2, 3, 4, 5, 6, 7, 8};

// Per-servo trim offsets — adjust until each servo holds still at rest.
// Positive = bump PWM up, negative = bump PWM down.
int servoTrim[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

int midpoint(int idx) {
  return LINSERVOMID + servoTrim[idx];
}

// --- Trim calibration state ---
int trimServo = -1; // which servo is being trimmed (-1 = none)

double readEncoder(int idx) {
  return map(pulseIn(ENCODER_PINS[idx], HIGH, 3000), 23, 1048, 0, 360);
}

void park(int idx) {
  pwm.setPWM(SERVO_CHANNELS[idx], 0, midpoint(idx));
}

void drive720(int idx, int direction) {
  int ch = SERVO_CHANNELS[idx];
  int mid = midpoint(idx);

  Serial.print(">>> Driving servo ");
  Serial.print(idx);
  Serial.print(" (PWM ch ");
  Serial.print(ch);
  Serial.print(", encoder pin ");
  Serial.print(ENCODER_PINS[idx]);
  Serial.print(") ");
  Serial.print(direction > 0 ? "CW" : "CCW");
  Serial.println(" 720°");

  const double TARGET    = 720.0 * direction;
  const double TOLERANCE = 2.0;

  double origin = readEncoder(idx);
  double prev   = origin;
  double offset = 0;

  while (true) {
    double raw   = readEncoder(idx);
    double delta = raw - prev;
    if (delta >  180) offset -= 360;
    if (delta < -180) offset += 360;
    prev = raw;

    double position = (raw + offset) - origin;
    double error    = TARGET - position;

    if (abs(error) < TOLERANCE) {
      pwm.setPWM(ch, 0, mid);
      Serial.print(">>> Servo ");
      Serial.print(idx);
      Serial.println(" done.");
      return;
    }

    double power = error * 0.08;
    if (abs(power) < 5) power = (power > 0) ? 8 : -8;
    power = constrain(power, -50, 50);

    pwm.setPWM(ch, 0, mid - (int)power);
    delay(10);
  }
}

void printTrimStatus(int idx) {
  Serial.print("Servo ");
  Serial.print(idx);
  Serial.print(" trim = ");
  Serial.print(servoTrim[idx]);
  Serial.print(" (midpoint = ");
  Serial.print(midpoint(idx));
  Serial.println(")");
}

void printHelp() {
  Serial.println("Commands:");
  Serial.println("  0-8   Drive that servo CW 720°");
  Serial.println("  r0-r8 Drive that servo CCW 720°");
  Serial.println("  a     Drive all servos CW sequentially");
  Serial.println("  b     Drive all servos CCW sequentially");
  Serial.println("  t0-t8 Enter trim mode for that servo");
  Serial.println("  +/-   Adjust trim ±1 (in trim mode)");
  Serial.println("  s     Kill signal to trim servo");
  Serial.println("  q     Exit trim mode");
  Serial.println("  p     Print all trim values");
  Serial.println("  ?     Show this help");
}

void setup() {
  Serial.begin(9600);
  pwm.begin();
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);

  for (int i = 0; i < NUM_SERVOS; i++) park(i);
  delay(500);

  printHelp();
}

void loop() {
  if (!Serial.available()) return;

  String input = Serial.readStringUntil('\n');
  input.trim();
  if (input.length() == 0) return;

  char c = input[0];

  // --- Trim mode entry: t0 through t8 ---
  if (c == 't' && input.length() >= 2) {
    int idx = input[1] - '0';
    if (idx >= 0 && idx < NUM_SERVOS) {
      trimServo = idx;
      park(trimServo);
      Serial.print("Trim mode: servo ");
      Serial.println(trimServo);
      printTrimStatus(trimServo);
      Serial.println("Use +/- to adjust, s to kill signal, q to exit.");
    }
    return;
  }

  // --- Trim mode adjustments ---
  if (trimServo >= 0) {
    if (c == '+') {
      servoTrim[trimServo]++;
      park(trimServo);
      printTrimStatus(trimServo);
    }
    else if (c == '-') {
      servoTrim[trimServo]--;
      park(trimServo);
      printTrimStatus(trimServo);
    }
    else if (c == 's') {
      pwm.setPWM(SERVO_CHANNELS[trimServo], 0, 0);
      Serial.println("Signal killed.");
    }
    else if (c == 'q') {
      park(trimServo);
      Serial.print("Exited trim mode. Final: ");
      printTrimStatus(trimServo);
      trimServo = -1;
    }
    return;
  }

  // --- CCW command: r0 through r8 ---
  if (c == 'r' && input.length() >= 2) {
    int idx = input[1] - '0';
    if (idx >= 0 && idx < NUM_SERVOS) {
      drive720(idx, -1);
    }
    return;
  }

  // --- Normal commands ---
  if (c >= '0' && c <= '8') {
    drive720(c - '0', 1);
  }
  else if (c == 'a' || c == 'A') {
    Serial.println(">>> Running all servos CW sequentially...");
    for (int i = 0; i < NUM_SERVOS; i++) {
      drive720(i, 1);
      delay(250);
    }
    Serial.println(">>> All done.");
  }
  else if (c == 'b' || c == 'B') {
    Serial.println(">>> Running all servos CCW sequentially...");
    for (int i = 0; i < NUM_SERVOS; i++) {
      drive720(i, -1);
      delay(250);
    }
    Serial.println(">>> All done.");
  }
  else if (c == 'p') {
    for (int i = 0; i < NUM_SERVOS; i++) printTrimStatus(i);
  }
  else if (c == '?') {
    printHelp();
  }
}
