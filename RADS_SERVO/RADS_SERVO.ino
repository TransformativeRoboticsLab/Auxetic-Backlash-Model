#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>
const int X = 3; // X axis
const int Y = 3; // Y axis
const int SERVONUM = 2; // Z axis
int state[SERVONUM][X][Y]

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();

#define SERVOMIN  150 // Minimum pulse length count (out of 4096)
#define SERVOMAX  600 // Maximum pulse length count (out of 4096)
#define SERVO_FREQ 50 // Analog servos run at ~50 Hz updates
#define NUM_SERVOS 16 // Total number of servo channels

/*
The convention is that servos 1-9 correspond to the linear stage in this way:

1  2  3
4  5  6
7  8  9

within these, the number refers to the linear distance from the bottom (in mm)

servos 10-16 correspond to the expanding stage in this way

10 11 12
X  13 X 
14 15 16

within these, the number refers to the angle between the two pieces

Note that two edges don't have connected servos; this is because we don't have enough channels.

*/

void setup() {
  Serial.begin(9600);
  Serial.println("16 channel simultaneous servo test!");

  pwm.begin();
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);

  delay(10);
}

void loop() {
  //update current position
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    parseState(line);
  }
  updateServoPositions();
  // Move all 16 servos from 0° to 180° simultaneously
  Serial.println("Moving all servos to 180 degrees");
  for (uint16_t pulselen = SERVOMIN; pulselen <= SERVOMAX; pulselen++) {
    for (uint8_t servonum = 0; servonum < NUM_SERVOS; servonum++) {
      pwm.setPWM(servonum, 0, pulselen);
    }
    delay(5); // Small delay for smooth motion
  }

  delay(1000); // Pause at 180°

  // Move all 16 servos from 180° back to 0° simultaneously
  Serial.println("Moving all servos to 0 degrees");
  for (uint16_t pulselen = SERVOMAX; pulselen >= SERVOMIN; pulselen--) {
    for (uint8_t servonum = 0; servonum < NUM_SERVOS; servonum++) {
      pwm.setPWM(servonum, 0, pulselen);
    }
    delay(5); // Small delay for smooth motion
  }

  delay(1000); // Pause at 0°
}

void updateServoPositions(){
  int servoIndex = 0
  //updating linear servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        
        pwm.setPwm(servoIndex, state[servo][x][y])
      }
    }
  //updating expanding servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        pwm.setPwm(servoIndex, state[servo][x][y])
      }
    }
}

void setServo(int servoIndex, int )

void parseState(String line) {
  line.trim();

  int idx = 0;
  int last = 0;
  int values[18];

  // Tokenize manually (lighter than strtok)
  for (int i = 0; i < line.length(); i++) {
    if (line[i] == ',' || i == line.length() - 1) {
      String token = line.substring(last, i);
      last = i + 1;

      token.trim();

      if (idx == 0 && token != "S") return;
      if (token == "E") break;

      if (idx > 0 && idx <= 18) {
        values[idx - 1] = token.toInt();
      }

      idx++;
    }
  }

  if (idx < 20) return;  // malformed packet

  // Rebuild 3D array
  int k = 0;
  for (int z = 0; z < SERVONUM; z++) {
    for (int x = 0; x < X; X++) {
      for (int y = 0; y < Y; Y++) {
        state[z][x][y] = values[k++];
      }
    }
  }

  Serial.println("State updated.");
}