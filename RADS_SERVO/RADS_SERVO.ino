#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>
const int X = 3; // X axis
const int Y = 3; // Y axis
const int SERVONUM = 2; // linear or expanding servo; 0 is linear, 1 is expanding

const int ENCODERPINS[] = {9, 10, 2, 3, 4, 5, 6, 7, 8}; // channels 0..8 used for linear motors

const int LINEAR_SERVO_LAYOUT[X][Y] = {
  {0, 1, 2},
  {3, 4, 5},
  {6, 7, 8}
};

const int EXPANDING_SERVO_LAYOUT[X][Y] = {
  {9, 10, 11},
  {-1, 12, -1},
  {13, 14, 15}
};

int targetStates[SERVONUM][X][Y];
int linCurrentStates[X][Y]; //keeps track of the current states of the linear actuator
unsigned long lastEncoderPrintMs = 0;

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();

#define SERVOMIN  205 // Minimum pulse length count (out of 4096)
#define SERVOMAX  410 // Maximum pulse length count (out of 4096)
#define LINSERVOMIN 200 //Min pulse len (conservative)
#define LINSERVOMID 288 //point of linear servo middle
#define LINSERVOMAX 376 // Max pulse len (conservative)
#define SERVO_FREQ 50 // Analog servos run at ~50 Hz updates
#define NUM_SERVOS 16 // Total number of servo channels
#define ENCODERMIN 500
#define ENCODERMAX 2500

/*
The convention is that channels 0-8 correspond to the linear stage in this way:

0  1  2
3  4  5
6  7  8

within these, the number refers to the linear distance from the bottom (in mm)

channels 9-15 correspond to the expanding stage in this way

9  10 11
X  12 X 
13 14 15

within these, the number refers to the angle between the two pieces

Note that two edges don't have connected servos; this is because we don't have enough channels.
*/

void setup() {
  Serial.begin(9600);

  pwm.begin();
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);

  // Startup state: all linear channels at neutral PWM, all expanding channels at 0 deg.
  for (int x = 0; x < X; x++) {
    for (int y = 0; y < Y; y++) {
      int linearChannel = LINEAR_SERVO_LAYOUT[x][y];
      pwm.setPWM(linearChannel, 0, LINSERVOMID);
      linCurrentStates[x][y] = 0;
    }
  }

  for (int x = 0; x < X; x++) {
    for (int y = 0; y < Y; y++) {
      int expandingChannel = EXPANDING_SERVO_LAYOUT[x][y];
      if (expandingChannel >= 0) {
        rotateMove(expandingChannel, 0);
      }
      targetStates[1][x][y] = 0;
    }
  }
}

void loop() {
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    parseState(line);
  }
  updateServoPositions();
}

void updateServoPositions(){
  // updating linear servos sequentially
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        int servoLabel = LINEAR_SERVO_LAYOUT[x][y];
        int currentState = linCurrentStates[x][y];
        int targetState = targetStates[0][x][y];
        int deltaMM = targetState - currentState;
        int deltaAngle = deltaMM * 500; //converts from input mm to angle, which is what linear_move does
        Serial.print("START_LINEAR,");
        Serial.print(servoLabel);
        Serial.print(",from=");
        Serial.print(currentState);
        Serial.print(",to=");
        Serial.println(targetState);
        linearMove(servoLabel, deltaAngle);
        linCurrentStates[x][y] = targetState;
        Serial.print("FINISH_LINEAR,");
        Serial.println(servoLabel);
      }
    }

  //updating expanding servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        int servoLabel = EXPANDING_SERVO_LAYOUT[x][y];
        if (servoLabel >= 0) {
          Serial.print("START_EXPANDING,");
          Serial.print(servoLabel);
          Serial.print(",angle=");
          Serial.println(targetStates[1][x][y]);
          rotateMove(servoLabel, targetStates[1][x][y]);
          Serial.print("FINISH_EXPANDING,");
          Serial.println(servoLabel);
        }
      }
    }
}

// delta is in degrees
void linearMove(int servoNum, int delta) {
  double encoderValueInitial = readEncoder(servoNum);
  double encoderValue = readEncoder(servoNum) - encoderValueInitial;
  double previousEncoderValue = encoderValue;
  while (abs(encoderValue - delta) > 1) { //error value, may need changing
    int pwmWidth = map(encoderValue - delta, -90, 90, LINSERVOMIN, LINSERVOMAX);
    double power = (delta - encoderValue) * 0.08; //amount of power to supply to the servo (as difference from center) e.g 20 power is center (288) + 20 = 308
    if (abs(power) < 5) { //if error is too small, motor will move and will never actually hit desired position
      power = ((power > 0) - (power < 0)) * 8; //assigns power to -10 or 10 (essencially minimum)
    }
    power = constrain(power, -50, 50); //makes sure motor doesn't spin too fast
    int constrainedPwm = LINSERVOMID - power;
    pwm.setPWM(servoNum, 0, constrainedPwm);

    previousEncoderValue = encoderValue;
    encoderValue = readEncoder(servoNum) - encoderValueInitial;

    if (encoderValue - previousEncoderValue > 90) { //if encoderValue suddenly jumps up a lot, assume it wrapped from 0 to 360; subtract 360 from both values
      encoderValueInitial += 360;
    } else if (encoderValue - previousEncoderValue < -90 ) { //else if encoderValue suddenly jumps down a lot, assume it wrapped from 360 to 0, add 360 to both values
      encoderValueInitial -= 360;
    }

    encoderValue = readEncoder(servoNum) - encoderValueInitial;
    delay(10); //just so it doesn't kill itself
  }

  pwm.setPWM(servoNum, 0, LINSERVOMID);
}

double readEncoder(int servoNum) {
  return map(pulseIn(ENCODERPINS[servoNum], HIGH, 3000), 23, 1048, 0, 360);
}

void rotateMove(int servoNum, int angle) {
  int pulseLen = map(angle, 0, 180, SERVOMIN, SERVOMAX);
  pwm.setPWM(servoNum, 0, pulseLen);
}

void parseState(String line) {
  line.trim();

  if (line.length() == 0) {
    return;
  }

  int last = 0;
  int values[18];
  int valueCount = 0;
  bool sawStart = false;
  bool sawEnd = false;

  // Tokenize manually (lighter than strtok)
  for (int i = 0; i <= line.length(); i++) {
    if (i == line.length() || line[i] == ',') {
      String token = line.substring(last, i);
      last = i + 1;

      token.trim();
      if (token.length() == 0) {
        continue;
      }

      if (!sawStart) {
        if (token != "S") {
          return;
        }
        sawStart = true;
        continue;
      }

      if (token == "E") {
        sawEnd = true;
        break;
      }

      if (valueCount >= 18) {
        return;
      }
      values[valueCount++] = token.toInt();
    }
  }

  if (!sawStart) {
    return;
  }
  if (!sawEnd) {
    return;
  }
  if (valueCount != 18) {
    return;
  }

  // Rebuild 3D array
  int k = 0;
  for (int z = 0; z < SERVONUM; z++) {
    for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        targetStates[z][x][y] = values[k++];
      }
    }
  }

  updateServoPositions();
}