#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>
#include <ServoInput.h>
const int X = 3; // X axis
const int Y = 3; // Y axis
const int SERVONUM = 2; // linear or expanding servo; 0 is linear, 1 is expanding

const int ENCODERPINS[] = {2, 2, 2, 2, 2, 2, 2, 2, 2}; //index is motor, value is pin

int targetStates[SERVONUM][X][Y];
int linCurrentStates[X][Y]; //keeps track of the current states of the linear actuator

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
  Serial.println("servo test!");

  pwm.begin();
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);

  delay(10);
  pwm.setPWM(1, 0, 288);
  delay(3000);
  linearMove(1, 1600);
}

void loop() {
  Serial.print("completed");
  //update current position
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    parseState(line);
  }
  //Serial.println(readEncoder(1));
  //linearMove(1, 120);
  //updateServoPositions();
}

void updateServoPositions(){
  int servoIndex = 0;
  //updating linear servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        int currentState = linCurrentStates[x][y];
        int targetState = targetStates[0][x][y];
        int deltaMM = targetState - currentState;
        int deltaAngle = deltaMM * 500;
        linearMove(servoIndex, deltaAngle);
        linCurrentStates[x][y] = targetState;
        servoIndex++;
      }
    }
  //updating expanding servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        rotateMove(servoIndex, targetStates[1][x][y]);
        servoIndex++;
      }
    }
}

// delta is in degrees
void linearMove(int servoNum, int delta) {
  int encoderValueInitial = readEncoder(servoNum);
  int encoderValue = readEncoder(servoNum) - encoderValueInitial;
  int previousEncoderValue = encoderValue;
  while (abs(encoderValue - delta) > 1) { //error value, may need changing
    //set motor speed based off of error; covers too extreme motor values or too small ones 
    int pwmWidth = map(encoderValue - delta, -90, 90, LINSERVOMIN, LINSERVOMAX);
    double power = (delta - encoderValue) * 0.08; //amount of power to supply to the servo (as difference from center) e.g 20 power is center (288) + 20 = 308
    if (abs(power) < 5) { //if error is too small, motor will move and will never actually hit desired position
      power = ((power > 0) - (power < 0)) * 5; //assigns power to -10 or 10 (essencially minimum)
    }
    power = constrain(power, -50, 50); //makes sure motor doesn't spin too fast
    int constrainedPwm = LINSERVOMID - power;
    pwm.setPWM(servoNum, 0, constrainedPwm);

    previousEncoderValue = encoderValue;
    encoderValue = readEncoder(servoNum) - encoderValueInitial;

    //handles wrapping: idea is that if previous value is 0, current value is 359, it wrapped around. thus, the true angle is actually -1
    //to accomplish this, we subtract 360 from the encoderValue
    //this way, the next time we read, (say it keeps spinning) and the value read is 300, it will update and say it's -59 or smth like that. 

    if (encoderValue - previousEncoderValue > 90) { //if encoderValue suddenly jumps up a lot, assume it wrapped from 0 to 360; subtract 360 from both values
      encoderValueInitial += 360;
    } else if (encoderValue - previousEncoderValue < -90 ) { //else if encoderValue suddenly jumps down a lot, assume it wrapped from 360 to 0, add 360 to both values
      encoderValueInitial -= 360;
    }

    encoderValue = readEncoder(servoNum) - encoderValueInitial;

    Serial.println(power);
    Serial.println(readEncoder(servoNum));
    Serial.println(encoderValue);
    Serial.println(previousEncoderValue);
    Serial.println(encoderValueInitial);
    Serial.println();

    
    delay(10); //just so it doesn't kill itself
  }
  
  pwm.setPWM(servoNum, 0, LINSERVOMID);
}

int readEncoder(int servoNum) {
  return map(pulseIn(ENCODERPINS[servoNum], HIGH, 3000), 23, 1048, 0, 360);
}

void rotateMove(int servoNum, int angle) {
  int pulseLen = map(angle, 0, 180, SERVOMIN, SERVOMAX);
  pwm.setPWM(servoNum, 0, pulseLen);
}

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
    for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        targetStates[z][x][y] = values[k++];
      }
    }
  }

  Serial.println("State updated.");
}