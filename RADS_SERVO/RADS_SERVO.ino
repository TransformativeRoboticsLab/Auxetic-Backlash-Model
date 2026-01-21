#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>
#include <ServoInput.h>
const int X = 3; // X axis
const int Y = 3; // Y axis
const int SERVONUM = 2; // linear or expanding servo; 0 is linear, 1 is expanding

ServoInputPin<1> encoder0;
ServoInputPin<2> encoder2;
ServoInputPin<3> encoder3;
ServoInputPin<4> encoder4;
ServoInputPin<5> encoder5;
ServoInputPin<6> encoder6;
ServoInputPin<7> encoder7;
ServoInputPin<8> encoder8;

ServoInputPinBase* encoders[] = {&encoder0, &encoder1, &encoder2, &encoder3, &encoder4, &encoder5, &encoder6 , &encoder7, &encoder8}

int targetStates[SERVONUM][X][Y];
int linCurrentStates[X][Y]; //keeps track of the current states of the linear actuator

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();

#define SERVOMIN  205 // Minimum pulse length count (out of 4096)
#define SERVOMAX  410 // Maximum pulse length count (out of 4096)
#define LINSERVOMIN 200 //Min pulse len (conservative)
#define LINSERVOMID 288 //point of linear servo middle
#define LINSERVOMAX 350 // Max pulse len (conservative)
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

  for (auto* s : inputs) s->attach();

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
  //updateServoPositions();
  pwm.setPWM(1, 0, 288); //292 stops, 284 stops; call it 288
/*

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
*/
}

void updateServoPositions(){
  int servoIndex = 0;
  //updating linear servos
  for (int x = 0; x < X; x++) {
      for (int y = 0; y < Y; y++) {
        int currentState = linCurrentStates[x][y];
        int targetState = targetStates[0][x][y];
        int delta = targetState - currentState;
        linearMove(servoIndex, delta);
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

void linearMove(int servoNum, int delta) {

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