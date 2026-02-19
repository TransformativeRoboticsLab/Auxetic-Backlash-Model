#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();

#define SERVOMIN  150 // Minimum pulse length count (out of 4096)
#define SERVOMAX  600 // Maximum pulse length count (out of 4096)
#define SERVO_FREQ 50 // Analog servos run at ~50 Hz updates
#define NUM_SERVOS 16 // Total number of servo channels

void setup() {
  Serial.begin(9600);
  Serial.println("16 channel simultaneous servo test!");

  pwm.begin();
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);

  delay(10);
}

void loop() {
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