// AD8232 ECG → Arduino Uno/Nano → Serial for MATLAB
// Sampling rate: 360 Hz

const int ECG_PIN   = A0;    // Analog input
const int LO_PLUS   = 10;    // Leads-off Detection +
const int LO_MINUS  = 11;    // Leads-off Detection -

const unsigned long SAMPLE_INTERVAL_US = 1000000UL / 360;  // ~2777 microseconds
unsigned long lastSample = 0;

void setup() {
  // Use a high baud rate for real-time DSP to avoid buffer lag
  Serial.begin(115200);
  
  pinMode(LO_PLUS,  INPUT);
  pinMode(LO_MINUS, INPUT);
  
  // Note: Arduino Uno/Nano resolution is fixed at 10-bit (0-1023)
  // unlike the ESP32's 12-bit (0-4095).
}

void loop() {
  unsigned long now = micros();
  if (now - lastSample >= SAMPLE_INTERVAL_US) {
    lastSample = now;

    // Leads-off detection: if either LO pin is HIGH, electrode disconnected
    if (digitalRead(LO_PLUS) == 1 || digitalRead(LO_MINUS) == 1) {
      Serial.println(0);
    } else {
      Serial.println(analogRead(ECG_PIN));
    }
  }
}
