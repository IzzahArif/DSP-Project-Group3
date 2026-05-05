ECG Signal Processing: From Simulation to Real-Time Hardware
This project explores the digital signal processing (DSP) life cycle of Electrocardiogram (ECG) signals, starting with offline analysis of medical databases and progressing to a real-time acquisition system using an Arduino Nano and AD8232 sensor.

Project Structure
📂 Task_1a: Simulated Analysis
Processing of MIT-BIH database records to establish a baseline for filtering and peak detection.

Algorithms: Implemented Notch filters (50Hz) and Bandpass filters to remove powerline noise and baseline wander.

Metrics: Automated calculation of Heart Rate (BPM) and Heart Rate Variability (HRV).

Results: View Simulated Results.

📂 Task_1b: Real-Time Hardware Acquisition
Live data acquisition and processing using physical hardware.

Hardware Setup: Arduino Nano integrated with an AD8232 ECG sensor module.

Acquisition: Real-time transmission of analog signals via Serial Communication to MATLAB.

Visualization: A live dashboard displaying the refined PQRST waveform and real-time BPM.

Media: Watch the Hardware Demo.
