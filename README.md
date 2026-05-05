# ECG Signal Processing: Simulation & Real-Time Hardware

This project explores the digital signal processing (DSP) life cycle of ECG signals.

## 📁 Task 1a: Simulated Analysis
Processing of MIT-BIH database records to establish a baseline for filtering.
* **Algorithms**: Notch filters (50Hz) and Bandpass filters.
* **Metrics**: Automated calculation of Heart Rate (BPM) and HRV.
* **Results**: [View Task 1a Results](./Results/Task_1a/)

## 📁 Task 1b: Real-Time Hardware
Live data acquisition using an Arduino Nano and AD8232 sensor.
* **Setup**: Physical electrodes connected to an Arduino Nano.
* **Acquisition**: Real-time Serial Communication to MATLAB.
* **Dashboard**: Live PQRST waveform and real-time BPM.
* **Demo**: [Watch Hardware Video](./Results/Task_1b/dsp%20project.mp4)

---

## Technical Results
| Feature | Simulated (1a) | Hardware (1b) |
| :--- | :--- | :--- |
| **Source** | MIT-BIH Database | Live Electrode Input |
| **BPM** | 74.4 BPM | Real-time Variable |
| **Filtering** | Cleaned Powerline Noise | Real-time Artifact Removal |
