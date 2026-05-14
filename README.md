# ECG Signal Processing: Simulation & Real-Time Hardware
**National University of Sciences and Technology (NUST)** **School of Electrical Engineering and Computer Science (SEECS)**

This project documents the design and implementation of a Digital Signal Processing (DSP) system for ECG signals, moving from database simulations to a live hardware acquisition dashboard.

---

## 📁 Project Structure

### [Task 1a: Simulated Analysis](./Scripts/Task_1a/)
* **Source:** MIT-BIH Arrhythmia Database.
* **Focus:** Offline processing, 50Hz Notch filtering, and automated BPM/HRV calculation.
* **Results:** [View Simulation Plots](./Results/Task_1a/)

### [Task 1b: Hardware Verification](./Scripts/Task_1b/)
* **Setup:** Initial integration of Arduino Nano and AD8232 ECG sensor.
* **Goal:** Verifying data transmission and basic waveform integrity.
* **Media:** [Watch Setup Video](./Results/Task_1b/)

### [Task 2: Advanced Real-Time Streams](./Scripts/Task_2/)
We implemented three parallel DSP approaches to optimize real-time performance:
* **Stream A**: Standard acquisition with fixed-window filtering.
* **Stream B**: Optimized R-Peak detection and dynamic BPM calculation.
* **Stream C**: Advanced Heart Rate Variability (HRV) analysis and a live Dashboard UI.
* **Visuals:** [Compare Stream Outputs](./Results/Task_2/)

---

## 🛠️ Hardware Specifications
* **Microcontroller**: Arduino Uno
* **Sensor**: AD8232 Heart Rate Monitor Module
* **Electrodes**: 3-Lead ECG System (RA, LA, LL)
* **Interface**: MATLAB Serial Communication (9600 Baud)

## 📊 Technical Highlights
| Feature | Simulation (1a) | Real-Time (2) |
| :--- | :--- | :--- |
| **Data Source** | MIT-BIH Records | Live Human Input |
| **Noise Filtering** | Digital IIR Notch | Real-Time Hybrid Filter |
| **Visualization** | Static Subplots | Live Scrolling Dashboard |
| **Heart Rate** | Automated Batch | Dynamic Real-Time Update |

---

## 🚀 How to Run
1. **Hardware Connection**: Connect the AD8232 to the Arduino Nano and attach electrodes.
2. **Arduino**: Upload the `.ino` sketch found in the `Scripts/Task_2` sub-folders.
3. **MATLAB**: Run the `StreamC_RealTime_DSP.m` script to launch the live dashboard.

---
*DSP Project Group 3 - SEECS NUST*
