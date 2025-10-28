# Lactic Acidosis Detection in Critical Care

This repository contains a Python script for detecting episodes of lactic or metabolic acidosis in critically ill patients using the MIMIC IV electronic health records database. The goal of the project is to demonstrate how to query large healthcare datasets, define rule‑based phenotypes using lab results and medication administration times, and evaluate those phenotypes against diagnostic codes.

## Overview

The included script (`acidosis_detection.py`) connects to Google BigQuery to retrieve MIMIC IV hospital lab events and medication administration records. It identifies episodes where lactate levels exceed 4.0 mmol/L or blood pH drops below 7.35 and checks whether these abnormalities persist for more than 120 minutes. Episodes are then cross‑referenced with the administration times of vasoactive drugs (e.g., vasopressin, dopamine, norepinephrine) to determine whether the patient received therapy during the abnormal period. The script also looks up ICD‑10 codes for lactic acidosis (e.g., E872) to allow evaluation of how many phenotyped cases align with charted diagnoses.

By automating this pipeline you can build cohorts for downstream modelling, monitor adherence to treatment protocols, or explore how soon lactic acidosis resolves after intervention. It is designed for portfolio demonstration rather than clinical use.

## Installation

1. Clone this repository and navigate into its folder.
2. Install the Python dependencies:
   ```
   pip install -r requirements.txt
   ```
3. Ensure you have access to the [MIMIC IV](https://physionet.org/content/mimiciv/) dataset via Google BigQuery and set the `GOOGLE_CLOUD_PROJECT` environment variable to your own project ID:
   ```
   export GOOGLE_CLOUD_PROJECT=<your-bigquery-project-id>
   ```
   Authentication to BigQuery can be handled by the `google-cloud-bigquery` library (e.g., via application default credentials).

## Usage

Run the script from the command line:
```
python acidosis_detection.py
```

By default it will:
- fetch lab events for lactate and pH measurements,
- retrieve vasoactive drug administration times,
- identify prolonged abnormal episodes,
- evaluate subjects against ICD‑10 diagnoses for lactic acidosis,
- print summary statistics on the number of identified subjects and the proportion with documented lactic acidosis.

You can customise the query functions or criteria by editing `acidosis_detection.py`.

## Project Structure

- `acidosis_detection.py` – core script implementing the BigQuery queries, episode detection and evaluation logic.
- `requirements.txt` – list of required Python packages.
- `README.md` – this document.

## Disclaimer

This project is intended for educational and demonstrative purposes. It is not a substitute for clinical judgement and should not be used to make decisions about patient care.
