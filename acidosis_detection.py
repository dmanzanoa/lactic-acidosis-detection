"""
Lactic acidosis detection pipeline.

This module provides a concise, reusable pipeline for analysing blood‑gas data
from the MIMIC‑IV critical care database.  It demonstrates how structured
electronic health records can be queried via Google BigQuery to identify
episodes of lactic or metabolic acidosis that persist for more than two
hours while a patient is receiving vasoactive medications.  It also
illustrates how to cross‑check the resulting cohort against discharge
diagnosis codes.  The script is deliberately self‑contained and suitable
for use in a research or portfolio setting without referring to any
university assignments.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Iterable, List, Tuple
import os

import numpy as np
import pandas as pd

try:
    # Import BigQuery only if available.  This allows the module to be
    # imported on systems without GCP credentials for testing or linting.
    from google.cloud import bigquery
except ImportError:  # pragma: no cover
    bigquery = None  # type: ignore


def run_query(client: bigquery.Client, query: str) -> pd.DataFrame:
    """Execute a SQL query and return the results as a DataFrame.

    Parameters
    ----------
    client : bigquery.Client
        An authenticated BigQuery client instance.
    query : str
        A SQL query string to execute.

    Returns
    -------
    pandas.DataFrame
        The query results loaded into a Pandas DataFrame.

    Raises
    ------
    RuntimeError
        If the BigQuery client library is not installed.
    """
    if bigquery is None:
        raise RuntimeError(
            "google‑cloud‑bigquery is not installed. Please install it with"
            " `pip install google‑cloud‑bigquery` and authenticate with GCP."
        )
    job = client.query(query)
    return job.result().to_dataframe()


def get_lactate_ph_items(client: bigquery.Client) -> Tuple[List[int], List[int]]:
    """Retrieve item identifiers for lactate and arterial pH measurements.

    The MIMIC‑IV database stores laboratory item metadata in
    ``physionet-data.mimiciv_hosp.d_labitems``.  This function returns two
    lists of item identifiers: one for lactate and one for arterial pH.
    """
    lactate_query = """
        SELECT itemid
        FROM `physionet-data.mimiciv_hosp.d_labitems`
        WHERE LOWER(label) LIKE '%lactate%'
    """
    ph_query = """
        SELECT itemid
        FROM `physionet-data.mimiciv_hosp.d_labitems`
        WHERE LOWER(label) LIKE '%ph%' AND LOWER(category) LIKE '%blood%'
    """
    lactate_df = run_query(client, lactate_query)
    ph_df = run_query(client, ph_query)
    return lactate_df["itemid"].tolist(), ph_df["itemid"].tolist()


def fetch_lab_events(
    client: bigquery.Client,
    lactate_items: Iterable[int],
    ph_items: Iterable[int],
) -> pd.DataFrame:
    """Fetch lactate and arterial pH laboratory events from MIMIC‑IV.

    The returned DataFrame contains ``subject_id``, ``hadm_id``, ``charttime``,
    ``itemid`` and a ``measure_value`` column (renamed from ``valuenum``).
    """
    items_str = ",".join(str(i) for i in set(list(lactate_items) + list(ph_items)))
    query = f"""
        SELECT
            le.subject_id,
            le.hadm_id,
            le.charttime,
            le.itemid,
            le.valuenum AS measure_value
        FROM `physionet-data.mimiciv_hosp.labevents` le
        WHERE le.itemid IN ({items_str})
          AND le.valuenum IS NOT NULL
    """
    return run_query(client, query)


def fetch_vasoactive_drugs(client: bigquery.Client) -> pd.DataFrame:
    """Retrieve administration times for common vasoactive medications.

    Vasoactive drugs are extracted from the
    ``physionet-data.mimiciv_hosp.emar`` table by filtering on a curated
    list of medication names.  Only administration records with a known
    schedule time are returned.
    """
    medications = [
        'dopamine', 'milrinone', 'vasopressin', 'nitroglycerin', 'nitroprusside',
        'epinephrine', 'norepinephrine', 'levophed', 'dobutamine', 'hydralazine',
        'labetalol', 'methylene blue', 'terlipressin', 'angiotensin ii',
    ]
    meds_str = ",".join(f"'{m}'" for m in medications)
    query = f"""
        SELECT subject_id, medication, scheduletime
        FROM `physionet-data.mimiciv_hosp.emar`
        WHERE LOWER(medication) IN ({meds_str})
    """
    return run_query(client, query)


def identify_acidosis_episodes(
    bp: pd.DataFrame,
    vadrugs: pd.DataFrame,
    lactate_items: Iterable[int],
    ph_items: Iterable[int],
    min_duration: timedelta = timedelta(minutes=120),
) -> List[int]:
    """Identify subjects with prolonged lactic or metabolic acidosis episodes.

    The lab events are grouped by subject.  Each subject's events are
    sorted chronologically and contiguous periods of abnormal values are
    detected.  An abnormal value is defined as lactate > 4.0 mmol/L or
    arterial pH ≤ 7.35.  If an abnormal episode lasts at least ``min_duration``
    and the patient receives at least one vasoactive drug during that
    interval, the subject is added to the output list.
    """
    candidates: List[int] = []
    for subject_id, data in bp.groupby('subject_id'):
        data = data.sort_values('charttime').reset_index(drop=True)
        episodes: List[Tuple[datetime, datetime]] = []
        start: datetime | None = None
        for _, record in data.iterrows():
            itemid = record['itemid']
            value = record['measure_value']
            is_abnormal = (
                (itemid in lactate_items and value > 4.0) or
                (itemid in ph_items and value <= 7.35)
            )
            if is_abnormal:
                if start is None:
                    start = record['charttime']
            else:
                if start is not None:
                    end = record['charttime']
                    if end - start >= min_duration:
                        episodes.append((start, end))
                    start = None
        subject_drugs = vadrugs[vadrugs['subject_id'] == subject_id]
        for start, end in episodes:
            meds_in_window = subject_drugs[
                (subject_drugs['scheduletime'] > start) &
                (subject_drugs['scheduletime'] < end)
            ]
            if not meds_in_window.empty:
                candidates.append(subject_id)
                break
    return candidates


def evaluate_against_diagnosis(
    client: bigquery.Client,
    subjects: Iterable[int],
    icd_code: str = 'E872',
) -> pd.DataFrame:
    """Compare subject cohort against discharge diagnosis codes.

    Returns a DataFrame with ``subject_id`` and a boolean ``has_code``
    indicating whether the subject has a matching ICD‑10 code in
    ``diagnoses_icd``.
    """
    if not subjects:
        return pd.DataFrame(columns=['subject_id', 'has_code'])
    subj_str = ",".join(str(s) for s in set(subjects))
    diag_query = f"""
        SELECT subject_id, icd_code
        FROM `physionet-data.mimiciv_hosp.diagnoses_icd`
        WHERE subject_id IN ({subj_str}) AND icd_code LIKE '{icd_code}'
    """
    diag_df = run_query(client, diag_query)
    diag_df['has_code'] = True
    all_df = pd.DataFrame({'subject_id': list(set(subjects))})
    result = all_df.merge(diag_df[['subject_id', 'has_code']], how='left', on='subject_id')
    result['has_code'].fillna(False, inplace=True)
    return result


def main() -> None:  # pragma: no cover
    """Command‑line entry point for the acidosis detection pipeline.

    This function authenticates with BigQuery using the
    ``GOOGLE_CLOUD_PROJECT`` environment variable, queries the necessary
    tables, identifies episodes, and prints a summary of how many
    subjects were flagged and how many carried an ICD‑10 code for lactic
    acidosis.
    """
    if bigquery is None:
        raise SystemExit(
            "google‑cloud‑bigquery is not installed. Please install it"
            " and authenticate before running this script."
        )
    project_id = os.environ.get('GOOGLE_CLOUD_PROJECT')
    if not project_id:
        raise SystemExit(
            "The environment variable GOOGLE_CLOUD_PROJECT must be set to your GCP project ID."
        )
    client = bigquery.Client(project=project_id)

    lactate_items, ph_items = get_lactate_ph_items(client)
    lab_events = fetch_lab_events(client, lactate_items, ph_items)
    drugs = fetch_vasoactive_drugs(client)
    subjects = identify_acidosis_episodes(lab_events, drugs, lactate_items, ph_items)
    results = evaluate_against_diagnosis(client, subjects, 'E872')
    total_subjects = len(results)
    num_with_code = results['has_code'].sum()
    accuracy = (num_with_code / total_subjects) * 100 if total_subjects > 0 else 0.0
    print(f"Number of subjects identified: {total_subjects}")
    print(f"Number with lactic acidosis diagnosis: {num_with_code}")
    print(f"Percentage with diagnosis: {accuracy:.2f}%")


if __name__ == '__main__':  # pragma: no cover
    main()
