"""Authorized read-only audit of ten all-missing targets against original assays.

Does not fill values, fit models, select replacement molecules, or change sources.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import openpyxl
import pandas as pd

ROOT = Path('/Users/saljh8/Dropbox/Collaborations/Grimes/Human-MS-impute')
WORKBOOK = Path('/Users/saljh8/Downloads/43018_2026_1175_MOESM3_ESM.xlsx')
A3 = Path(__file__).resolve().parents[3]


def main():
    state = json.loads((A3 / 'rna2lipid/integrity/decision_state.json').read_text())
    question = next(q for q in state['questions'] if q['id'] == 'AML_ten_all_missing_target_sources')
    if not question.get('further_source_audit_authorized') or not question.get('user_answer'):
        raise RuntimeError('Original-assay diagnostic requires actual user authorization')
    targets = question['target_ids']
    spec = importlib.util.spec_from_file_location('original_ms_extraction', ROOT / 'code/build_unique_ms_tables.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    report = {'authorization': question['user_answer'], 'targets': {t: [] for t in targets},
              'workbook': str(WORKBOOK), 'workbook_sha256': hashlib.sha256(WORKBOOK.read_bytes()).hexdigest(),
              'sources_changed': False, 'targets_excluded': False, 'fit_executed': False}
    for cached in (True, False):
        wb = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=cached)
        for sheet in ('Table21', 'Table22'):
            rows = wb[sheet].iter_rows(values_only=True)
            headers = [str(v) if v is not None else f'_c{i}' for i, v in enumerate(next(rows))]
            name_i = headers.index('Metabolite')
            case_i = [i for i, header in enumerate(headers) if original.norm(header)]
            for row_n, row in enumerate(rows, start=2):
                # Sparse XLSX rows omit trailing blank cells; represent those
                # absent cells as None, as the original dataframe reader does.
                source_row_width = len(row)
                row = tuple(row) + (None,) * max(0, len(headers) - len(row))
                target = str(row[name_i]).strip()
                if target not in report['targets']:
                    continue
                cells = [row[i] for i in case_i]
                if cached:
                    numeric = pd.to_numeric(pd.Series(cells), errors='coerce')
                    record = {'sheet': sheet, 'excel_row': row_n,
                              'source_row_width': source_row_width,
                              'source_name': row[name_i], 'case_columns': len(case_i),
                              'normalized_case_ids': [original.norm(headers[i]) for i in case_i],
                              'numeric_case_values': int(numeric.notna().sum()),
                              'nonblank_case_cells': [{'column': headers[i], 'value': str(row[i])}
                                                      for i in case_i if row[i] is not None],
                              'noncase_metadata': {headers[i]: str(value) for i, value in enumerate(row)
                                                   if i not in case_i and value is not None}}
                    report['targets'][target].append(record)
                else:
                    matching = next(r for r in report['targets'][target]
                                    if r['sheet'] == sheet and r['excel_row'] == row_n)
                    matching['formula_case_cells'] = [{'column': headers[i], 'formula': row[i]}
                                                     for i in case_i if isinstance(row[i], str)
                                                     and row[i].startswith('=')]
        wb.close()
    for target, records in report['targets'].items():
        if not records:
            raise RuntimeError(f'Required target absent from original assay rows: {target}')
    report['recoverable_numeric_observations'] = sum(r['numeric_case_values']
                                                  for records in report['targets'].values() for r in records)
    out = Path(__file__).with_name('all_missing_assay_row_audit.json')
    out.write_text(json.dumps(report, indent=2) + '\n')
    for target, records in report['targets'].items():
        print(target, [(r['sheet'], r['excel_row'], r['numeric_case_values'],
                        len(r['formula_case_cells'])) for r in records])
    print('Recoverable numeric observations:', report['recoverable_numeric_observations'])


if __name__ == '__main__':
    main()
