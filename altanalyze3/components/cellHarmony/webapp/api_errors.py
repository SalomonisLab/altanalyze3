"""User-facing request validation, without framework dumps or submitted values."""
LABELS = {
    'analysis_mode': 'Analysis', 'min_genes': 'Min genes', 'min_counts': 'Min counts',
    'min_cells': 'Min cells', 'mit_percent': 'Mito %', 'align_cutoff': 'Minimum alignment score',
    'existing_state_column': 'Existing cell states / clusters', 'max_k': 'Max K',
    'species': 'Species', 'reference': 'Reference', 'files': 'Input files', 'sample_names': 'Sample names',
}


def validation_message(errors):
    messages = []
    for error in errors:
        field = next((str(v) for v in reversed(error.get('loc', ())) if isinstance(v, str) and v not in {'body', 'query', 'path'}), '')
        label = LABELS.get(field, field.replace('_', ' ').capitalize() or 'A required setting')
        kind = error.get('type', '')
        if field == 'analysis_mode':
            message = 'Choose Supervised, Unsupervised, Both or Use Existing.'
        elif field == 'max_k':
            message = 'Max K must be a whole number of 2 or more, or left blank.'
        elif kind == 'missing':
            message = f'{label} is required.'
        elif kind in {'int_parsing', 'int_type', 'int_from_float'}:
            message = f'{label} must be a whole number.'
        elif kind in {'float_parsing', 'float_type', 'finite_number'}:
            message = f'{label} must be a valid number.'
        elif kind in {'greater_than_equal', 'greater_than', 'less_than_equal', 'less_than'}:
            relation, bound = {'greater_than_equal': ('at least', 'ge'), 'greater_than': ('greater than', 'gt'),
                               'less_than_equal': ('at most', 'le'), 'less_than': ('less than', 'lt')}[kind]
            message = f"{label} must be {relation} {error.get('ctx', {}).get(bound)}."
        else:
            message = f'Check {label.lower()} and try again.'
        if message not in messages:
            messages.append(message)
    return ' '.join(messages) or 'Check the selected settings and try again.'
