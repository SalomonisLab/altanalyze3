from fastapi.testclient import TestClient
from altanalyze3.components.cellHarmony.webapp.app import create_app


def test_invalid_upload_mode_has_readable_validation_message(tmp_path):
    app=create_app({'JOB_STORAGE':str(tmp_path),'ISOLATE_JOBS':False})
    with TestClient(app) as client:
        response=client.post('/api/jobs',data={'species':'human','reference':'fixture','analysis_mode':'unsupported',
            'sample_names':'one'},files={'files':('input.h5ad',b'not read','application/octet-stream')})
        assert response.status_code==422
        assert response.json()['detail']=='Choose Supervised, Unsupervised, Both or Use Existing.'
        assert 'string_pattern_mismatch' not in response.text and 'File "' not in response.text
    app.state.job_runner.executor.shutdown(wait=True)


def test_invalid_qc_does_not_expose_submitted_values(tmp_path):
    app=create_app({'JOB_STORAGE':str(tmp_path),'ISOLATE_JOBS':False})
    with TestClient(app) as client:
        response=client.post('/api/jobs/unused/qc',json={'min_genes':'/private/invalid_value','max_k':1})
        assert response.status_code==422
        message=response.json()['detail']
        assert 'Min genes must be a whole number.' in message
        assert 'Max K must be a whole number of 2 or more, or left blank.' in message
        assert '/private' not in response.text
    app.state.job_runner.executor.shutdown(wait=True)


def test_existing_h5_rejected_before_creating_or_cancelling_jobs(tmp_path,monkeypatch):
    app=create_app({'JOB_STORAGE':str(tmp_path),'ISOLATE_JOBS':False})
    def forbidden(*args): raise AssertionError('invalid format must not cancel an existing job')
    monkeypatch.setattr(app.state.job_runner,'cancel_job',forbidden)
    with TestClient(app) as client:
        response=client.post('/api/jobs',data={'species':'human','analysis_mode':'existing',
            'sample_names':'one','previous_job_id':'a'*32},files={'files':('BOS.h5',b'unread counts','application/octet-stream')})
        assert response.status_code==400
        assert 'saved clusters or cell types' in response.json()['detail']
        assert 'Supervised, Unsupervised or Both' in response.json()['detail']
        assert not list(tmp_path.glob('*/job.json'))
    app.state.job_runner.executor.shutdown(wait=True)
