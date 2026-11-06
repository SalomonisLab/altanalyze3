"""Launch the included web app and optional GGUF model service together.

After installing AltAnalyze3[web,biology-summary], pass an installed model:
python -m altanalyze3.components.cellHarmony.webapp.local --model /path/model.gguf
Weights are never downloaded implicitly. Model environment settings are passed
only to this web process; no machine-specific paths are written to the repo.
"""
import argparse
import importlib.util
import os
import platform
import signal
import socket
import subprocess
import sys
import time
import urllib.request
from .biology_model import validate_model, __file__ as model_script



def model_environment(port, model_name='Qwen2.5-7B-Instruct-Q4_K_M'):
    return {'SCALABLE_SUMMARY_URL': f'http://127.0.0.1:{port}/v1/chat/completions',
            'SCALABLE_SUMMARY_MODEL': model_name,
            'SCALABLE_SUMMARY_PROTOCOL': 'openai'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True)
    parser.add_argument('--model-name', default='Qwen2.5-7B-Instruct-Q4_K_M')
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('--model-port', type=int, default=8086)
    parser.add_argument('--threads', type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument('--gpu-layers', type=int, default=-1 if platform.system() == 'Darwin' else 0)
    args = parser.parse_args()
    if not importlib.util.find_spec('llama_cpp'):
        parser.error('Model runtime is missing. Install the optional biology-summary dependencies; see README Local setup.')
    try:
        model = validate_model(args.model)
    except ValueError as exc:
        parser.error(str(exc))
    if args.port == args.model_port or not all(1 <= p <= 65535 for p in (args.port, args.model_port)) or args.threads < 1:
        parser.error('Use distinct valid ports and threads >= 1.')
    for port in (args.port, args.model_port):
        with socket.socket() as probe:
            try:
                probe.bind(('127.0.0.1', port))
            except OSError:
                parser.error(f'Port {port} is already in use. Choose another port or stop the existing service.')
    model_process = subprocess.Popen([sys.executable, model_script,
        '--model', str(model), '--model-name', args.model_name, '--port', str(args.model_port), '--threads', str(args.threads),
        '--gpu-layers', str(args.gpu_layers)])
    web_process = None
    def stop(*_):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop)
    try:
        deadline = time.monotonic() + 30
        while True:
            if model_process.poll() is not None:
                raise RuntimeError('Model service did not start; check the model path, runtime and port above.')
            try:
                with urllib.request.urlopen(f'http://127.0.0.1:{args.model_port}/health', timeout=1):
                    break
            except OSError:
                if time.monotonic() > deadline:
                    raise RuntimeError('Model service startup timed out.') from None
                time.sleep(.2)
        env = dict(os.environ, **model_environment(args.model_port, args.model_name))
        web_process = subprocess.Popen([sys.executable, '-m', 'uvicorn',
            'altanalyze3.components.cellHarmony.webapp.app:app', '--host', '127.0.0.1', '--port', str(args.port)], env=env)
        web_process.wait()
    except KeyboardInterrupt:
        pass
    finally:
        for process in (web_process, model_process):
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill(); process.wait()


if __name__ == '__main__':
    main()
