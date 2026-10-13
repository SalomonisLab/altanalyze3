"""Local GGUF completion sidecar; model failures cannot terminate the web app.

Run in an environment with llama-cpp-python installed. No model downloads,
dataset reads or inference dependencies are added to the analysis container.
"""
from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import select
import signal
import subprocess
import sys
import threading
import time

MAX_BODY = 256 * 1024


class ModelWorker:
    def __init__(self, options):
        self.options = options
        self.process = None
        self.lock = threading.Lock()
        self.last_used = time.monotonic()

    def stop(self):
        process, self.process = self.process, None
        if process is not None:
            process.kill()
            process.wait(timeout=5)

    def complete(self, payload):
        with self.lock:
            if self.process is None or self.process.poll() is not None:
                self.process = subprocess.Popen(
                    [sys.executable, str(Path(__file__).resolve()), '--worker'],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
            self.last_used = time.monotonic()
            request = {'options': self.options, 'payload': payload}
            process = self.process
            try:
                process.stdin.write(json.dumps(request) + '\n')
                process.stdin.flush()
                if not select.select([process.stdout], [], [], 230)[0]:
                    raise TimeoutError('Local interpretation exceeded generation timeout')
                line = process.stdout.readline()
                if not line:
                    raise RuntimeError('Local interpretation process exited')
                reply = json.loads(line)
                if 'error' in reply:
                    raise RuntimeError(reply['error'])
                return reply
            except Exception:
                self.stop()
                raise
            finally:
                self.last_used = time.monotonic()

    def reap_idle(self):
        # Return model memory after inactivity; the next request reloads lazily.
        while True:
            time.sleep(30)
            with self.lock:
                if time.monotonic() - self.last_used > self.options['idle_seconds']:
                    self.stop()


def worker_main():
    from llama_cpp import Llama
    engine = None
    for line in sys.stdin:
        try:
            request = json.loads(line)
            options, payload = request['options'], request['payload']
            if engine is None:
                engine = Llama(model_path=options['model_path'], n_ctx=options['ctx_size'],
                    n_threads=options['threads'], n_threads_batch=options['threads'],
                    n_gpu_layers=options['gpu_layers'], verbose=False)
            # Independent jobs must not inherit a prior conversation/KV context.
            engine.reset()
            messages = payload.get('messages')
            if not isinstance(messages, list) or not messages:
                raise ValueError('messages must be a nonempty list')
            reply = engine.create_chat_completion(messages=messages,
                temperature=0.1, max_tokens=min(int(payload.get('max_tokens', 3000)), 3000),
                response_format=payload.get('response_format', {'type':'json_object'}))
            reply['model'] = options['model_name']
        except Exception as exc:
            reply = {'error': f'{type(exc).__name__}: {exc}'}
        print(json.dumps(reply), flush=True)


def serve(options):
    worker = ModelWorker(options)

    class Handler(BaseHTTPRequestHandler):
        def send_json(self, status, payload):
            raw = json.dumps(payload).encode()
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

        def do_GET(self):
            if self.path != '/health':
                self.send_json(404, {'error': 'Unknown route'})
                return
            self.send_json(200, {'model': options['model_name'], 'runtime': 'llama-cpp-python',
                'loaded': worker.process is not None and worker.process.poll() is None,
                'ctx_size': options['ctx_size'], 'gpu_layers': options['gpu_layers']})

        def do_POST(self):
            if self.path != '/v1/chat/completions':
                self.send_json(404, {'error': 'Unknown route'})
                return
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if length <= 0 or length > MAX_BODY:
                    self.send_json(413, {'error': 'Request exceeds model context transport budget'})
                    return
                payload = json.loads(self.rfile.read(length))
                if payload.get('model') != options['model_name']:
                    self.send_json(400, {'error': 'Requested model is not installed on this sidecar'})
                    return
                if payload.get('stream'):
                    self.send_json(400, {'error': 'Use asynchronous scALABLE summary polling'})
                    return
                self.send_json(200, worker.complete(payload))
            except (ValueError, TypeError):
                self.send_json(400, {'error': 'Invalid completion request'})
            except Exception as exc:
                self.send_json(503, {'error': f'{type(exc).__name__}: {exc}'})

    threading.Thread(target=worker.reap_idle, daemon=True).start()
    server = ThreadingHTTPServer((options['host'], options['port']), Handler)
    server.daemon_threads = True
    try:
        server.serve_forever()
    finally:
        worker.stop()
        server.server_close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--model', dest='model_path')
    parser.add_argument('--model-name', default='Qwen2.5-7B-Instruct-Q4_K_M')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8086)
    parser.add_argument('--ctx-size', type=int, default=16384)
    parser.add_argument('--threads', type=int, default=6)
    parser.add_argument('--gpu-layers', type=int, default=0)
    parser.add_argument('--idle-seconds', type=int, default=300)
    options = vars(parser.parse_args())
    if options.pop('worker'):
        worker_main()
        return
    if not options['model_path'] or not Path(options['model_path']).is_file():
        parser.error('--model must name an installed GGUF file (first shard for split models)')
    if options['ctx_size'] < 4096 or options['threads'] < 1 or options['idle_seconds'] < 30:
        parser.error('ctx-size >= 4096, threads >= 1 and idle-seconds >= 30 are required')
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
    serve(options)


if __name__ == '__main__':
    main()
