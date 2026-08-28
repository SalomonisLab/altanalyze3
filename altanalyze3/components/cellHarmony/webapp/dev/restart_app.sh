#!/bin/bash
# Restart the cellHarmony-web server on 127.0.0.1:8000 without racing the old
# process: kill it, wait for the port to be free, start, wait for a fresh 200.
LOG=/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/server_8000.log
pkill -f "uvicorn altanalyze3.components.cellHarmony.webapp.app:app" 2>/dev/null
for i in $(seq 1 60); do
  lsof -nP -iTCP:8000 -sTCP:LISTEN >/dev/null 2>&1 || break
  perl -e 'select undef,undef,undef,0.5'
done
if lsof -nP -iTCP:8000 -sTCP:LISTEN >/dev/null 2>&1; then echo "port 8000 still held"; exit 1; fi
cd /Users/saljh8/Documents/GitHub/altanalyze3
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 nohup /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m uvicorn altanalyze3.components.cellHarmony.webapp.app:app \
  --host 127.0.0.1 --port 8000 --log-level info >> "$LOG" 2>&1 &
for i in $(seq 1 60); do
  code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 3 http://127.0.0.1:8000/ 2>/dev/null)
  if [ "$code" = "200" ]; then echo "server up, pid $(pgrep -f 'uvicorn altanalyze3.components.cellHarmony.webapp.app' | tail -1)"; exit 0; fi
  perl -e 'select undef,undef,undef,0.5'
done
echo "server did not come up"; exit 1
