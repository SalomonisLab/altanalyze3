from __future__ import annotations

from . import create_app

app = create_app()

if __name__ == "__main__":  # pragma: no cover - manual launch
    app.run(host="0.0.0.0", port=8000, debug=True)
