"""
Lecture 1 Pre-Survey + Live Q&A Server
Run: python server.py
  - Survey:  http://<ip>:8000/lecture1.html  (students)
  - Live:    http://<ip>:8000/live.html      (students)
  - Teacher: http://localhost:8000/teacher.html
"""

import http.server
import json
import os
import socket
import threading
import time
from datetime import datetime

PORT = 8000
RESPONSES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "responses")
LIVE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "live_responses")

# --- In-memory live Q&A state ---
live_state = {
    "question": None,       # { id, text_zh, text_en, type, options, timestamp }
    "answers": {},           # { student_name: answer_value }
    "question_counter": 0,
    "history": [],           # list of { question, answers } for past questions
}
live_lock = threading.Lock()


class SurveyHandler(http.server.SimpleHTTPRequestHandler):

    def __init__(self, *args, **kwargs):
        try:
            super().__init__(*args, directory=os.path.dirname(os.path.abspath(__file__)), **kwargs)
        except ConnectionResetError:
            pass

    def _send_json(self, code, data):
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(body)

    def _read_body(self):
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length)
        return json.loads(raw.decode("utf-8"))

    # ---- routing ----
    def do_GET(self):
        if self.path == "/api/question":
            self._handle_get_question()
        elif self.path == "/api/responses":
            self._handle_get_responses()
        elif self.path == "/api/history":
            self._handle_get_history()
        else:
            super().do_GET()

    def do_POST(self):
        if self.path == "/submit":
            self._handle_survey_submit()
        elif self.path == "/api/question":
            self._handle_push_question()
        elif self.path == "/api/answer":
            self._handle_live_answer()
        elif self.path == "/api/clear":
            self._handle_clear_question()
        elif self.path == "/api/save-history":
            self._handle_save_history()
        else:
            self._send_json(404, {"status": "error", "message": "Not found"})

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    # ---- survey submit (unchanged) ----
    def _handle_survey_submit(self):
        try:
            data = self._read_body()
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send_json(400, {"status": "error", "message": "Invalid JSON"})
            return

        name = data.get("name", "").strip()
        if not name:
            self._send_json(400, {"status": "error", "message": "Name is required"})
            return

        safe_name = name.replace("/", "_").replace("\\", "_").replace("..", "_")
        filepath = os.path.join(RESPONSES_DIR, f"{safe_name}.json")

        os.makedirs(RESPONSES_DIR, exist_ok=True)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)

        now = datetime.now().strftime("%H:%M:%S")
        print(f"  [{now}] Survey from: {name}")
        self._send_json(200, {"status": "ok"})

    # ---- live Q&A: teacher pushes a question ----
    def _handle_push_question(self):
        try:
            data = self._read_body()
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send_json(400, {"status": "error", "message": "Invalid JSON"})
            return

        with live_lock:
            # Archive current question if it has answers
            if live_state["question"] and live_state["answers"]:
                live_state["history"].append({
                    "question": live_state["question"],
                    "answers": dict(live_state["answers"]),
                })

            live_state["question_counter"] += 1
            live_state["question"] = {
                "id": live_state["question_counter"],
                "text_zh": data.get("text_zh", ""),
                "text_en": data.get("text_en", ""),
                "type": data.get("type", "textarea"),        # textarea | single-choice | multi-choice
                "options": data.get("options", []),
                "timestamp": datetime.now().isoformat(),
            }
            live_state["answers"] = {}

        now = datetime.now().strftime("%H:%M:%S")
        print(f"  [{now}] Live Q#{live_state['question_counter']}: {data.get('text_zh', '')[:40]}")
        self._send_json(200, {"status": "ok", "id": live_state["question_counter"]})

    # ---- live Q&A: student polls current question ----
    def _handle_get_question(self):
        with live_lock:
            q = live_state["question"]
            count = len(live_state["answers"])
        if q:
            self._send_json(200, {"question": q, "answer_count": count})
        else:
            self._send_json(200, {"question": None, "answer_count": 0})

    # ---- live Q&A: student submits answer ----
    def _handle_live_answer(self):
        try:
            data = self._read_body()
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send_json(400, {"status": "error", "message": "Invalid JSON"})
            return

        name = data.get("name", "").strip()
        answer = data.get("answer", "")
        qid = data.get("question_id")

        if not name:
            self._send_json(400, {"status": "error", "message": "Name is required"})
            return

        with live_lock:
            if not live_state["question"] or live_state["question"]["id"] != qid:
                self._send_json(410, {"status": "error", "message": "Question has changed"})
                return
            live_state["answers"][name] = answer

        now = datetime.now().strftime("%H:%M:%S")
        print(f"  [{now}] Live answer from: {name}")
        self._send_json(200, {"status": "ok"})

    # ---- live Q&A: teacher polls all answers ----
    def _handle_get_responses(self):
        with live_lock:
            q = live_state["question"]
            answers = dict(live_state["answers"])
        self._send_json(200, {"question": q, "answers": answers})

    # ---- live Q&A: teacher clears current question ----
    def _handle_clear_question(self):
        with live_lock:
            if live_state["question"] and live_state["answers"]:
                live_state["history"].append({
                    "question": live_state["question"],
                    "answers": dict(live_state["answers"]),
                })
            live_state["question"] = None
            live_state["answers"] = {}
        self._send_json(200, {"status": "ok"})

    # ---- live Q&A: get history ----
    def _handle_get_history(self):
        with live_lock:
            history = list(live_state["history"])
        self._send_json(200, {"history": history})

    # ---- live Q&A: save all history to disk ----
    def _handle_save_history(self):
        with live_lock:
            # Include current question if active
            history = list(live_state["history"])
            if live_state["question"]:
                history.append({
                    "question": live_state["question"],
                    "answers": dict(live_state["answers"]),
                })

        os.makedirs(LIVE_DIR, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        filepath = os.path.join(LIVE_DIR, f"session_{ts}.json")
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(history, f, ensure_ascii=False, indent=2)

        print(f"  Live session saved to: {filepath}")
        self._send_json(200, {"status": "ok", "file": filepath})

    def log_message(self, format, *args):
        # Suppress noisy GET logs
        if args and "POST" in str(args[0]):
            super().log_message(format, *args)


def get_local_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"


if __name__ == "__main__":
    os.makedirs(RESPONSES_DIR, exist_ok=True)

    ip = get_local_ip()

    print()
    print("=" * 56)
    print("  Electric Motor Control - Course Server")
    print("=" * 56)
    print()
    print(f"  Survey (students): http://{ip}:{PORT}/lecture1.html")
    print(f"  Live   (students): http://{ip}:{PORT}/live.html")
    print(f"  Teacher panel:     http://localhost:{PORT}/teacher.html")
    print()

    try:
        import qrcode
        qr = qrcode.QRCode(box_size=1, border=1)
        qr.add_data(f"http://{ip}:{PORT}/live.html")
        qr.make(fit=True)
        qr.print_ascii(invert=True)
        print()
    except ImportError:
        print("  (pip install qrcode for terminal QR code)")
        print()

    print(f"  Responses: {RESPONSES_DIR}")
    print("  Press Ctrl+C to stop.")
    print()

    server = http.server.HTTPServer(("0.0.0.0", PORT), SurveyHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  Server stopped.")
        server.server_close()
