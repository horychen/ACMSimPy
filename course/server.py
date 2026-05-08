"""
Electric Motor Control - Course Server
Run: python server.py
  - The startup menu auto-discovers lecture*.html,
    homework*.html, homework*_solution.html,
    and codingProject*.html under this folder.
  - Home:       http://<ip>:8000/index.html      (students)
  - Teacher:    http://localhost:8000/teacher.html
"""

import cgi
import http.server
import json
import os
import re
import socket
import threading
import time
from datetime import datetime

PORT = 8000
COURSE_DIR = os.path.dirname(os.path.abspath(__file__))
RESPONSES_DIR = os.path.join(COURSE_DIR, "responses")
LIVE_DIR = os.path.join(COURSE_DIR, "live_responses")
QUIZ_DIR = os.path.join(COURSE_DIR, "quiz_responses")
HOMEWORK_DIR = os.path.join(COURSE_DIR, "homework_submissions")


def build_course_pages():
    pages = []

    if os.path.exists(os.path.join(COURSE_DIR, "index.html")):
        pages.append(("Home", "index.html"))

    discovered = []
    for filename in os.listdir(COURSE_DIR):
        if not filename.endswith(".html") or filename == "index.html":
            continue

        m = re.fullmatch(r"lecture(\d+)\.html", filename)
        if m:
            number = int(m.group(1))
            discovered.append(((10, number, 0), f"Lecture {number}", filename))
            continue

        m = re.fullmatch(r"homework(\d+)(?:(_solution))?\.html", filename)
        if m:
            number = int(m.group(1))
            is_solution = 1 if m.group(2) else 0
            label = f"Homework {number}" + (" Solution" if is_solution else "")
            discovered.append(((20, number, is_solution), label, filename))
            continue

        m = re.fullmatch(r"codingProject(\d+)\.html", filename)
        if m:
            number = int(m.group(1))
            discovered.append(((30, number, 0), f"Coding Project {number}", filename))
            continue

    discovered.sort(key=lambda item: item[0])
    pages.extend((label, filename) for _, label, filename in discovered)

    tail_pages = [
        ("Live Q&A", "live.html"),
        ("Teacher Panel", "teacher.html"),
        ("Q&A Demo", "qa-demo.html"),
    ]
    for label, filename in tail_pages:
        if os.path.exists(os.path.join(COURSE_DIR, filename)):
            pages.append((label, filename))

    return pages


COURSE_PAGES = build_course_pages()

# --- In-memory live Q&A state ---
live_state = {
    "question": None,       # { id, text_zh, text_en, type, options, timestamp }
    "answers": {},           # { student_name: answer_value }
    "question_counter": 0,
    "history": [],           # list of { question, answers } for past questions
}
live_lock = threading.Lock()

# --- Online presence tracking ---
online_users = {}  # { name: last_heartbeat_timestamp }
online_lock = threading.Lock()
ONLINE_TIMEOUT = 15  # seconds

# --- File write locks ---
survey_lock = threading.Lock()
quiz_lock = threading.Lock()


class SurveyHandler(http.server.SimpleHTTPRequestHandler):

    def __init__(self, *args, **kwargs):
        try:
            super().__init__(*args, directory=COURSE_DIR, **kwargs)
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
        elif self.path == "/api/online":
            self._handle_get_online()
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
        elif self.path == "/submit-quiz":
            self._handle_quiz_submit()
        elif self.path == "/submit-homework":
            self._handle_homework_submit()
        elif self.path == "/api/heartbeat":
            self._handle_heartbeat()
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

        with survey_lock:
            os.makedirs(RESPONSES_DIR, exist_ok=True)
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)

        now = datetime.now().strftime("%H:%M:%S")
        print(f"  [{now}] Survey from: {name}")
        self._send_json(200, {"status": "ok"})

    # ---- quiz submit ----
    def _handle_quiz_submit(self):
        try:
            data = self._read_body()
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send_json(400, {"status": "error", "message": "Invalid JSON"})
            return

        name = data.get("name", "").strip()
        if not name:
            self._send_json(400, {"status": "error", "message": "Name is required"})
            return

        lecture = data.get("lecture", "unknown")
        safe_name = name.replace("/", "_").replace("\\", "_").replace("..", "_")
        filepath = os.path.join(QUIZ_DIR, f"{lecture}_{safe_name}.json")

        with quiz_lock:
            os.makedirs(QUIZ_DIR, exist_ok=True)
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)

        now = datetime.now().strftime("%H:%M:%S")
        score = data.get("score", "?")
        total = data.get("total", "?")
        print(f"  [{now}] Quiz {lecture} from: {name} ({score}/{total})")
        self._send_json(200, {"status": "ok"})

    # ---- homework file upload ----
    def _handle_homework_submit(self):
        content_type = self.headers.get("Content-Type", "")
        if "multipart/form-data" not in content_type:
            self._send_json(400, {"status": "error", "message": "Expected multipart/form-data"})
            return

        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={"REQUEST_METHOD": "POST", "CONTENT_TYPE": content_type},
        )

        name = form.getfirst("name", "").strip()
        if not name:
            self._send_json(400, {"status": "error", "message": "Name is required"})
            return

        homework = form.getfirst("homework", "unknown")
        safe_name = name.replace("/", "_").replace("\\", "_").replace("..", "_")
        student_dir = os.path.join(HOMEWORK_DIR, homework, safe_name)
        os.makedirs(student_dir, exist_ok=True)

        files_field = form["files"]
        if not isinstance(files_field, list):
            files_field = [files_field]

        saved = []
        for item in files_field:
            if item.filename:
                safe_fn = os.path.basename(item.filename).replace("..", "_")
                dest = os.path.join(student_dir, safe_fn)
                with open(dest, "wb") as f:
                    f.write(item.file.read())
                saved.append(safe_fn)

        now = datetime.now().strftime("%H:%M:%S")
        print(f"  [{now}] Homework {homework} from: {name} ({len(saved)} file(s))")
        self._send_json(200, {"status": "ok", "files": saved})

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

    # ---- heartbeat: student pings presence ----
    def _handle_heartbeat(self):
        try:
            data = self._read_body()
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send_json(400, {"status": "error"})
            return
        name = data.get("name", "").strip()
        if name:
            with online_lock:
                online_users[name] = time.time()
        self._send_json(200, {"status": "ok"})

    # ---- online: get active user count and names ----
    def _handle_get_online(self):
        now = time.time()
        with online_lock:
            active = {n: t for n, t in online_users.items() if now - t < ONLINE_TIMEOUT}
            # Clean up stale entries
            online_users.clear()
            online_users.update(active)
        self._send_json(200, {
            "count": len(active),
            "users": sorted(active.keys()),
        })

    def log_message(self, format, *args):
        # Suppress noisy GET logs and heartbeats
        if args and "POST" in str(args[0]) and "heartbeat" not in str(args[0]):
            super().log_message(format, *args)


def get_local_ip():
    """Get the machine's real LAN IP, skipping VPN/VMware/virtual adapters."""
    import subprocess
    try:
        # Parse ipconfig to find the WLAN adapter's IPv4 address
        out = subprocess.check_output("ipconfig", encoding="gbk", errors="replace")
        lines = out.splitlines()
        in_wlan = False
        for i, line in enumerate(lines):
            # Match adapter header lines (no leading spaces)
            if not line.startswith(" ") and ("WLAN" in line or "Wi-Fi" in line):
                in_wlan = True
            elif not line.startswith(" ") and line.strip() and in_wlan:
                in_wlan = False  # hit next adapter header
            elif in_wlan and "IPv4" in line:
                ip = line.split(":")[-1].strip()
                if ip and not ip.startswith("127."):
                    return ip
    except Exception:
        pass
    # Fallback: list all IPs, prefer 10.x or 192.168.x, skip 198.18 (VPN)
    try:
        addrs = socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET)
        candidates = sorted(set(a[4][0] for a in addrs))
        for ip in candidates:
            if ip.startswith("10.") or ip.startswith("172.") or ip.startswith("192.168."):
                if not ip.endswith(".1"):  # skip VMware gateway-style .1 addresses
                    return ip
        for ip in candidates:
            if ip.startswith("10.") or ip.startswith("192.168."):
                return ip
    except Exception:
        pass
    return "127.0.0.1"


if __name__ == "__main__":
    os.makedirs(RESPONSES_DIR, exist_ok=True)
    os.makedirs(QUIZ_DIR, exist_ok=True)
    os.makedirs(HOMEWORK_DIR, exist_ok=True)

    ip = get_local_ip()

    print()
    print("=" * 56)
    print("  Electric Motor Control - Course Server")
    print("=" * 56)
    print()
    for label, page in COURSE_PAGES:
        base = "localhost" if page == "teacher.html" else ip
        print(f"  {label:<18} http://{base}:{PORT}/{page}")
    print()

    try:
        import qrcode
        qr = qrcode.QRCode(box_size=1, border=1)
        qr.add_data(f"http://{ip}:{PORT}/index.html")
        qr.make(fit=True)
        qr.print_ascii(invert=True)
        print()
    except ImportError:
        print("  (pip install qrcode for terminal QR code)")
        print()

    print(f"  Survey responses:  {RESPONSES_DIR}")
    print(f"  Quiz responses:    {QUIZ_DIR}")
    print(f"  Homework uploads:  {HOMEWORK_DIR}")
    print("  Press Ctrl+C to stop.")
    print()

    server = http.server.ThreadingHTTPServer(("0.0.0.0", PORT), SurveyHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  Server stopped.")
        server.server_close()
