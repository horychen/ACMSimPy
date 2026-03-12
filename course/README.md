# ACMSimPy Course Platform

A local web-based teaching toolkit for the graduate-level Electric Motor Control course. No internet required — runs entirely on your laptop over campus WiFi.

## Quick Start

```bash
cd course
python server.py
```

The server prints three URLs:

```
Survey (students): http://192.168.x.x:8000/lecture1.html
Live   (students): http://192.168.x.x:8000/live.html
Teacher panel:     http://localhost:8000/teacher.html
```

Share the student URL with the class (write on board, project QR code, etc).

## Pages

| Page | URL | Who | Purpose |
|------|-----|-----|---------|
| **Landing** | `/index.html` | Anyone | Animated navigation page with course tour |
| **Pre-Survey** | `/lecture1.html` | Students | Background assessment before first lecture |
| **Live Q&A** | `/live.html` | Students | Answer real-time questions during lecture |
| **Teacher Panel** | `/teacher.html` | Instructor | Push questions, view stats, monitor online users |
| **Q&A Demo** | `/qa-demo.html` | Anyone | Interactive interview-style demo (original index.html) |

## Typical Lecture Flow

### Before Class: Pre-Survey

1. Run `python server.py`
2. Students open `lecture1.html` on their phones (same WiFi network)
3. They fill out 14 questions across 3 sections (About You / Conceptual / Expectations)
4. Responses are saved to `responses/<name>.json`
5. Feed all JSONs to Claude for class profile analysis

### During Class: Live Q&A

1. Open `teacher.html` on your laptop
2. Students open `live.html` on their phones and enter their name
3. Compose a question — pick a type:
   - **Open-ended**: students type free text
   - **Single choice**: students pick one option
   - **Multi-choice**: students pick multiple options
4. Use quick templates: `Yes/No`, `A/B/C/D`, `1-5 Scale`
5. Click **Push to students** — the question appears on all phones instantly
6. Watch the **Statistics** tab update live:
   - Choice questions: bar chart with vote counts and percentages
   - Text questions: response count, average length, frequent terms
7. Switch to **All Answers** tab to see individual responses
8. The **online badge** (top-right) shows how many students are connected — click it to see names
9. Push the next question when ready (previous one is auto-archived)

### After Class: Save & Analyze

1. Click **Save session to disk** in teacher panel
2. Session history is written to `live_responses/session_<timestamp>.json`
3. Feed the JSON to Claude with a prompt like:

```
Analyze these live Q&A responses for my electric motor control class.
Identify misconceptions, knowledge gaps, and engagement patterns.
Suggest what to cover in the next lecture.
```

## File Structure

```
course/
  server.py              # Python HTTP server (all endpoints)
  index.html             # Landing page with animated tour
  lecture1.html           # Pre-survey questionnaire (mobile-first)
  live.html              # Student live Q&A page (mobile-first)
  teacher.html           # Teacher control panel with statistics
  qa-demo.html           # Original Q&A interview demo
  responses/             # Pre-survey JSON files (one per student)
  live_responses/        # Saved live session JSON files
  spec-lecture1-presurvey.md  # Design spec for the pre-survey
```

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| POST | `/submit` | Submit pre-survey response |
| POST | `/api/question` | Teacher pushes a new question |
| GET | `/api/question` | Students poll for current question |
| POST | `/api/answer` | Student submits an answer |
| GET | `/api/responses` | Teacher gets all answers for current question |
| POST | `/api/clear` | Clear current question (archives it) |
| GET | `/api/history` | Get all past questions and answers |
| POST | `/api/save-history` | Save session to disk |
| POST | `/api/heartbeat` | Student presence ping (every 5s) |
| GET | `/api/online` | Get online user count and names |

## Requirements

- Python 3.8+ (uses only stdlib: `http.server`, `json`, `socket`, `threading`)
- Optional: `pip install qrcode` for terminal QR code on startup

## Network Setup

Both your laptop and students' phones must be on the **same WiFi network**.

If phones can't connect, check Windows Firewall:

```powershell
# Run in admin PowerShell
New-NetFirewallRule -DisplayName "Course Server" -Direction Inbound -LocalPort 8000 -Protocol TCP -Action Allow
```

## Notes

- **No internet needed** — all assets are inline, no CDN dependencies
- **Mobile-first** — all student pages are optimized for phone screens
- **Duplicate submissions** — pre-survey overwrites by student name; live Q&A overwrites by name per question
- **localStorage** — pre-survey saves progress locally; live Q&A remembers student name
- **Online tracking** — students are counted as "online" if they pinged within the last 15 seconds
