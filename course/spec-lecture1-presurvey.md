# Lecture 1 Pre-Survey: Electric Motor Control Baseline Assessment

## Overview

A mobile-friendly, single-page web app served via Python's `http.server`. Students scan a QR code (or enter a short URL) on their phones to answer a series of questions about electric machines and control theory. Answers are collected as JSON files on the instructor's laptop for AI-assisted analysis before the next lecture.

**Class size**: Under 20 graduate students.

## Constraints

- **Mobile-first**: Most students will use their phones on day one — no laptop required.
  - Large touch targets (min 44px), readable font sizes (min 16px to prevent iOS zoom), no hover-dependent interactions.
  - Viewport fits phone screens (360px-430px width). Tablet and desktop are secondary.
  - Soft keyboard must not obscure the current input — auto-scroll into view.
- **No internet required**: The instructor's laptop runs a local Python HTTP server on the campus WiFi. All assets are inline (no CDN fonts, no external dependencies).
- **No build step**: Pure HTML + CSS + vanilla JS in a single file, served as-is.
- **No database**: Each submission is POSTed as JSON to a tiny Python handler that writes one `.json` file per student into a `responses/` folder.

## Architecture

```
course/
  server.py            # Python 3 http.server + POST handler
  lecture1.html         # The questionnaire (single file, all inline)
  responses/            # Created automatically, one JSON per submission
    张三.json
```

### server.py

- Extend `http.server.SimpleHTTPRequestHandler`.
- Handle `POST /submit`:
  - Parse JSON body.
  - Write to `responses/` using student name as filename (e.g., `张三.json`).
  - **Overwrite on duplicate**: if the same name submits again, the latest submission replaces the previous file.
  - Respond `200 OK` with `{ "status": "ok" }`.
- Serve static files from `course/` directory.
- Print each submission to console as confirmation (name + timestamp).
- Bind to `0.0.0.0` so phones on the same network can access it.
- On startup, print the local IP address and port, e.g., `http://192.168.1.42:8000/lecture1.html`.

### lecture1.html

Adopt the visual style of `index.html` (dark theme, amber accent, card-based flow, progress bar) with these modifications:

- **No Google Fonts** — use system font stack: `-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif` and `Georgia, serif` for display.
- **Touch-optimized**: larger buttons, more padding, `font-size: 16px` on all inputs.
- **Bilingual**: Questions in Chinese (primary) with English subtitle in smaller/muted text.
- **localStorage auto-save**: Each section completion saves answers to `localStorage`. If the student reloads or loses connection, they resume where they left off.

## UX Flow: Hybrid Section Cards

Instead of one-question-per-screen (too many taps) or a single scrollable form (too mundane), use **3 animated section cards**. Each card contains a group of related questions. Transitions between cards use the same `fadeUp` animation from `index.html`.

```
[========>              ] Section 2 of 3   67%

  Section 1 complete

  +------------------------------------+
  | PART 2: CONCEPTUAL BASELINE        |
  |                                    |
  | Q04: IM vs PMSM difference?        |
  | [textarea                       ]  |
  |       还不了解 / I don't know yet   |
  |                                    |
  | Q05: How is torque produced?       |
  | [textarea                       ]  |
  |       还不了解 / I don't know yet   |
  |                                    |
  | ...                                |
  |                                    |
  |              [Back]  [Continue]     |
  +------------------------------------+
```

### Progress Bar

- Shows "Section X of 3" and percentage.
- Fills proportionally: Section 1 = 0-33%, Section 2 = 33-67%, Section 3 = 67-100%.
- Uses amber gradient with glow, same as `index.html`.

### Section Transitions

- When "Continue" is tapped, the current card fades out + slides up, next card fades in.
- A compact summary of completed sections appears above the current card (collapsed, showing section name + checkmark).

## Questions

Total: 14 questions across 3 sections.

### Section 1: About You (4 questions)

| #  | Text (ZH) | Text (EN) | Type | Options / Notes |
|----|-----------|-----------|------|-----------------|
| 00 | 你的姓名 | Your name | text input | **Required**. Used as filename for response JSON. |
| 01 | 你之前接触过电机控制吗？ | Have you had prior exposure to motor control? | single-choice | 从未 / Never, 课程学过 / Studied in class, 做过项目 / Hands-on project, 工作经验 / Industry experience |
| 02 | 你用过哪些编程语言？ | Which programming languages have you used? | multi-choice | Python, MATLAB, C/C++, Julia, 其他 / Other |
| 03 | 你对以下哪些数学工具比较熟悉？ | Which mathematical tools are you comfortable with? | multi-choice | 微分方程 / ODEs, 线性代数 / Linear Algebra, 拉普拉斯变换 / Laplace Transform, 复数与相量 / Complex Phasors, 状态空间 / State-Space |

**Validation**: Name (Q00) is required before "Continue" is enabled. Q01-Q03 can be left unanswered.

### Section 2: Conceptual Baseline (8 questions)

All textarea questions show a subtle "还不了解 / I don't know yet" link below the input. Tapping it fills the answer with "还不了解" and visually marks the question as skipped (muted styling). The student can still clear it and type an answer.

| #  | Text (ZH) | Text (EN) | Type | Details |
|----|-----------|-----------|------|---------|
| 04 | 请简述交流感应电机和永磁同步电机的主要区别。 | Briefly describe the main difference between an AC induction motor and a PMSM. | textarea + skip | Tests fundamental machine knowledge |
| 05 | 电磁转矩是如何产生的？用你自己的话解释。 | How is electromagnetic torque produced? Explain in your own words. | textarea + skip | Open-ended, reveals depth |
| 06 | Clarke变换和Park变换的目的是什么？ | What is the purpose of the Clarke and Park transforms? | textarea + skip | Core FOC prerequisite |
| 07 | 在矢量控制（FOC）中，为什么要把 id 控制为零？（可多选） | In FOC, why do we typically regulate id to zero? (select all that apply) | **multi-choice** + skip | 减少铜损 / Reduce copper loss, 最大转矩电流比 / Max torque-per-ampere, 简化控制 / Simplify control, 不确定 / Not sure |
| 08 | SVPWM 和 SPWM 有什么区别？ | What is the difference between SVPWM and SPWM? | textarea + skip | Covered in ACMSimPy ep3 |
| 09 | PID 控制器的三个参数各自起什么作用？ | What role does each of the three PID parameters play? | textarea + skip | Fundamental controls knowledge |
| 10 | 什么是Bode图？它能告诉我们什么？ | What is a Bode plot and what does it tell you? | textarea + skip | Frequency-domain analysis basics |
| 11 | 什么是增益裕度和相位裕度？它们为什么重要？ | What are gain margin and phase margin? Why do they matter? | textarea + skip | Stability analysis basics |

**Validation**: All questions in this section are optional (can skip all). "Continue" is always enabled.

### Section 3: Expectations (2 questions)

| #  | Text (ZH) | Text (EN) | Type | Details |
|----|-----------|-----------|------|---------|
| 12 | 你最希望从这门课学到什么？ | What do you most want to learn from this course? | textarea | Helps tailor course content |
| 13 | 你更喜欢哪种学习方式？ | How do you prefer to learn? | multi-choice | 阅读理论 / Reading theory, 运行仿真 / Running simulations, 动手实验 / Hands-on experiments, 讨论与协作 / Discussion & collaboration |

**Validation**: Both optional. "Submit" button replaces "Continue" on this section.

## Submission Flow

1. Student fills in sections, progress bar advances per section.
2. On the final section (Section 3), the action button says "Submit" instead of "Continue".
3. Tapping "Submit" sends `POST /submit` with JSON payload:
   ```json
   {
     "name": "张三",
     "timestamp": "2026-03-06T14:30:22.000Z",
     "answers": {
       "q00_name": "张三",
       "q01_experience": "课程学过",
       "q02_languages": ["Python", "MATLAB"],
       "q03_math": ["微分方程", "线性代数"],
       "q04_im_vs_pmsm": "感应电机没有永磁体...",
       "q05_torque": "还不了解",
       "q06_clarke_park": "...",
       "q07_foc_id_zero": ["最大转矩电流比"],
       "q08_svpwm_spwm": "还不了解",
       "q09_pid": "P是比例...",
       "q10_bode": "还不了解",
       "q11_margins": "还不了解",
       "q12_expectations": "想学FOC的实现",
       "q13_learning_style": ["运行仿真", "动手实验"]
     }
   }
   ```
4. Server writes `responses/张三.json` (overwrites if exists) and responds `200 OK`.
5. Student sees a "Thank you / 感谢你的参与" confirmation screen with a checkmark animation.
6. `localStorage` is cleared after successful submission.

## Error Handling

- **Network failure on submit**: Show a retry button with the message "提交失败，请重试 / Submission failed, please retry". Do not clear the form.
- **Empty name**: "Continue" button on Section 1 stays disabled until name is non-empty.
- **Server not reachable**: If the initial page loads (HTML is cached), but POST fails, the retry flow handles it. If the page can't load at all, that's a network setup issue outside the app's scope.

## AI Analysis (Post-Collection)

After the lecture, the instructor feeds all `responses/*.json` to an AI (Claude) with a prompt like:

> Analyze these student survey responses for a graduate-level electric motor control course. Identify:
> 1. Overall class background distribution (prior experience, programming skills, math comfort)
> 2. Common misconceptions revealed in the conceptual answers (Q04-Q11)
> 3. Knowledge gaps that need extra attention in the course — especially around frequency-domain analysis (Bode plots, stability margins) and coordinate transforms (Clarke/Park)
> 4. Student expectations and interests (Q12) and preferred learning styles (Q13)
> 5. Suggested topics to emphasize or de-emphasize based on the class profile
> 6. How many students selected "还不了解" for each conceptual question — this shows which topics are completely unknown vs partially understood

This analysis shapes the teaching plan for the remaining lectures.

## Implementation Priority

1. `server.py` — minimal, working POST handler with overwrite-by-name
2. `lecture1.html` — mobile-first questionnaire with 3 section cards, localStorage save, bilingual
3. QR code convenience (optional: print URL on startup, or use `qrcode` package to generate a terminal QR code)
