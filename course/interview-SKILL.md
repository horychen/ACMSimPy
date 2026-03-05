---
name: interview
description: Interactive deep-dive interview to replace plan mode. Reads a spec/requirements file and conducts a thorough multi-round interview using AskUserQuestion, then writes the finalized spec. Use when the user wants to flesh out requirements, design decisions, or implementation details through structured Q&A before writing code.
disable-model-invocation: false
argument-hint: <spec-file-path>
---

# Interview Skill — Interactive Requirements Deep-Dive

You are conducting an **interactive interview** to deeply understand what the user wants before any code is written. This replaces plan mode with a more conversational, thorough approach.

## Input

The user provides a file path as `$ARGUMENTS`. This is typically a spec file, requirements doc, or design brief.

## Process

### Phase 1: Read & Analyze

1. **Read the provided file** at path `$ARGUMENTS` using the Read tool.
2. Silently analyze the content. Identify:
   - Stated requirements (what's explicit)
   - Ambiguities (what's unclear or could be interpreted multiple ways)
   - Gaps (what's missing but needed for implementation)
   - Assumptions (what the spec assumes without stating)
   - Trade-offs (where choices have consequences)
   - Edge cases (what could go wrong or be overlooked)
   - Technical depth (what needs more engineering detail)
   - UX/UI concerns (user experience implications)
   - Architecture decisions (structural choices that are hard to change later)
   - Integration points (how this connects to existing systems)

### Phase 2: Interview

Conduct a **multi-round deep-dive interview** using the `AskUserQuestion` tool. Follow these rules:

- **Ask 2-4 questions per round** (use the multi-question capability of AskUserQuestion)
- **Continue for as many rounds as needed** until ALL ambiguities are resolved — do NOT rush or cut short
- **Never ask obvious questions** whose answers are clearly stated in the spec
- **Be specific, not generic** — reference exact parts of the spec, exact variable names, exact UI elements
- **Go deep** — ask follow-up questions based on previous answers; drill into surprising or unexpected answers
- **Cover all angles**: technical implementation, UI/UX, edge cases, error handling, performance, security, accessibility, testing strategy, deployment, maintenance
- **Challenge assumptions** — if something seems like it could go wrong, ask about it
- **Propose trade-offs** — present options with pros/cons and let the user choose
- **Use previews** (the `markdown` field on options) when comparing concrete alternatives like code approaches, UI layouts, or architecture patterns
- **Respect the user's expertise** — ask questions that make them think, not questions that waste their time

#### Question Categories (cycle through these across rounds):

1. **Clarification** — "The spec says X, but does that mean A or B?"
2. **Edge Cases** — "What should happen when...?"
3. **Trade-offs** — "We could do X (fast but rigid) or Y (flexible but complex). Which matters more?"
4. **Technical Depth** — "For the Z component, should we use approach A or B? Here's what each implies..."
5. **UX/UI** — "How should the user experience this? What happens on error/loading/empty states?"
6. **Integration** — "How does this interact with existing system X?"
7. **Constraints** — "Are there performance/size/compatibility requirements not mentioned?"
8. **Priority** — "If we had to cut scope, which features are must-have vs nice-to-have?"
9. **Future-proofing** — "Is this likely to change? Should we design for extensibility here?"
10. **Validation** — "Let me confirm my understanding: [summary]. Is this correct?"

#### Interview Flow:

- **Round 1**: High-level architecture & major ambiguities
- **Round 2-3**: Drill into specific features and components
- **Round 4-5**: Edge cases, error handling, UX details
- **Round 6+**: Validate understanding, tie up loose ends
- **Final Round**: Present a concise summary of all decisions and confirm

### Phase 3: Write the Spec

After the interview is complete:

1. **Read the original file again** to get the latest content
2. **Rewrite or update the spec file** at the same path (`$ARGUMENTS`) incorporating ALL interview decisions
3. The updated spec should:
   - Preserve the original structure where possible
   - Add new sections for newly discovered requirements
   - Mark resolved ambiguities with clear decisions
   - Include technical implementation notes where relevant
   - Add edge case handling specifications
   - Note any trade-off decisions and their rationale
   - Be detailed enough that a developer can implement without further questions

## Important Rules

- **NEVER skip the interview** — even if the spec seems complete, there are always hidden assumptions
- **NEVER write code during the interview** — this is purely about understanding requirements
- **NEVER make decisions for the user** — always ask, even if you think you know the answer
- **DO use the user's language** (Chinese if they write in Chinese, English if English)
- **DO keep a mental running tally** of all decisions made so far
- **DO tell the user how many more rounds you estimate** at the start and as you go
- **DO end the interview gracefully** — summarize all decisions before writing

## Output

The final deliverable is the **updated spec file** at `$ARGUMENTS`, enriched with all interview findings.
