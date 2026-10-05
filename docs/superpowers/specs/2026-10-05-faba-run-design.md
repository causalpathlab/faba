# `faba run`: set up and run the pipeline in the terminal

> **Superseded on merge with 0.16.2:** only Shift+Enter opens the preview and starts the run; `G` and `y` do nothing (as in qc). `c` copies the command on every screen and `p` prints it and leaves.

Status: design, approved in conversation; not implemented.
Replaces: `faba all` (and its aliases `pipeline`, `full`, `magic`).
Builds on: the `faba qc` view (#12) for the shared pop-up, list and Shift+Enter helpers.

## Goal

Someone who does not know faba's flags can go from BAMs to a faba output directory without
reading `--help`. The command that ran is saved as a script that reproduces it exactly. Scripts
and scheduled jobs keep a plain command line, the same way `faba qc` does.

One `faba run` session sets up **one** pipeline run: one set of BAMs into one output directory.

## Out of scope

- A queue of several runs.
- Choosing which modality filters the later steps (by gene, or by site overlap within a window).
  The count step stays the filter, by genes and cells, as today. This idea is recorded separately
  with its open questions.
- Launching `faba qc` from the run view.
- A crate shared with `senna run`.

## 1. Command line

- `faba all` and its aliases are removed. `faba run` takes every flag `all` took, under the same
  names, and runs the same pipeline (SNP → count → [depth] → ATOI → m6A → APA).
- `faba run --batch-process BAM... -g GFF -f GENOME -o OUT [flags]` runs straight through, as
  `faba all` did. Under `--batch-process` the BAMs, `-g`, `-f` and `-o` are required; a missing
  one is an error naming it. Clap itself treats them as optional so the view can start without
  them.
- Plain `faba run` opens the setup view. Anything given on the command line pre-fills it.
- Without a terminal on stdin and stdout, and without `--batch-process`, `faba run` stops with an
  error naming `--batch-process`.
- The whole-run record stays `pipeline_summary.json`, beside each step's `{step}.run.json`; only
  its job name changes from `all` to `run`. Nothing else about the outputs changes, so `faba qc`
  reads them as before.
- Run-wide flags get the clap help heading **Common**: threads, sparse backend, `--no-zip`, cell /
  gene / UMI barcode tags, UMI dedup, minimum base and mapping quality. Each step's flags get a
  heading named after the step. `faba run --help` is grouped the same way the view is. Verbose
  logging is a top-level `faba` flag and is passed through to the run.
- Docs move from `faba all` to `faba run`: README command table and examples, the methods
  write-up (§9 and every mention), `--help` texts.

## 2. The view

Four screens. Tab / Shift+Tab move between them, 1–4 jump to one, `q` quits (refused while a run
is going, with the reason in the footer).

### 2.1 Inputs

- The main area is a file browser filtered to `.bam`. Arrows move, Enter opens a folder, ← goes
  up, **Space toggles a BAM in or out of the run**.
- A selected BAM is tagged **fg** (foreground: signal, passed positionally) or **bg**
  (background: control, passed as `--control-bam`). New selections are fg; `b` flips the tag.
- A selected BAM shows the batch name its outputs will be prefixed with, as the pipeline
  derives it. A BAM with no index beside it is marked.
- Rows beside the browser: GFF, genome, known SNPs (optional), output directory, threads.
  Enter on a file row opens the browser as a pop-up to pick one file of that kind. Enter on the
  output row types it. Threads is the same field as on the Flags screen.
- The output directory must be new or empty (the rule `faba qc` uses). It is suggested as
  `faba_out` beside the first selected BAM, numbered when that exists.

### 2.2 Steps

- A checklist: SNP, count, depth, ATOI, m6A, APA. Space toggles. It sets the `--skip-*` flags.
- m6A shows "needs a bg BAM" and cannot be on until a BAM is tagged bg.
- Turning depth on asks for its resolution in kb (`--depth-resolution-kb`); off clears it.
- A line states that count's genes and cells restrict the later steps.

### 2.3 Flags

- Every other `run` flag, read from the clap definition. No flag is named in the form code.
- Groups follow the clap help headings, **Common** first, then the steps in pipeline order.
  A step's group is dimmed while the step is off.
- A help pane shows the selected flag's long help and its default.
- ↑/↓ (PgUp/PgDn, Home/End) move. Space toggles an on/off flag. ←/→ cycle a choice. Enter types
  a value. `/` filters by name. `a` shows the hidden (advanced) flags. `r` resets one, `R` all.
- The command is checked by clap as it is edited. A value clap rejects is marked on its row with
  clap's message condensed to one line.
- Flags owned by the Inputs and Steps screens are not repeated here.

### 2.4 Run

- Reachable once a run has started.
- The log, the last 2000 lines, with live progress bars; End follows the tail, arrows scroll.
- A status line: the running step, elapsed time, and at the end success or the exit code.
- `s` asks, then interrupts the run (as Ctrl-C would); a second `s` kills it.
- After the run ends, the log stays and `q` leaves.

### 2.5 Preview and confirm

- Shift+Enter, or `G` (for terminals that cannot report Shift+Enter), from any screen opens the
  preview pop-up:
  - the exact `faba run --batch-process …` command, one argument per line;
  - the problems that block the run, if any;
  - where the script will be saved.
- Problems that block the run: no fg BAM; GFF or genome missing; output not new or empty; a value
  clap rejects; m6A on with no bg BAM.
- In the preview: Shift+Enter or `y` saves the script and starts the run; `c` copies the command;
  Esc goes back. Plain Enter does nothing.

## 3. Running

- On confirm: create the output directory, write the script, start the run, switch to Run.
- The script is `{out}/faba_run.cmd.sh`, created only if it does not exist, mode 755:

  ```bash
  #!/usr/bin/env bash
  # Made by `faba run` (faba X.Y.Z). Run it again with: bash faba_run.cmd.sh
  set -euo pipefail
  cd "$(dirname "$0")"
  if [ -f pipeline_summary.json ]; then
    echo "pipeline_summary.json exists; move the outputs away to run again" >&2
    exit 1
  fi
  "${FABA:-faba}" run --batch-process \
    /abs/a.bam /abs/b.bam \
    --control-bam /abs/c.bam \
    -g /abs/genes.gff -f /abs/genome.fa -o . \
    --flag value
  ```

  Inputs are absolute paths. `-o .` writes beside the script. Only flags changed from their
  defaults are written. Arguments are shell-quoted.
- The run is a child process of this same binary (`std::env::current_exe()`), running the
  script's command with the output directory as its working directory.
- The child's stderr goes through a pseudo-terminal so faba's progress bars render; where a
  pseudo-terminal cannot be opened, a pipe. Output is split on `\r` and `\n` into log lines and
  progress updates, kept in a shared buffer that the view reads on each tick.
- A failed or stopped run leaves the script, so it can be fixed and run again.

## 4. Code layout

Home: `src/pipeline/`, beside the existing pipeline code.

- `src/pipeline/args.rs`: the BAMs, `-g`, `-f`, `-o` optional to clap; `--batch-process`; help
  headings.
- `src/pipeline/run.rs`: the batch path, unchanged apart from the checks above.
- `src/pipeline/tui.rs`: the view's entry point, app state, screen and key routing, on
  data-beans `Screen` / `run_screen`; polls the child from `tick()`.
- `src/pipeline/tui/form.rs`: the clap-generic flag form: fields from `clap::Arg` (kind, default,
  help, hidden, choices, heading), argv of changed fields, the clap check and the blamed row.
- `src/pipeline/tui/inputs.rs`: the BAM browser with selection and fg/bg tags; the single-file
  rows; the output rule.
- `src/pipeline/tui/steps.rs`: the checklist to and from `--skip-*` and `--depth-resolution-kb`.
- `src/pipeline/tui/script.rs`: quoting, the script text, the create-only write.
- `src/pipeline/tui/child.rs`: spawning through a pseudo-terminal or pipe, the line/progress
  splitter, interrupt then kill.
- `src/pipeline/tui/draw.rs`: the four screens and the preview.
- The pop-up, list-window and Shift+Enter helpers move from the qc view into one shared module
  used by both views.
- New dependency: `rustix` (unix only) for the pseudo-terminal and signals.

## 5. Testing

- Form: every `run` flag appears in its group; changed values give an argv that clap parses back
  to the same arguments (compared through serde); a rejected value is blamed on its row.
- Inputs and Steps: Space and `b` give the right positional and `--control-bam` lists; the steps
  give the right `--skip-*` flags; m6A is refused without a bg BAM; the output rule.
- Script: exact text; quoting of paths with spaces and quotes; never overwrites; run with bash in
  a temporary directory, the guard refuses when `pipeline_summary.json` exists.
- Child: a small shell command's lines, progress updates and exit code are captured; stop works.
- App: key routing per screen; the preview's problem list; every screen drawn at a large and a
  tiny terminal size without panicking.
- End to end: the existing pipeline tests call `faba run --batch-process` instead of `faba all`.
