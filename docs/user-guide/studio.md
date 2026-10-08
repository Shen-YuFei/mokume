# Mokume Studio

Mokume Studio is an optional local web workbench for the Rust-backed `mokume`
wheel. It adds a project browser, typed workflow forms, run history, logs,
artifacts, and an optional AI assistant. The native workflows remain fully
usable without configuring an AI provider.

## Install and start

```bash
pip install "mokume[studio]"
mokume studio
```

Studio listens only on `127.0.0.1`. By default it selects the first available
port starting at 8765 and opens the browser. Use an exact port or suppress the
browser launch when needed:

```bash
mokume studio --port 9000
mokume studio --no-browser
```

The command does not accept a project path. Choose **File > Open Folder** after
the page opens. The launch URL contains a one-time token; after it is exchanged,
the token is removed from the address bar and an HTTP-only local session cookie
is used instead.

## Workbench flow

1. Open a local project folder.
2. Select a native Mokume workflow, or inspect a protein matrix and SDRF.
3. Enter and validate the parameters.
4. Run directly, or review and approve the final parameter set proposed in
   Assistant **Agent** mode.
5. Follow progress in **Run History**, inspect logs, and open registered
   artifacts.

Studio reads project inputs without modifying them. Each run gets a new output
directory inside the selected project, and paths that escape the project root
are rejected. The active project cannot be changed while a run is in progress.
Cancelling a run terminates its worker rather than merely hiding the task.

For dataset inspection, the protein matrix must have one unique protein-ID
column followed by numeric sample-intensity columns. Supply the matching SDRF,
two distinct contrast labels, and the input scale explicitly. Vendor matrices
with extra annotation columns should first be converted to this analysis-matrix
shape.

## Menus

The menu bar keeps only actions that correspond to the proteomics workbench:

- **File** — Open Folder, Close Project, Exit Mokume Studio
- **Analysis** — Inspect Dataset, Validate Parameters, Run Analysis, Cancel Run,
  Run History
- **View** — Sidebar, Assistant, Bottom Panel, Artifacts, Full Screen
- **Help** — Documentation, Keyboard Shortcuts, System Status, About Mokume

Unavailable actions are disabled according to the current project and run
state. Menus support arrow-key navigation, Home/End, Escape, and the shortcuts
shown beside each action. On narrow screens, the file browser and Assistant
open as mutually exclusive drawers.

## Optional Assistant

Open **Provider** in the Assistant panel to configure OpenAI, Anthropic, or an
explicit OpenAI-compatible endpoint. By default, the API key is held only in
the current Studio server process and is not written to the project or Studio
database. Select **Persist Studio configuration** before saving to write the
configuration and API key to `mokume-studio-providers.json` at the Mokume Git
root. When running from an installed wheel without a source checkout, Studio
uses the per-user Mokume configuration directory selected by `platformdirs`
instead (for example, `~/.config/mokume/` on Linux). Source checkouts add the
file to `.gitignore`, and the configuration file is created with user-only
permissions. The file contains the API key in plain text, so do not force-add
it to version control. The authenticated provider configuration endpoint returns
the key to the local browser so reopening the form can display it; responses are
marked `no-store`. Unselecting persistence and saving removes the file while
keeping the configuration for the current Studio session.

The two modes have different authority:

- **Ask** reads the persisted inspection summary and explains the dataset or
  available methods. It has no compute tool.
- **Agent** is the only assistant mode with write authority. It can prepare only
  policy-bounded candidates and pauses at a final card containing the exact
  parameters and input snapshot; computation starts only after explicit
  approval.

The model receives the user's message and the derived inspection context, not
raw matrix rows. Scientific inspection, policy selection, and evaluation still
run through Mokume's deterministic `RecommendationService`. Without ground
truth, comparisons remain `exploratory_unranked` and do not declare a winner.

## Recovery and provenance

Run state, approval decisions, logs, and artifact metadata are stored in the
local Studio state directory. Each completed run records its canonical
parameters, input identities, software version, and output checksums. If Studio
restarts during a compute run or dataset inspection, the unfinished record is
marked interrupted or failed instead of being left permanently active.
