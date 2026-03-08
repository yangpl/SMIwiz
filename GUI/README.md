# SMIwiz GTK Launcher

This is a GTK3 interface to run `SMIwiz` with editable `inputpar.txt`. Make sure you install GTK first.

## Build

```bash
cd GUI
make
```

## Run

```bash
cd GUI
./gui
```

## What It Does

- Select project directory, run directory, and `SMIwiz` binary.
- Load and edit `run_dir/inputpar.txt`.
- Edit common parameters in a form (`n1`, `n2`, `n3`, `d1`, `d2`, `d3`, `order`,`nt`, `dt`, etc.).
- Sync controls:
  - `Text -> Form` parses `inputpar.txt` into form fields.
  - `Form -> Text` writes form values back to `inputpar.txt` text.
- Select a mode (`0..10`) from the GUI.
- Save `inputpar.txt` (auto-updates `mode=` and applies form values).
- Run `mpirun -n <ranks> <bin> $(cat inputpar.txt)` in the selected run directory.
- View live stdout/stderr logs in the right-side terminal panel.
- Stop a running job with `SIGTERM`.
- Optional plot hook after successful run (`Auto plot on success`, default command `bash plot.sh`).
