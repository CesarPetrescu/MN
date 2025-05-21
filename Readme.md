# Interactive Circuit Canvas

> **A visual, LU‑powered sandbox for designing and solving linear circuits in real time**

![screenshot](docs/screenshot.png)

---

## Table of Contents

1. [Key Features](#key-features)
2. [Quick Start](#quick-start)
3. [Installation](#installation)
4. [Configuration](#configuration)
5. [Using the Canvas](#using-the-canvas)
6. [Keyboard & Mouse Cheat‑Sheet](#keyboard--mouse-cheat-sheet)
7. [Behind the Scenes – MNA & LU](#behind-the-scenes)
8. [Troubleshooting](#troubleshooting)
9. [Contributing](#contributing)
10. [License](#license)

---

## Key Features

* **Drag‑and‑drop schematic editor** – place R, C, L, DC & AC sources, and GND with live grid‑snap.
* **Instant Solve (`S`)** – one‑keystroke Modified Nodal Analysis with pivoted LU decomposition.
* **Inline measurement overlay** – voltages & currents rendered directly on the schematic (phasor or animated time‑domain view).
* **Matrix inspectors** – heat‑map, full text view, and a dedicated pop‑up for `P·Y = L·U`.
* **Fully configurable** GUI scaling, Gmin, DC‑short resistance, and window geometry via `config.cfg`.
* **Zero‑install runtime** – pure Python ≥ 3.8, only `pygame`, `numpy`, and `scipy`.

---

## Quick Start

```bash
# clone or download this repo
$ python -m venv venv && source venv/bin/activate   # optional
$ pip install -r requirements.txt                  # pygame, numpy, scipy
$ python app.py                                    # launch ✨
```

Press **`T`** (or **L / H / ?**) in‑app to view all controls.

---

## Installation

### Prerequisites

* Python 3.8 +
* pip or any PEP 517‑compatible installer

### From PyPI *(coming soon)*

```bash
pip install circuit‑canvas
```

### From source

```bash
git clone https://github.com/your‑org/circuit‑canvas.git
cd circuit‑canvas
pip install -r requirements.txt
python app.py
```

---

## Configuration

All tunables live in **`config.cfg`** (auto‑generated on first run). If the file is missing or a key is absent the built‑in defaults are used.

```ini
[GUI]
size_scale   = 1.25   # global multiplier for everything                                                               
font_scale   = 1.50   # text‑only multiplier
base_win_w   = 1220   # unscaled window width  ( canvas + palette )
base_win_h   =  820   # unscaled window height ( canvas + info bar )
base_pal_w   =  220   # palette width
base_info_h  =  120   # info‑bar height

[SIMULATION]
gmin_default  = 1e‑12 # conductance added to every node for numerical stability
r_l_dc_short  = 1e‑6  # substitute R for 0 H inductors / DC ind. shorts
```

Edit, save, and restart to apply.

---

## Using the Canvas

### 1 • Add components

Drag icons from the right‑hand palette onto the grid. They snap automatically; release to place. Pins that drop near other pins auto‑wire.

### 2 • Wire

Click **pin → pin**. Middle‑drag shows a preview; click again to confirm. Right‑click a component to delete it, or press **`W`** to enter *Wire‑Delete* mode.

### 3 • Edit values

* Mouse‑wheel over a part multiplies/divides by 10.
* **Shift + wheel** uses √10 for finer steps.
* Press **`C`** and click a component to type an exact value (`1k`, `47u`, `5V@45`).

### 4 • Solve & view results

Press **`S`**. Node voltages and branch currents appear immediately. For AC (set with **`F`**), toggle between phasor and animated views with **`A`**.

### 5 • Deep dive (optional)

* **`M`** cycles: Overlay → Heat‑map → `Y·x = b` text → Circuit info.
* **`U`** pops up the LU factorization viewer.
* **`T`** shows a full control legend.

---

## Keyboard & Mouse Cheat‑Sheet

| Action            | Keys / Mouse               | Notes                               |
| ----------------- | -------------------------- | ----------------------------------- |
| Place component   | *Drag from palette*        | R, C, L, VDC, VAC, GND              |
| Wire pins         | **L‑click pin → pin**      | Esc cancels current wire            |
| Move              | **Drag selected**          | Grid‑snaps                          |
| Solve             | **S**                      | Re‑runs MNA & LU                    |
| Toggle info pages | **M**                      | 0‒3                                 |
| AC animation      | **A**                      | Only if F > 0 Hz                    |
| Edit value        | **C** then click           | Wheel = coarse / Shift+wheel = fine |
| Set frequency     | **F**                      | `1k`, `60` etc. DC=0 Hz             |
| Delete part(s)    | **Right‑click** or **Del** | Ctrl‑A selects all                  |
| Wire‑delete mode  | **W**                      | Hover wire, L‑click                 |
| Legend            | **T / L / H / ?**          | Esc closes                          |
| LU viewer         | **U**                      | Esc closes                          |

---

## Behind the Scenes

The canvas builds a **Modified Nodal Analysis** system every time you press `S`:

```
Y · x = b        # N + M linear equations
```

* **Nodes → N** unknown voltages (ground = 0 V).
* **Voltage sources → M** extra current unknowns & KVL rows.

We solve using **pivoted LU decomposition** (`scipy.linalg.lu_factor`/`lu_solve`):

```
P · Y = L · U     →    L · (U · x) = P · b
```

Forward‑ then backward‑substitution gives `x`.  Gmin and an optional DC‑short resistor for inductors guarantee the matrix is nonsingular in pathological cases.

---

## Troubleshooting

| Symptom                            | Likely Cause                                                       | Fix                                                  |
| ---------------------------------- | ------------------------------------------------------------------ | ---------------------------------------------------- |
| **"No GND component"**             | Forgetting the reference node                                      | Drop a GND anywhere.                                 |
| **"Singular matrix"**              | Floating nodes, series voltage loops, multiple unconnected grounds | Ensure every pin can reach GND, check wiring.        |
| **Instantaneous AC values freeze** | `A` toggled off, or frequency = 0 Hz                               | Press `A` or set a non‑zero freq with `F`.           |
| **GUI too big/ small**             | High‑DPI display                                                   | Tweak `size_scale` and `font_scale` in `config.cfg`. |

---

## Contributing

Pull requests are welcome!  Please open an issue first to discuss major changes.  Dev setup:

```bash
pip install -r requirements-dev.txt  # black, flake8, mypy, etc.
pre‑commit install                   # run linters automatically
```

---

## License

This project is licensed under the **MIT License** – see [`LICENSE`](LICENSE) for details.
