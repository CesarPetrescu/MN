---
marp: true
theme: uncover # Or your preferred theme
paginate: true
header: 'Circuit Canvas'
footer: 'MNA & LU'
size: 16:9
---

<!-- _class: lead -->
# Interactive Circuit Canvas
Python & Pygame Tool
Circuit Design & Analysis

---

## What It Is
Interactive circuit application.
Design with components.
Wire connections easily.
Analyze DC/AC circuits.
Visualize solutions directly.

---

## Core: MNA Intro
Modified Nodal Analysis.
Models circuit behavior.
Identifies circuit nodes.
Uses reference ground (0V).
Finds unknown node voltages.

---

## MNA: KCL Step
Kirchhoff's Current Law.
Applied at each node.
Currents sum to zero.
Express current via voltage.
Forms N linear equations.

---

## MNA: Voltage Sources
Handles ideal voltage sources.
Adds new current unknowns.
Adds new voltage equations.
`V_pos - V_neg = V_src`.
This "modifies" nodal analysis.

---

## MNA: System `Yx = b`
MNA creates linear system.
`Y * x = b`.
`Y`: Admittance/MNA Matrix.
`x`: Unknown Voltages/Currents.
`b`: Known Source Values.

---

## Solving: LU Intro
`Yx = b` needs solving.
LU Decomposition is used.
Factors matrix `Y`.
Into `L` and `U`.
`L`: Lower triangular.
`U`: Upper triangular.

---

## LU: Pivoting
Pivoting ensures stability.
Reorders `Y` matrix rows.
Uses `P` (Permutation Matrix).
So, `P @ Y = L @ U`.
Improves solution accuracy.

---

## LU: Solving Steps
1. Factor: `P @ Y = L @ U`.
2. Let `z = U @ x`.
3. Solve `Lz = Pb` (Forward).
4. Solve `Ux = z` (Backward).
Efficient for `x`.

---

## LU Viewer ('T')
Popup shows LU matrices.
1. Original `Y` Matrix.
2. Permutation Matrix `P`.
3. Lower Factor `L`.
4. Upper Factor `U`.
Aids learning and debug.

---

## UI: Mouse Basics
Drag palette items: Place.
Click pin-pin: Wire.
Drag part: Move.
R-Click part: Delete.
Wheel on part: Value.

---

## UI: Key Shortcuts 1
`S`: Solve the circuit.
`M`: Info Pages (cycle).
`C`: Edit Values (toggle).
`F`: Set Frequency.
`W`: Delete Wires (toggle).

---

## UI: Key Shortcuts 2
`T`: LU Viewer (toggle).
`L`: Legend/Help (toggle).
`Del`: Delete selected.
`Ctrl+A`: Select all.
`Esc`: Cancel / Close.

---

<!-- _class: lead -->
## Conclusion
Canvas: Simulates circuits.
MNA forms equations.
LU solves reliably.
UI for interaction.
Visuals aid understanding.