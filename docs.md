# Interactive Circuit Canvas: Overview

This document provides an overview of the Interactive Circuit Canvas, a Python and Pygame-based tool for designing, analyzing, and understanding electronic circuits.

## What is it?

The Interactive Circuit Canvas is an application designed for:

*   **Circuit Design:**
    *   Users can drag and drop standard components: Resistors (R), Capacitors (C), Inductors (L), DC Voltage Sources (VDC), AC Voltage Sources (VAC), and Ground (GND).
*   **Circuit Wiring:**
    *   Connections are made by clicking from one component pin to another.
*   **Circuit Analysis:**
    *   Supports both DC (zero frequency) and AC (frequency domain) steady-state analysis.
*   **Solution Visualization:**
    *   Directly overlays calculated node voltages and component currents onto the schematic after solving.
    *   Supports displaying either the **phasor** (magnitude and phase) of AC values or an **animated instantaneous** value over time.
*   **Educational Insight:**
    *   Allows users to inspect the underlying Modified Nodal Analysis (MNA) matrices (Y, b).
    *   Features a viewer for the LU decomposition (P, L, U matrices) used in solving the system.
    *   Configuration options via `config.cfg` allow scaling the GUI and adjusting simulation parameters (like Gmin).

The tool is built using Python, Pygame for the graphical interface, NumPy for numerical operations, and SciPy for linear algebra functions (specifically LU decomposition).

## Configuration (`config.cfg`)

The application reads settings from a `config.cfg` file in the same directory. This file allows customization of GUI scaling and simulation parameters without modifying the source code.

```ini
[GUI]
# GUI Scaling Factors
# Adjust these to scale the overall window size and font sizes
# Default: 1.25 (25% larger than base)
size_scale = 1.25
# Default: 1.5 (50% larger font size)
font_scale = 1.5

# Base Window Dimensions (before scaling)
# These values will be multiplied by size_scale
# WIN_W = CANVAS_W + PAL_W
# WIN_H = CANVAS_H + INFO_H
base_win_w = 1220
base_win_h = 820
base_pal_w = 220
base_info_h = 120

[SIMULATION]
# Simulation Constants
# Gmin (Minimum conductance added to diagonals for stability)
# Default: 1e-12
gmin_default = 1e-12
# R_L_DC_Short (Resistance used for 0H inductors or inductors at DC)
# Default: 1e-6
r_l_dc_short = 1e-6
```

If `config.cfg` is missing or values cannot be read, the application will use the default values specified in the code.

## Core Concept: Modified Nodal Analysis (MNA) Explained

Modified Nodal Analysis (MNA) is a powerful algorithm used to systematically generate the set of equations that describe an electronic circuit's behavior. Instead of manually applying Kirchhoff's laws and Ohm's law in an ad-hoc manner, MNA provides a structured recipe.

**How MNA Works:**

1.  **Node Identification:**
    *   First, every distinct connection point in the circuit is identified as a "node."
    *   One of these nodes is chosen as the **reference node (ground)**, and its voltage is defined as 0 Volts. All other node voltages are measured relative to this ground.
    *   Let's say there are `N` non-ground nodes. This means we will have `N` unknown node voltages (e.g., `v_1, v_2, ..., v_N`).

2.  **Kirchhoff's Current Law (KCL) Equations:**
    *   For each of the `N` non-ground nodes, KCL is applied. KCL states that the sum of all currents entering a node must equal the sum of all currents leaving it (or, the algebraic sum of currents at a node is zero).
    *   The current through passive components (R, L, C) is expressed in terms of the node voltages at their terminals and their impedance/admittance. For example, the current through a resistor `R` connected between node `a` and node `b` is `(v_a - v_b) / R`.
    *   This step yields `N` linear equations involving the unknown node voltages.

3.  **Handling Voltage Sources (The "Modified" Part):**
    *   **Ideal Voltage Sources:** If a circuit contains ideal voltage sources, simply applying KCL isn't enough because the current *through* an ideal voltage source is initially unknown.
    *   MNA addresses this by:
        *   Introducing a **new unknown variable** for the current flowing through each voltage source. If there are `M` voltage sources, this adds `M` new unknowns.
        *   Adding a **new equation** for each voltage source. This equation directly relates the node voltages at the terminals of the voltage source to its specified value (e.g., `v_positive_terminal - v_negative_terminal = V_source_value`).
    *   This "modification" to basic nodal analysis is why it's called "Modified Nodal Analysis."

**The Resulting System:**
After applying these steps, MNA produces a system of `N + M` linear algebraic equations with `N + M` unknowns (the `N` node voltages and the `M` currents through voltage sources). This system can be written in matrix form:

**`Y * x = b`**

*   **`Y` (Admittance Matrix / MNA Matrix):**
    *   This is an `(N+M) x (N+M)` matrix.
    *   The top-left `N x N` submatrix primarily contains terms related to admittances (1/impedance) of passive components connected to the nodes.
    *   The remaining parts of the `Y` matrix incorporate the contributions and constraints from the voltage sources (often as +1 or -1 entries linking node voltage equations to voltage source current equations).
    *   A small conductance (`Gmin`) is added to the diagonal elements corresponding to non-ground nodes for numerical stability, especially in circuits with loops of voltage sources or similar structures that can lead to singular matrices.
*   **`x` (Vector of Unknowns):**
    *   An `(N+M) x 1` column vector.
    *   The first `N` elements are the unknown node voltages (`v_1` to `v_N`).
    *   The next `M` elements are the unknown currents flowing through the voltage sources (`i_vs1` to `i_vsM`).
*   **`b` (Vector of Knowns / Excitations):**
    *   An `(N+M) x 1` column vector.
    *   Elements corresponding to KCL equations are typically set by independent current sources connected to the nodes (though the current canvas only includes voltage sources currently).
    *   Elements corresponding to voltage source equations are set by the specified values of the independent voltage sources (voltage values for VDC/VAC sources).

By systematically constructing `Y` and `b` based on the circuit's topology and component values, MNA provides a standardized way to describe any linear circuit.

## Solving `Yx = b`: LU Decomposition Explained

Once the MNA has formulated the system `Yx = b`, the next step is to solve for the vector `x`, which contains all the unknown node voltages and voltage source currents. For a system of linear equations, LU Decomposition is a common and efficient method, especially when dealing with the dense matrices that often arise from MNA.

**What is LU Decomposition?**

LU Decomposition is a matrix factorization technique. The idea is to decompose the original matrix `Y` into the product of two (or three, with pivoting) simpler matrices:

*   **`L` (Lower Triangular Matrix):** A square matrix where all the entries *above* the main diagonal are zero. In the context of the `scipy.linalg.lu_factor` function, `L` is typically a *unit* lower triangular matrix, meaning its main diagonal elements are all 1.
*   **`U` (Upper Triangular Matrix):** A square matrix where all the entries *below* the main diagonal are zero.

**Pivoting for Numerical Stability:**

Directly decomposing `Y` into `LU` can be numerically unstable or even impossible if, for example, a diagonal element becomes zero during the factorization process. To avoid this, **pivoting** is used. Pivoting involves reordering the rows (and sometimes columns) of the `Y` matrix to ensure that the largest possible (in magnitude) element is used as the "pivot" element at each step of the factorization. This significantly improves the numerical stability and accuracy of the solution.

When pivoting is included, the factorization becomes:

**`P @ Y = L @ U`**

*   **`P` (Permutation Matrix):** This matrix represents the row interchanges performed during pivoting. It's essentially an identity matrix with its rows reordered. Multiplying `Y` by `P` (i.e., `P @ Y`) effectively applies these row swaps to `Y`.

**How LU Decomposition Helps Solve `Yx = b`:**

1.  **Factorization:** First, the `Y` matrix is decomposed into `P`, `L`, and `U`. This is the most computationally intensive step but only needs to be done once for a given `Y`. The `scipy.linalg.lu_factor(Y)` function performs this, returning `L` and `U` often packed into a single matrix, and `P` as a pivot index array.

2.  **Substitution:** With the factorization `P @ Y = L @ U`, the original system `Yx = b` can be rewritten.
    Since `Y = Pᵀ @ L @ U` (using `P⁻¹ = Pᵀ` for permutation matrices), we have:
    `Pᵀ @ L @ U @ x = b`
    To get the original system `Yx=b`, we pre-multiply `Y` by `P` to get `P @ Y = L @ U`. The equation becomes `P @ (Y @ x) = P @ b`. Since `P @ Y = L @ U`, we substitute to get `L @ U @ x = P @ b`.

    This system is solved in two stages using the triangular nature of `L` and `U`:
    *   **Forward Substitution:** Let `z = U @ x`. Then the system becomes `Lz = P @ b`. Since `L` is lower triangular, `z` can be solved for easily and efficiently using forward substitution (solving for `z_1`, then `z_2` using `z_1`, and so on).
    *   **Backward Substitution:** Once `z` is known, we solve `Ux = z`. Since `U` is upper triangular, `x` can be solved for easily and efficiently using backward substitution (solving for `x_last`, then `x_last-1` using `x_last`, and so on).

The `scipy.linalg.lu_solve((lu, piv), b)` function takes the packed LU factors and pivot information from `lu_factor` and the `b` vector, and performs these forward and backward substitution steps to find `x`.

**Why LU Decomposition?**

*   **Efficiency for Multiple `b` Vectors:** If you need to solve `Yx = b` for the same `Y` but different `b` vectors (e.g., changing source values without changing circuit topology), the expensive LU factorization is done only once. Subsequent solves with new `b` vectors only require the much faster forward and backward substitution steps.
*   **Numerical Stability:** With pivoting, LU decomposition is a robust method for solving a wide range of linear systems.
*   **Determinant Calculation:** The determinant of `Y` can be easily calculated from the diagonal elements of `U` and the sign of the permutation `P`. (det(Y) = det(P⁻¹) * det(L) * det(U)). This can be useful for checking if the matrix is singular.

In the Interactive Circuit Canvas, LU decomposition provides a reliable engine for solving the MNA equations.

## User Interface (UI) and User Experience (UX) via Shortcuts & Interactions

The application is designed for interactive use, primarily through mouse actions and keyboard shortcuts. Understanding these interactions is key to using the canvas effectively.

**General Mouse Interactions:**

*   **Drag from Palette (Left-Click & Drag):**
    *   **UI:** The right-hand side of the screen contains a palette of circuit components (R, C, L, VDC, VAC, GND).
    *   **UX:** To add a component to the circuit, left-click on its icon in the palette, drag it onto the main canvas area (the gridded section), and release the mouse button. A "ghost" image of the component follows the mouse, snapping to the grid, to indicate where it will be placed.
*   **Place Component (Release Drag on Canvas):**
    *   **UX:** When the dragged component is released on the canvas, it's officially added to the circuit. The component will snap its center to the nearest grid intersection. If pins of the new component are close enough to existing pins of other components, wires may be automatically created within a small radius.
*   **Select Component (Left-Click on Component on Canvas):**
    *   **UI:** Clicking a component on the canvas selects it, usually indicated by a highlight or border.
    *   **UX:** Selected components can then be moved or deleted. By default, clicking a new component deselects any previously selected ones.
*   **Multi-Select Component (Ctrl + Left-Click on Component):**
    *   **UX:** Holding the `Ctrl` key while clicking components allows multiple components to be selected or deselected individually without affecting the selection state of other components. This is useful for moving or deleting a group.
*   **Move Component(s) (Left-Click & Drag Selected Component):**
    *   **UX:** Once a component (or multiple components) is selected, left-click and drag any of the selected components. All selected components will move together with the mouse, snapping to the grid. Wires connected to the moved components remain attached to their respective pins.
*   **Wire Mode (Left-Click on Pin, then Left-Click on another Pin):**
    *   **UI:** Component pins are small circles on the component bodies. A temporary line follows the mouse after the first pin is clicked.
    *   **UX:** To create a wire, left-click on a pin of one component. Then move the cursor to a pin on a *different* component and left-click again. A permanent wire is created connecting these two pins. Clicking the same pin again or a pin on the same component cancels wire mode.
*   **Delete Component (Right-Click on Component):**
    *   **UX:** Right-clicking on a component on the canvas immediately deletes it from the circuit, along with any wires connected to it. There is no confirmation, so this action is immediate.
*   **Adjust Component Value (Mouse Wheel over Component):**
    *   **UX:** Hovering the mouse cursor over a component (R, C, L, VDC, VAC) and scrolling the mouse wheel will change its value. Scrolling up multiplies the value by a factor (default 10), and scrolling down divides by that factor.
*   **Fine Adjust Component Value (Shift + Mouse Wheel over Component):**
    *   **UX:** Holding the `Shift` key while scrolling the mouse wheel over a component allows for finer adjustments to its value (multiplying/dividing by a smaller factor, default sqrt(10)).

**Keyboard Shortcuts & Their UI/UX Implications:**

*   **`S` - Solve Circuit:**
    *   **UI:** The status bar updates with the solution status (e.g., "Solved (LU)", "Singular matrix"). If successful and on the "Solution Overlay" page (0), calculated voltages/currents appear on the schematic.
    *   **UX:** This is the primary action to analyze the current circuit configuration. It triggers the MNA formulation and LU solving process. Results are then available for display and visualization. Stops AC animation.
*   **`M` - Cycle Info Pages:**
    *   **UI:** The content in the info bar at the bottom cycles through different pages of information (Heatmap, Y\*x=b Text, Circuit Info). A small page indicator (e.g., "Pg 1/3") appears. Page 0 (Overlay) is the default when `M` is not active.
    *   **UX:** Allows the user to inspect detailed data from the last solution (if any) without cluttering the main canvas. Stops AC animation.
*   **`E` - Circuit Info Page:**
    *   **UI:** Switches directly to the "Circuit Info" page (Page 3) in the info bar.
    *   **UX:** Provides a direct shortcut to view component details, connections, and solver information. Stops AC animation.
*   **`T` / `L` / `H` / `?` - Toggle Legend/Help:**
    *   **UI:** A large pop-up overlay appears, displaying a list of mouse and keyboard controls.
    *   **UX:** Provides a quick reference for all available actions. Press any of these keys again or `Esc` to close the legend. Stops AC animation.
*   **`U` - Toggle LU Decomposition Viewer:**
    *   **UI:** A large pop-up overlay appears, displaying the P, Y, L, and U matrices from the last LU decomposition (if available). The main canvas is obscured.
    *   **UX:** Provides a detailed look at the matrices involved in solving `Yx=b`. Requires the circuit to have been solved at least once. Press `U` again or `Esc` to close the viewer. This is primarily an educational/debugging tool. Stops AC animation.
*   **`A` - Toggle AC Steady-State Animation:**
    *   **UI:** Toggles the display mode on the Solution Overlay (Page 0) between static phasors (Magnitude∠Phase) and animating instantaneous values (real number changing over time). The status bar indicates the animation state. Voltage/current text color changes (e.g., Cyan for animating voltage, Magenta for animating current).
    *   **UX:** Provides a visual representation of how the steady-state AC values change over the cycle. Only works if the circuit is solved and the frequency is > 0 Hz. Pressing `A` when not applicable (not solved, freq=0) gives a status message. Starts the animation time from 0.0s. Stops animation on most other interactions (solve, change page, edit, drag, prompt, etc.).
*   **`W` - Toggle Wire Delete Mode:**
    *   **UI:** Status bar indicates "Wire Delete ON/OFF". When ON, hovering over wires may highlight them for deletion.
    *   **UX:**
        *   When **ON**: Left-clicking on a wire deletes it.
        *   When **OFF** (default): Left-clicking on pins initiates wiring. Stops AC animation.
*   **`Del` / `Backspace` - Delete Selected Component(s):**
    *   **UX:** If one or more components are currently selected (highlighted), pressing `Delete` or `Backspace` removes them from the circuit and re-solves. Stops AC animation.
*   **`Ctrl + A` - Select All:**
    *   **UI:** All components on the canvas become selected (highlighted).
    *   **UX:** A quick way to select every component, useful for moving the entire circuit or deleting everything.
*   **`Esc` - Cancel / Close:**
    *   **UI:** Closes active pop-ups (Legend, LU Viewer). Closes Info Pages by returning to Page 0. Cancels active text input prompts (Edit Value, Set Frequency). Deselects wire start pin if in wiring mode. Stops AC animation.
    *   **UX:** A general-purpose "undo" or "back out" key for many operations. If no specific action is active, it often deselects all components.

**Info Bar Interactions:**

*   **Status Display:**
    *   **UI:** The bottom-most bar.
    *   **UX:** Continuously provides feedback: current action (e.g., "Dragging Resistor..."), results of operations (e.g., "Solved (LU)"), error messages (e.g., "Singular matrix"), or tooltips when hovering over palette items. Includes information about the simulation type (DC/AC), frequency, and animation state (On/Off).
*   **Input Prompts:**
    *   **UI:** When `C` (Edit Value) or `F` (Set Frequency) is used, the info bar transforms into a text input field showing a prompt.
    *   **UX:** The user types directly into this area. A blinking cursor appears. `Enter` submits the input, `Backspace` deletes characters, `Esc` cancels the input mode and restores the previous status. Inputting values supports SI prefixes and phase for AC sources (e.g., `10k`, `2.2u`, `5V@45`).

This combination of direct mouse manipulation and keyboard shortcuts aims to provide an efficient workflow for both circuit construction and analysis, catering to users who prefer graphical interaction and those who leverage keyboard power.

## Information & Visualization Pages (Accessed via M / E keys)

*   **Page 0: Solution Overlay (Default)**
    *   **Content:** Displays calculated node voltages and component currents directly on the circuit schematic.
    *   **Visualization:** For AC circuits (Freq > 0), this page displays either the **phasor** (Magnitude∠Phase) or the **instantaneous value** over time, toggled by the 'A' key. For DC circuits (Freq = 0), it shows the DC voltage/current values (real numbers).

*   **Page 1: MNA Matrix Heatmap**
    *   **Content:** A visual representation of the absolute magnitude of the elements in the MNA matrix (`Y`) and the right-hand side vector (`b`). Darker cells represent larger magnitude values.
    *   **Visualization:** Provides a quick way to see the structure and relative magnitudes of the MNA equations.

*   **Page 2: Y\*x=b Matrix Text View**
    *   **Content:** A formatted text dump of the `Y` matrix and `b` vector, showing the complex values for each element and the labels for the unknown variables (`x`).
    *   **Visualization:** Allows detailed inspection of the generated equations. Elements are displayed using SI prefixes and polar/rectangular complex formats for readability.

*   **Page 3: Circuit Info**
    *   **Content:** Text summary including the number of nodes, sources, matrix dimension, frequency, Gmin, R_L_DC_Short, a list of solved variable values (node voltages and source currents with magnitude/phase), and a list of each component's name, value, and connected nets.
    *   **Visualization:** Provides a summary of the circuit configuration and the calculated results in a list format.

## LU Decomposition Viewer (Shortcut: 'U')

*   **Content:** Displays the matrices resulting from the LU decomposition (`P`, `L`, `U`) of the MNA matrix (`Y`). The original `Y` matrix is also shown for comparison.
*   **Visualization:** Illustrates the `P @ Y = L @ U` factorization. This is useful for understanding how the linear system is solved internally and for diagnosing issues related to matrix conditioning (though the display is simplified for space).
*   **Access:** Activated by pressing the 'U' key (requires the circuit to be solved). Closed by pressing 'U' or 'Esc'.

## Installation

1.  Make sure you have Python 3.7+ and pip installed.
2.  Ensure Pygame, NumPy, and SciPy are installed:
    ```bash
    pip install pygame numpy scipy
    ```
3.  Save the provided Python code block as `app.py`.
4.  Save the provided Configuration block as `config.cfg` in the same directory.
5.  Save this Documentation block as `README.md` in the same directory.

## Usage

1.  Open your terminal or command prompt.
2.  Navigate to the directory where you saved `app.py` and `config.cfg`.
3.  Run the script:
    ```bash
    python app.py
    ```
The GUI window should open. Refer to the Legend (press `T`, `L`, `H`, or `?`) for controls.

## Dependencies

*   pygame
*   numpy
*   scipy

## Conclusion

The Interactive Circuit Canvas aims to be a user-friendly and educational tool for simulating basic electronic circuits. By combining a graphical interface with robust MNA-based solving (using LU decomposition) and providing multiple layers of visualization and data inspection (including the LU factors themselves and live AC animation), it facilitates both practical circuit analysis and a deeper understanding of the underlying mathematical principles.