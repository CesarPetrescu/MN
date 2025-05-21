
#!/usr/bin/env python
# =======================================================================
# app.py  –  Interactive Circuit Canvas  (2025‑04‑23  ·  v11.2 - Live AC View, Legend T Key, Fixes)
# =======================================================================
#
# Adds a pop-up overlay to view LU decomposition (Y, P, L, U) - press 'U' (changed from T)
# Includes Gmin for stability and refined DC/AC component models.
# GUI Scaling and Base Window Size are now read from config.cfg.
# Live AC voltage/current display added for steady-state values - press 'A'.
# Legend/Help shortcut is now 'T' (also L, H, ?)
#
# Shortcuts (Full legend: T key)
#   S=Solve, M=Info Pages, C=Edit Vals, F=Set Freq, W=Del Wires, U=LU View, A=Animate AC, T=Legend
#   Drag Palette=Place, Click Pin-Pin=Wire, R-Click=Delete Part
# -----------------------------------------------------------------------
from __future__ import annotations
import sys, cmath, math, itertools, re
import pygame as pg
import numpy  as np
import scipy.linalg
import configparser # Import configparser for configuration file reading
import os           # Import os for path checking
from collections import defaultdict, deque

# Make NumPy print options more readable for debugging
np.set_printoptions(precision=3, suppress=True, linewidth=200)

# --- Configuration Loading ---
CONFIG_FILE = "config.cfg"

# Default values (fallback if config file is missing or invalid)
DEFAULT_SIZE_SCALE = 1.25
DEFAULT_FONT_SCALE = 1.5
DEFAULT_BASE_WIN_W = 1220
DEFAULT_BASE_WIN_H = 820
DEFAULT_BASE_PAL_W = 220
DEFAULT_BASE_INFO_H = 120
DEFAULT_GMIN = 1e-12
DEFAULT_RL_SHORT = 1e-6

# Variables to hold the loaded config values (initialized with defaults)
SIZE_SCALE = DEFAULT_SIZE_SCALE
FONT_SCALE = DEFAULT_FONT_SCALE
GMIN_DEFAULT = DEFAULT_GMIN
R_L_DC_SHORT = DEFAULT_RL_SHORT

def load_config():
    """Reads configuration from config.cfg and updates global settings."""
    global SIZE_SCALE, FONT_SCALE, GMIN_DEFAULT, R_L_DC_SHORT # Declare globals we intend to modify
    config = configparser.ConfigParser()
    config_loaded = False

    if os.path.exists(CONFIG_FILE):
        try:
            config.read(CONFIG_FILE)

            # --- GUI Settings ---
            if 'GUI' in config:
                gui_section = config['GUI']
                # Use getfloat and getint with defaults for robustness
                SIZE_SCALE = gui_section.getfloat('size_scale', DEFAULT_SIZE_SCALE)
                FONT_SCALE = gui_section.getfloat('font_scale', DEFAULT_FONT_SCALE)
                # Note: Base dimensions are read but used *after* this function returns
                # when defining the final scaled window constants.
                # print(f"Config [GUI] loaded from {CONFIG_FILE}") # Optional: Verbose config loading
                config_loaded = True
            else:
                print(f"Warning: '{CONFIG_FILE}' found but no '[GUI]' section. Using default GUI settings.")

            # --- Simulation Settings ---
            if 'SIMULATION' in config:
                sim_section = config['SIMULATION']
                GMIN_DEFAULT = sim_section.getfloat('gmin_default', DEFAULT_GMIN)
                R_L_DC_SHORT = sim_section.getfloat('r_l_dc_short', DEFAULT_RL_SHORT)
                # print(f"Config [SIMULATION] loaded from {CONFIG_FILE}") # Optional: Verbose config loading
                config_loaded = True # Mark as loaded even if only SIMULATION section exists

            # Optional: Print a summary of loaded scales
            # print(f"Active Scaling: SIZE_SCALE={SIZE_SCALE:.2f}, FONT_SCALE={FONT_SCALE:.2f}")

        except (configparser.Error, ValueError) as e:
            print(f"Error reading or parsing '{CONFIG_FILE}': {e}. Using default settings.")

    # Optional: Print a warning if config wasn't loaded at all
    # if not config_loaded:
    #      print(f"Warning: '{CONFIG_FILE}' not found or could not be fully loaded. Using default settings.")
         # Globals are already initialized to defaults, no action needed here

# Load the configuration immediately when the script starts
load_config()


# ───────── UI CONSTANTS (Scaled, derived from loaded config) ─────────
# Re-read config here to get base dimensions *after* the global SIZE_SCALE is set
# This is a bit clunky, ideally base dimensions would be loaded into separate globals
_config_parser = configparser.ConfigParser()
_config_parser.read(CONFIG_FILE) # Read again to access the config object

_base_win_w = _config_parser['GUI'].getint('base_win_w', DEFAULT_BASE_WIN_W) if 'GUI' in _config_parser else DEFAULT_BASE_WIN_W
_base_win_h = _config_parser['GUI'].getint('base_win_h', DEFAULT_BASE_WIN_H) if 'GUI' in _config_parser else DEFAULT_BASE_WIN_H
_base_pal_w = _config_parser['GUI'].getint('base_pal_w', DEFAULT_BASE_PAL_W) if 'GUI' in _config_parser else DEFAULT_BASE_PAL_W
_base_info_h = _config_parser['GUI'].getint('base_info_h', DEFAULT_BASE_INFO_H) if 'GUI' in _config_parser else DEFAULT_BASE_INFO_H


# Calculate final scaled dimensions using the potentially updated SIZE_SCALE
WIN_W, WIN_H, PAL_W, INFO_H = int(_base_win_w*SIZE_SCALE), int(_base_win_h*SIZE_SCALE), int(_base_pal_w*SIZE_SCALE), int(_base_info_h*SIZE_SCALE)

CANVAS_W, CANVAS_H          = WIN_W-PAL_W, WIN_H-INFO_H
FPS, GRID, PIN_R            = 60, int(20*SIZE_SCALE), int(4*SIZE_SCALE) # FPS is typically not scaled
GHOST_A                     = 120 # Alpha is not scaled

FONT_PRIMARY_NAME = "DejaVu Sans Mono, Segoe UI, Arial, sans-serif"
FONT_MONO_NAME = "DejaVu Sans Mono, Consolas, Courier New, monospace"
FONT_SYMBOL_NAME = "DejaVu Sans Mono, Segoe UI Symbol, Arial Unicode MS, sans-serif"

# ───────── SIMULATION CONSTANTS (Now global, potentially loaded) ─────────
# These are now globals set by load_config: GMIN_DEFAULT, R_L_DC_SHORT


# ───────── PART TYPES ───────────
class CType: R, C, L, VDC, VAC, GND = range(6)
LBL = {CType.R:"R", CType.C:"C", CType.L:"L", CType.VDC:"DC", CType.VAC:"AC", CType.GND:"GND"}
COL = {CType.R:(215,215,215), CType.C:(120,170,255), CType.L:(130,255,170),
       CType.VDC:(120,200,255), CType.VAC:(255,170,120), CType.GND:(150,150,255)}

# ───────── BASIC CLASSES ─────────
class Pin:
    def __init__(self, comp:'Comp', dx:int, dy:int, name:str=""):
        # Pin offsets use the global SIZE_SCALE, which is set by load_config
        self.c, self.dx, self.dy = comp, int(dx*SIZE_SCALE), int(dy*SIZE_SCALE)
        self.net = -1; self.name = name
    @property
    def pos(self): return self.c.x+self.dx, self.c.y+self.dy

class Wire:
    def __init__(self, a:Pin, b:Pin): self.a=a; self.b=b

class Comp:
    comp_counters = {}
    # Component base dimensions are calculated *after* SIZE_SCALE is potentially loaded
    # These will be assigned after the class definition
    BW,BH = 0, 0 # Placeholder values

    def __init__(self, ct:CType, x:int, y:int):
        self.ct,self.x,self.y,self.sel = ct,x,y,False
        self.val  = {CType.R:1e3, CType.C:1e-6, CType.L:1e-3, CType.VDC:5.0, CType.VAC:5.0}.get(ct,0.0)
        self.phase= 0.0 # Phase in degrees for VAC
        Comp.comp_counters[ct] = Comp.comp_counters.get(ct, 0) + 1
        self.name = f"{LBL[ct]}{Comp.comp_counters[ct]}"
        # Pin offset uses the global SIZE_SCALE, which is set by load_config
        pin_offset = int(30*SIZE_SCALE)
        if ct==CType.GND: self.pins = [Pin(self,0,0, name="ref")]
        else: self.pins = [Pin(self,-pin_offset,0, name="1"), Pin(self,pin_offset,0, name="2")]

    def rect(self):
        # Rect calculations use the global SIZE_SCALE, which is set by load_config
        gnd_w, gnd_h = int(10*SIZE_SCALE), int(12*SIZE_SCALE)
        return pg.Rect(self.x-gnd_w, self.y-gnd_h, gnd_w*2, gnd_h*2) if self.ct==CType.GND \
               else pg.Rect(self.x-Comp.BW//2,self.y-Comp.BH//2,Comp.BW,Comp.BH)

    def label(self):
        # This method doesn't directly use SIZE_SCALE, but value_to_str might implicitly
        # format numbers derived from scaled component values. Keep as is.
        if self.ct == CType.R: unit = "\u03A9"
        elif self.ct == CType.C: unit = "F"
        elif self.ct == CType.L: unit = "H"
        else: unit = ""
        if unit: return f"{value_to_str(self.val)}{unit}"
        if self.ct==CType.VDC: return f"{self.val:.3g}V"
        # Use phase in degrees for VAC label
        if self.ct==CType.VAC: return f"{self.val:.3g}V\u2220{self.phase:.0f}\u00B0"
        return ""
    def __repr__(self): return self.name

# Assign the scaled class attributes for Comp *after* the class is defined
Comp.BW = int(80*SIZE_SCALE)
Comp.BH = int(40*SIZE_SCALE)


# ───────── VALUE PARSING / FORMAT ─────────
_SI = {"p":1e-12,"n":1e-9,"u":1e-6,"µ":1e-6, "m":1e-3,"":1,"k":1e3,"meg":1e6,"M":1e6,"g":1e9,"G":1e9}
def parse_value(txt:str) -> float|None:
    """Parses a string with optional SI prefix and returns float."""
    txt=txt.strip().lower()
    # Remove common unit suffixes for parsing number+prefix
    txt = re.sub(r"\s*(v|a|hz|f|h|ohm|\u03A9)\s*$", "", txt, flags=re.IGNORECASE)
    # Remove phase part if present (after value)
    if "@" in txt: txt=txt.split("@")[0].strip()
    if "∠" in txt or "\u2220" in txt : txt=txt.split("∠")[0].split("\u2220")[0].strip()

    # Match number (integer, float, scientific notation) and optional prefix
    m=re.match(r"([-+]?[0-9.]+(?:e[-+]?\d+)?)\s*(p|n|u|µ|m|k|meg|M|g|G)?",txt)
    if not m: return None
    try: return float(m.group(1))*_SI.get(m.group(2) or "",1)
    except ValueError: return None

def value_to_str(v:float) -> str:
    """Converts a float value to a string with SI prefixes (u, m, k, etc.)."""
    if v == 0: return "0"
    # Handle very small/large values or values very close to zero
    abs_v = abs(v)
    if abs_v < 1e-15: return "0"
    if abs_v > 1e15:
         return f"{v:.2e}" # Fallback to scientific notation for extreme values


    prefixes = [("G",1e9),("M",1e6),("k",1e3),("",1),("m",1e-3),("u",1e-6),("n",1e-9),("p",1e-12)]
    # Sort prefixes from largest to smallest factor
    prefixes.sort(key=lambda item: item[1], reverse=True)

    for suf, fac in prefixes:
        # Find the largest prefix where the scaled value is >= 1
        if abs_v >= fac or (abs_v < 1 and fac == 1): # Also include "" prefix for values >= 1
            sv = v / fac
            # Determine format based on magnitude of scaled value
            if abs(sv) >= 100: fmt = "%.0f"
            elif abs(sv) >= 10: fmt = "%.1f"
            elif abs(sv) >= 1:  fmt = "%.2f"
            else: # If abs_v < 1 but no smaller prefix matched (only "" remains)
                 fmt = "%.3f" # For values < 1 but >= m

            str_v = (fmt % sv)
            # Remove trailing zeros and potentially decimal point for whole numbers
            if '.' in str_v and 'e' not in str_v.lower():
                 str_v = str_v.rstrip('0').rstrip('.')

            return f"{str_v}{suf}"

    # Fallback just in case (shouldn't be reached with the prefixes list)
    return f"{v:.2e}"


# ───────── CIRCUIT / SOLVER ─────────
class Circuit:
    def __init__(self):
        self.comps, self.wires = [], []
        self.lu_raw = None; self.piv = None
        self.P_matrix = None; self.L_matrix = None; self.U_matrix = None
        self.invalidate()

    def invalidate(self):
        """Clears simulation results and netlist mapping."""
        self.solution=None; self.Y=self.b=None
        self.lu_raw = None; self.piv = None;
        self.P_matrix = None; self.L_matrix = None; self.U_matrix = None
        self.node_cnt=0; self.nmap={}; self.vsrc_idx={}; self.idx_to_nodename = {}
        self.net_id_to_comp_pins = {} # Maps net ID to list of "Comp.pin" strings

    def add_comp(self,c):
        self.comps.append(c)
        self.invalidate() # Circuit changes -> results invalid

    def add_wire(self,a:Pin,b:Pin):
        # Prevent self-loop wires on the same component
        if a.c is b.c: return
        # Prevent duplicate wires (check both directions)
        if any((w.a is a and w.b is b) or (w.a is b and w.b is a) for w in self.wires): return
        self.wires.append(Wire(a,b))
        self.invalidate() # Circuit changes -> results invalid

    def delete_comp(self,c_to_delete):
        if c_to_delete in self.comps:
            self.comps.remove(c_to_delete)
        # Remove any wires connected to the deleted component
        self.wires=[w for w in self.wires if w.a.c is not c_to_delete and w.b.c is not c_to_delete]
        self.invalidate() # Circuit changes -> results invalid

    def delete_wire(self,w_to_delete):
        if w_to_delete in self.wires:
            self.wires.remove(w_to_delete)
            self.invalidate() # Circuit changes -> results invalid


    def _nets(self)->tuple[bool, str | None]:
        """Identifies connected nets (nodes) via BFS."""
        # print("\n--- Starting Netlisting ---") # Optional: Verbose
        self.net_id_to_comp_pins.clear() # Clear previous net mapping
        for p in (p_ for c_ in self.comps for p_ in c_.pins):
            p.net = -1 # Reset all pin net IDs

        def bfs(start_pin:Pin, net_id:int):
            """Performs Breadth-First Search to find all pins in a net."""
            if start_pin.net != -1: return # Already visited

            queue=[start_pin]; start_pin.net=net_id
            # Store which pins belong to this net ID
            if net_id not in self.net_id_to_comp_pins:
                 self.net_id_to_comp_pins[net_id] = []
            self.net_id_to_comp_pins[net_id].append(f"{start_pin.c.name}.{start_pin.name}")

            head = 0
            while head < len(queue):
                current_pin = queue[head]
                head += 1 # Dequeue

                # Find connected pins via wires
                for w in self.wires:
                    neighbor = None
                    if w.a is current_pin and w.b.net == -1: neighbor = w.b
                    elif w.b is current_pin and w.a.net == -1: neighbor = w.a

                    if neighbor:
                        neighbor.net = net_id # Assign net ID
                        # Store pin-to-net mapping
                        if net_id not in self.net_id_to_comp_pins:
                            self.net_id_to_comp_pins[net_id] = []
                        self.net_id_to_comp_pins[net_id].append(f"{neighbor.c.name}.{neighbor.name}")
                        queue.append(neighbor) # Enqueue neighbor

        nid=0 # Current net ID to assign

        # Start BFS from the first GND component found (designating it Net 0)
        gnd_comps_found = [c for c in self.comps if c.ct == CType.GND]
        if gnd_comps_found:
             # Only BFS from the first GND if its pin hasn't been assigned yet
             if gnd_comps_found[0].pins[0].net == -1:
                 bfs(gnd_comps_found[0].pins[0], nid)
                 nid+=1
        # If multiple GNDs, they should have been connected already by the above BFS
        # If they are not connected, they will be assigned different net IDs below

        # Start BFS from any unvisited pin (including other GNDs if unconnected)
        for p in (p_ for c_ in self.comps for p_ in c_.pins):
            if p.net==-1:
                bfs(p,nid)
                nid+=1

        # Validate GND connections
        gnd_comps=[c for c in self.comps if c.ct==CType.GND]
        if not gnd_comps: return False,"No GND component. Circuit cannot be solved."

        # Ensure Net 0 is the ground net. If the first GND happened to get a non-zero
        # net ID because other floating nets were found first, swap its net ID with 0.
        # This is important for MNA matrix construction assuming Node 0 is ground.
        actual_ground_net_id = gnd_comps[0].pins[0].net
        if actual_ground_net_id != 0:
            # Swap net IDs 0 and actual_ground_net_id for all pins
            for p in (p_ for c_ in self.comps for p_ in c_.pins):
                if p.net == actual_ground_net_id: p.net = 0
                elif p.net == 0: p.net = actual_ground_net_id

            # Update net_id_to_comp_pins mapping accordingly
            # Need to handle potential edge case where original net 0 had pins
            pins_on_orig_net0 = self.net_id_to_comp_pins.pop(0, [])
            pins_on_orig_gnd_net = self.net_id_to_comp_pins.pop(actual_ground_net_id, [])

            if pins_on_orig_gnd_net: self.net_id_to_comp_pins[0] = pins_on_orig_gnd_net
            # Re-assign pins that were originally in net 0 to the swapped net ID
            if pins_on_orig_net0: self.net_id_to_comp_pins[actual_ground_net_id] = pins_on_orig_net0


        # Final check for multiple disconnected GNDs
        gnd_nets = set(c.pins[0].net for c in gnd_comps)
        if len(gnd_nets) > 1:
             # Find the net IDs that are not 0 (the intended ground)
             disconnected_gnd_nets = [net_id for net_id in gnd_nets if net_id != 0]
             has_disconnected_gnds_with_connections = False
             for c in gnd_comps:
                 if c.pins[0].net != 0:
                     # Check if the net this disconnected GND is on has other components/pins
                     # A net_id might not be in net_id_to_comp_pins if it only contains one (GND) pin and it was popped
                     if c.pins[0].net in self.net_id_to_comp_pins and len(self.net_id_to_comp_pins[c.pins[0].net]) > 1: # More than just the GND pin itself
                         has_disconnected_gnds_with_connections = True
                         break
                     elif c.pins[0].net not in self.net_id_to_comp_pins and any(p.net == c.pins[0].net for other_comp in self.comps if other_comp is not c for p in other_comp.pins):
                         # This GND's net ID is present on other components but wasn't in the map (unlikely but for safety)
                         has_disconnected_gnds_with_connections = True
                         break


             if len(gnd_nets) > 1 and has_disconnected_gnds_with_connections:
                 return False,"Multiple disconnected GNDs. Ensure all GND symbols are connected to the same net."
             elif len(gnd_nets) > 1:
                 # This case means there are multiple GND symbols, but any non-zero GND symbols
                 # are floating points with nothing else connected.
                 # Allow this for now, but it's potentially confusing. Could add a stricter warning.
                 pass # print("Warning: Multiple floating GND symbols detected. Only the one on Net 0 is active.")

        # print("--- Netlisting Complete ---") # Optional: Verbose
        return True, None


    def solve(self, freq:float):
        """Builds and solves the MNA matrix."""
        # print(f"\n=== Solving Circuit at f = {freq:.2f} Hz ===") # Optional: Verbose
        ok, msg = self._nets() # Perform netlisting first
        if not ok:
            self.invalidate()
            return False, msg

        omega = 2 * math.pi * freq

        # Map unique net IDs (excluding ground, net 0) to MNA matrix indices (0 to node_cnt-1)
        self.nmap.clear(); self.idx_to_nodename.clear(); k_node_idx=0
        # Get unique net IDs present in the circuit, excluding 0 (ground)
        unique_net_ids = sorted(list(set(p.net for c_ in self.comps for p in c_.pins if p.net > 0)))

        for net_id_val in unique_net_ids:
            self.nmap[net_id_val] = k_node_idx # Map net ID to matrix column/row index
            # Create a user-friendly label for this node index
            conn_comps_list = list(set(cp_str.split('.')[0] for cp_str in self.net_id_to_comp_pins.get(net_id_val, [])))
            conn_str = ", ".join(conn_comps_list[:2])
            if len(conn_comps_list) > 2: conn_str += ",..."
            mna_name = f"V(N{net_id_val}[{conn_str}])" if conn_str else f"V(N{net_id_val})"
            self.idx_to_nodename[k_node_idx] = mna_name # Store index -> name mapping
            k_node_idx += 1

        self.node_cnt = k_node_idx # Total number of non-ground nodes

        # Find all voltage sources (VDC, VAC) for MNA current equations
        v_sources = [c for c in self.comps if c.ct in (CType.VDC, CType.VAC)]
        m_vsrc_count = len(v_sources) # Number of voltage sources

        # Total matrix dimension: (non-ground nodes) + V_sources
        dim = self.node_cnt + m_vsrc_count

        # Handle empty circuit
        if dim == 0:
            # Reset all simulation results
            self.Y=np.array([],dtype=complex).reshape(0,0); self.b=np.array([],dtype=complex)
            self.solution=np.array([]);
            self.lu_raw=self.piv=self.P_matrix=self.L_matrix=self.U_matrix=None
            return True, "Solved (empty circuit)"

        # Initialize MNA matrix Y and right-hand side vector b
        Y = np.zeros((dim, dim), dtype=complex)
        b = np.zeros(dim, dtype=complex)

        # Helper function to get the MNA matrix index for a net ID
        # Ground (net 0) maps to index -1, which is ignored in MNA node equations
        def idx(net_id_val):
            # Check if net_id_val is a valid key in self.nmap for non-zero nets
            if net_id_val > 0:
                return self.nmap.get(net_id_val, -2) # Return -2 if > 0 but not mapped (shouldn't happen if netlisting is correct)
            return -1 if net_id_val == 0 else -2 # Return -1 for ground (net 0), -2 for any other unexpected negative ID


        # Add Gmin to diagonal for stability (applies to non-ground nodes, which are indices 0 to node_cnt-1)
        for i in range(self.node_cnt):
             Y[i,i] += GMIN_DEFAULT # Use loaded GMIN_DEFAULT

        # Stamp components into the Y matrix
        for c_comp in self.comps:
            admittance = 0j # Default admittance

            # Calculate admittance for passive components (R, C, L)
            if c_comp.ct==CType.R:
                # Avoid division by zero for 0-ohm resistors, treat as infinite conductance
                admittance = 1/max(c_comp.val,1e-12)
            elif c_comp.ct==CType.C:
                # Capacitor impedance is 1/(jωC). Admittance is jωC.
                if c_comp.val == 0 or omega == 0: # Open circuit at DC (omega=0) or for 0F cap
                    admittance = 0j
                else:
                    admittance = 1j*omega*c_comp.val
            elif c_comp.ct==CType.L:
                # Inductor impedance is jωL. Admittance is 1/(jωL).
                if c_comp.val == 0 or omega == 0:
                    # Short circuit at DC (omega=0) or for 0H inductor. Use small resistance -> large conductance.
                    admittance = 1/R_L_DC_SHORT # Use loaded R_L_DC_SHORT conductance for short
                else:
                    # Avoid division by zero if omega or val is zero, already handled above
                    admittance = 1/(1j*omega*c_comp.val)

            # Stamp passive components (R, C, L) into the nodal portion of Y
            if c_comp.ct in (CType.R,CType.C,CType.L):
                if len(c_comp.pins)!=2: continue # Skip if not a 2-terminal component

                n_a,n_b=c_comp.pins[0].net,c_comp.pins[1].net # Get net IDs of the two pins
                idx_a,idx_b=idx(n_a),idx(n_b) # Get corresponding MNA indices (-1 if ground)

                # Apply stamps for admittance between nodes a and b
                # Admittance adds to diagonal terms, subtracts from off-diagonal terms
                if idx_a>=0: # If node a is not ground
                     Y[idx_a,idx_a]+=admittance
                if idx_b>=0: # If node b is not ground
                     Y[idx_b,idx_b]+=admittance
                if idx_a>=0 and idx_b>=0: # If both nodes are not ground
                     Y[idx_a,idx_b]-=admittance
                     Y[idx_b,idx_a]-=admittance


        # Stamp voltage sources (VDC, VAC)
        # These add rows/columns related to KVL and branch currents
        self.vsrc_idx.clear() # Clear previous voltage source index mapping

        # Iterate through voltage sources and stamp them
        for s_counter,c_vs in enumerate(v_sources):
            if len(c_vs.pins)!=2: continue # Skip if not a 2-terminal source

            # The voltage source branch current is the MNA unknown at index node_cnt + s_counter
            mna_curr_eq_row=self.node_cnt+s_counter # Row index for the KVL equation (or current variable)
            self.vsrc_idx[c_vs]=s_counter # Store component -> index mapping for later lookup

            # Create a user-friendly label for this current variable
            self.idx_to_nodename[mna_curr_eq_row]=f"I({c_vs.name}:{c_vs.pins[0].name}\u2192{c_vs.pins[1].name})"

            # Get MNA indices for the positive and negative terminals of the source
            # Source pins are typically ordered [negative, positive]
            node_neg_idx,node_pos_idx=idx(c_vs.pins[0].net),idx(c_vs.pins[1].net)

            # KVL equations (rows in Y from node_cnt to dim-1)
            # The equation is V_pos - V_neg = V_source_value
            # Contributions to Y matrix:
            # If positive node is not ground, add +1 coefficient for V_pos column
            if node_pos_idx>=0:
                Y[mna_curr_eq_row,node_pos_idx]+=1
            # If negative node is not ground, add -1 coefficient for V_neg column
            if node_neg_idx>=0:
                Y[mna_curr_eq_row,node_neg_idx]-=1

            # Branch current contributions to KCL equations (columns in Y from node_cnt to dim-1)
            # The current flows from positive to negative terminal internally.
            # It leaves the positive node (-I_vs) and enters the negative node (+I_vs).
            # Contributions to Y matrix:
            # If positive node is not ground, add +1 coefficient for I_vs column in the positive node's KCL row
            if node_pos_idx>=0:
                Y[node_pos_idx,mna_curr_eq_row]+=1
            # If negative node is not ground, add -1 coefficient for I_vs column in the negative node's KCL row
            if node_neg_idx>=0:
                Y[node_neg_idx,mna_curr_eq_row]-=1

            # Set the source value in the right-hand side vector (b)
            phasor_b=0j
            if c_vs.ct==CType.VAC:
                # AC source value is magnitude * exp(j * phase_radians)
                phasor_b=cmath.rect(c_vs.val,math.radians(c_vs.phase))
            elif c_vs.ct==CType.VDC:
                 # DC source value. Only contributes if omega is 0 (DC analysis).
                 # For AC analysis (omega > 0), DC sources are treated as shorts (0V).
                 phasor_b=cmath.rect(c_vs.val,0) if omega==0 else 0j

            b[mna_curr_eq_row]=phasor_b # Assign the source value to the RHS vector


        self.Y,self.b=Y,b # Store the built matrices

        # Optional: Verbose matrix output
        # print("  --- MNA System (Y x = b) ---")
        # print("  Y matrix:\n",Y)
        # print("  b vector:\n",b)
        # print("  Unknowns order:", [self.idx_to_nodename.get(i, f"x{i}") for i in range(dim)])


        # Solve the system Y * x = b for x (the unknown voltages and currents)
        try:
            # Handle a likely empty or all-zero matrix if no sources/components are present
            if np.all(np.abs(Y)<1e-18) and np.all(np.abs(b)<1e-18) and dim>0:
                self.solution=np.zeros_like(b);
                 # Still need LU/pivots for matrix views, even if solution is trivial
                self.lu_raw,self.piv=(np.eye(dim),np.arange(dim)) if dim>0 else (None,None)
            else:
                # Use scipy's LU factorization and solve for efficiency and numerical stability
                # check_finite=False avoids checks that can be slow, use with caution or if sure inputs are not NaN/inf
                self.lu_raw,self.piv=scipy.linalg.lu_factor(Y,check_finite=False)
                self.solution=scipy.linalg.lu_solve((self.lu_raw,self.piv),b,check_finite=False)

            # Store P, L, U matrices from the factorization for the LU view
            if dim>0 and self.lu_raw is not None and self.piv is not None:
                self.L_matrix=np.tril(self.lu_raw,k=-1)+np.eye(dim); # L matrix (lower triangular, diagonal=1)
                self.U_matrix=np.triu(self.lu_raw)              # U matrix (upper triangular)
                # Construct P matrix from the pivot array. np.eye(dim)[self.piv] does this.
                self.P_matrix=np.eye(dim)[self.piv]

            # Optional: Verbose solution output
            # print("  --- Solution Vector x ---")
            # if self.solution is not None:
            #     for i_sol, val_sol in enumerate(self.solution):
            #         print(f"    x[{i_sol}] ({self.idx_to_nodename.get(i_sol, '?')}): {val_sol:.3g}")
            # print("===============================")

            return True,"Solved (LU)"

        except(ValueError,scipy.linalg.LinAlgError,np.linalg.LinAlgError) as e:
            # Handle solver errors (e.g., singular matrix for an improperly connected circuit)
            self.solution=None; self.lu_raw=self.piv=self.P_matrix=self.L_matrix=self.U_matrix=None
            # print(f"  !!! SOLVER ERROR: {type(e).__name__} - {e}") # Optional: Verbose
            # print("===============================")
            if "singular matrix" in str(e).lower() or "lu decomposition" in str(e).lower() or "exactly singular" in str(e).lower():
                 return False,"Singular matrix (check circuit connectivity/GND)"
            return False,f"Solver Error ({type(e).__name__})"


# ───────── APP ────────────
class App:
    def __init__(self):
        pg.init()
        # Font sizes use the global FONT_SCALE, which is set by load_config
        try: self.font = pg.font.SysFont(FONT_PRIMARY_NAME, int(15 * FONT_SCALE))
        except: self.font = pg.font.Font(None, int(15 * FONT_SCALE))
        try: self.font_small = pg.font.SysFont(FONT_MONO_NAME, int(12 * FONT_SCALE))
        except: self.font_small = pg.font.Font(None, int(12 * FONT_SCALE))
        try: self.font_tiny = pg.font.SysFont(FONT_MONO_NAME, int(10 * FONT_SCALE))
        except: self.font_tiny = pg.font.Font(None, int(10 * FONT_SCALE))
        try: self.font_symbol = pg.font.SysFont(FONT_SYMBOL_NAME, int(15 * FONT_SCALE))
        except: self.font_symbol = self.font
        # Screen size uses the global WIN_W, WIN_H, which are set by load_config
        self.scr=pg.display.set_mode((WIN_W,WIN_H)); pg.display.set_caption("Circuit Canvas v11.2 (Live AC View, Legend T Key, Fixes)") # Updated version title
        self.clock=pg.time.Clock(); self.circ=Circuit()
        self.drag_type:CType|None = None; self.drag_inst:Comp|None = None
        self.offx=self.offy=0
        self.wire_start:Pin|None = None
        self.mode_edit = False; self.mode_wire_del = False
        self.show_page = 0; self.show_help = False; self.show_lu_popup = False

        # --- New AC Animation State ---
        self.animate_ac = False
        self.sim_time = 0.0 # Time in seconds for animation

        # ───────── NEW "LIVE GRAPH" STATE ──────────
        self.show_graph_popup = False          # toggled with G
        self.graph_history_s  = 5.0            # seconds to keep
        self.graph_samples_hz = 60             # sampling rate
        self._next_sample_t   = 0.0            # internal timer
        #  history[(comp, 'V'|'I')]  → deque of (t, value)
        self.history = defaultdict(
                lambda: deque(maxlen=int(self.graph_history_s *
                                     self.graph_samples_hz)))

        # Updated initial status with new shortcuts
        self.status = "Drag parts to start. (T=Legend, U=LU View, A=Animate AC)"
        self.freq = 0.0
        self.input_buffer = ""; self.prompt = ""; self.active_edit:Comp|None = None
        # Palette rects calculation uses global CANVAS_W, PAL_W, SIZE_SCALE
        self.pal_rects={}; pal_item_x = CANVAS_W + int(24*SIZE_SCALE)
        pal_item_y_start, pal_item_y_step = int(28*SIZE_SCALE), int(58*SIZE_SCALE)
        pal_item_w, pal_item_h = PAL_W - int(48*SIZE_SCALE), int(44*SIZE_SCALE)
        self.pal_tooltips = {CType.R:"Resistor",CType.C:"Capacitor",CType.L:"Inductor",CType.VDC:"DC Source",CType.VAC:"AC Source",CType.GND:"Ground"}
        for i, t in enumerate((CType.R,CType.C,CType.L,CType.VDC,CType.VAC,CType.GND)):
            self.pal_rects[t] = pg.Rect(pal_item_x, pal_item_y_start+i*pal_item_y_step, pal_item_w, pal_item_h)
        Comp.comp_counters = {ct_val: 0 for ct_val in CType.__dict__.values() if isinstance(ct_val, int)}

    def get_comp_at(self,x,y):
        """Returns the component at the given screen coordinates, or None."""
        # Iterate components in reverse order to select topmost visually
        for c in reversed(self.circ.comps):
            # Comp rect uses global SIZE_SCALE
            if c.rect().collidepoint(x,y): return c
        return None

    def get_pin_at(self,x,y):
        """Returns the pin at the given screen coordinates within a radius, or None."""
        # Squared distance for click detection radius
        click_r_sq = (PIN_R + int(3*SIZE_SCALE))**2
        for c in self.circ.comps:
            for p in c.pins:
                # Check if point (x,y) is within the click radius of the pin's position
                if (p.pos[0]-x)**2+(p.pos[1]-y)**2 <= click_r_sq: return p
        return None

    def get_wire_at(self,x,y):
        """Returns the wire closest to the given screen coordinates within a tolerance, or None."""
        # Squared tolerance distance for click detection on a line
        tol_sq = int(6*SIZE_SCALE)**2
        for w in self.circ.wires:
            ax,ay=w.a.pos; bx,by=w.b.pos; # Pin positions
            dx,dy=bx-ax,by-ay; # Vector from pin A to pin B

            # Calculate the squared length of the wire segment
            len_sq = dx*dx + dy*dy

            if len_sq == 0: # Handle zero-length wires (shouldn't happen with pin-to-pin)
                t = 0 # Point is projected onto pin A
            else:
                # Project the mouse point (x,y) onto the line segment AB
                # t is the parameter along the line segment, clipped between 0 and 1
                t = max(0, min(1, ((x-ax)*dx + (y-ay)*dy)/len_sq))

            # Calculate the coordinates of the projected point on the line segment
            cx, cy = ax + t*dx, ay + t*dy

            # Check if the distance from the mouse point to the projected point is within tolerance
            if (x-cx)**2+(y-cy)**2 < tol_sq: return w
        return None

    def run(self):
        """Main application loop."""
        running = True
        while running:
            # Get time delta in seconds since last frame for animation speed independence
            dt = self.clock.tick(FPS) / 1000.0
            self.handle_events()

            # --- Update Simulation Time for Animation ---
            # Only update time if animation is active and frequency is meaningful
            if self.animate_ac and self.freq > 0:
                 self.sim_time += dt
                 self._maybe_sample_live_waveforms()
                 # Optional: Reset sim_time periodically to prevent floating point issues over very long runs
                 # A few periods are enough for visual effect
                 # if self.freq > 0:
                 #    period = 1.0 / self.freq
                 #    if self.sim_time > period * 5: # Reset after 5 periods
                 #         self.sim_time %= period # Keep time within a few cycles

            self.draw()

    def handle_events(self):
        """Handles Pygame events (keyboard, mouse)."""
        mx, my = pg.mouse.get_pos()
        ctrl_pressed = pg.key.get_mods() & pg.KMOD_CTRL
        shift_pressed = pg.key.get_mods() & pg.KMOD_SHIFT
        for e in pg.event.get():
            if e.type == pg.QUIT:
                pg.quit()
                sys.exit()

            # Handle events for the live graph popup
            if self.show_graph_popup:
                if e.type == pg.KEYDOWN and e.key in (pg.K_ESCAPE, pg.K_g):
                    self.show_graph_popup = False
                    self.status = "Graph closed."
                    continue

            # Handle events specific to the LU popup first
            if self.show_lu_popup:
                # LU popup uses U and Escape keys to close
                if e.type == pg.KEYDOWN and e.key in (pg.K_ESCAPE, pg.K_u): # Changed T to U
                    self.show_lu_popup = False
                    self.status="LU View closed."
                    continue # Consume event
                # Also close LU popup on any mouse click
                if e.type == pg.MOUSEBUTTONDOWN :
                    self.show_lu_popup = False
                    self.status="LU View closed."
                    continue # Consume event

            # Handle input prompt separately
            if self.prompt and e.type == pg.KEYDOWN:
                if e.key in (pg.K_RETURN, pg.K_KP_ENTER): self.process_prompt_input()
                elif e.key == pg.K_BACKSPACE: self.input_buffer = self.input_buffer[:-1]
                elif e.key == pg.K_ESCAPE: self.reset_action(); self.status="Input cancelled."
                elif e.unicode.isprintable(): self.input_buffer += e.unicode
                continue # Consume event

            # Handle general keyboard shortcuts
            if e.type == pg.KEYDOWN:
                if e.key == pg.K_ESCAPE:
                    # Close overlays/modes first
                    if self.show_lu_popup: self.show_lu_popup=False; self.status="LU View closed."
                    elif self.show_help: self.show_help=False; self.status="Legend closed."
                    elif self.show_page != 0 : self.show_page = 0; self.status="View closed."
                    # --- Stop animation on Escape ---
                    elif self.animate_ac: self.animate_ac = False; self.status = "Animation stopped."
                    else: self.reset_action(); self.status="Action cancelled."

                elif e.key in (pg.K_DELETE, pg.K_BACKSPACE):
                    if self.delete_selection():
                         self.status = "Deleted selected."
                         self.animate_ac = False # Stop animation if deleting

                elif e.key == pg.K_s:
                    ok, msg = self.circ.solve(self.freq)
                    self.status = msg + (" (T=Legend, U=LU View, A=Animate AC)" if ok else "") # Updated status message with U for LU, T for Legend
                    if ok: self.show_page = 0 # Switch to overlay view on successful solve
                    # --- Stop animation after solve ---
                    self.animate_ac = False
                    self.sim_time = 0.0 # Reset time after solve

                elif e.key == pg.K_m: # Cycle through info pages
                    if self.circ.Y is not None and self.circ.Y.shape[0] > 0:
                        self.show_page = (self.show_page + 1) % 4
                        pages = ["Overlay", "Heatmap", "Y*x=b", "Circuit Info"]
                        self.status = f"Info: {pages[self.show_page]}"
                        # --- Stop animation when changing info page ---
                        self.animate_ac = False
                    else: self.status = "Solve first (S) or add components."

                elif e.key == pg.K_e: # Direct to Circuit Info page
                    if self.circ.Y is not None and self.circ.Y.shape[0] > 0:
                        self.show_page = 3; self.status = "Info: Circuit Info"
                         # --- Stop animation when changing info page ---
                        self.animate_ac = False
                    else: self.status = "Solve first (S) or add components."

                elif e.key == pg.K_u: # Toggle LU View popup (Changed from T to U)
                    if self.circ.L_matrix is not None:
                        self.show_lu_popup = not self.show_lu_popup
                        self.status = "LU Decomposition View" if self.show_lu_popup else "LU View closed."
                        if self.show_lu_popup:
                            self.show_help = False; self.show_page=0 # Close other overlays
                            self.animate_ac = False # Stop animation
                    else: self.status = "Solve first (S) to view LU factors."

                elif e.key == pg.K_g:                     # NEW – Graph popup
                    self.show_graph_popup = not self.show_graph_popup
                    if self.show_graph_popup:
                        self.show_help = False; self.show_lu_popup = False
                        self.status = "Live graph view (G to close)"
                    else:
                        self.status = "Graph closed."

                elif e.key == pg.K_c: # Toggle Edit Values mode
                    self.mode_edit = not self.mode_edit
                    if self.mode_edit: self.mode_wire_del=False; self.reset_action(keep_edit=True)
                    self.status = "Component Edit "+("ON" if self.mode_edit else "OFF")
                    self.animate_ac = False # Stop animation

                elif e.key == pg.K_f: # Set Frequency
                    self.reset_action(keep_freq=True) # Keep frequency setting intention
                    self.prompt=f"Frequency (now={value_to_str(self.freq)}Hz): "; self.input_buffer=""
                    self.status = self.prompt
                    self.animate_ac = False # Stop animation

                elif e.key == pg.K_w: # Toggle Delete Wires mode
                    self.mode_wire_del = not self.mode_wire_del
                    if self.mode_wire_del: self.mode_edit=False; self.reset_action(keep_wire_del=True)
                    self.status = "Wire Delete "+("ON" if self.mode_wire_del else "OFF")
                    self.animate_ac = False # Stop animation

                elif e.key == pg.K_a and not ctrl_pressed: # --- Toggle AC Animation ---
                    if self.circ.solution is not None and self.circ.solution.size > 0 and self.freq > 0:
                        self.animate_ac = not self.animate_ac
                        if self.animate_ac:
                            self.sim_time = 0.0 # Reset time when starting animation
                            self.status = "AC Animation ON (Steady-State)"
                            self.show_page = 0 # Switch to overlay page
                            self.show_lu_popup = False # Close LU popup
                            self.show_help = False # Close legend
                        else:
                            self.status = "AC Animation OFF (Steady-State Phasors)"
                    elif self.circ.solution is None or self.circ.solution.size == 0:
                         self.status = "Solve circuit first (S) to animate AC."
                    else: # self.freq == 0
                         self.status = "Cannot animate DC circuit (Freq=0). Set AC frequency (F)."
                elif e.key == pg.K_a and ctrl_pressed: # Ctrl+A (Select All)
                    for c_ in self.circ.comps: c_.sel = True; self.status = "Selected all."

                elif e.key in (pg.K_l, pg.K_h, pg.K_QUESTION, pg.K_t): # Toggle Legend (Added T)
                    self.show_help = not self.show_help
                    self.status = "Legend " + ("shown" if self.show_help else "hidden")
                    # --- Stop animation when showing legend ---
                    if self.show_help: self.show_lu_popup = False; self.show_page=0; self.animate_ac=False

            # Handle mouse events (only if no popup or prompt is active)
            if e.type == pg.MOUSEBUTTONDOWN and not self.show_lu_popup and not self.prompt:
                # --- Stop animation on mouse interaction ---
                if self.animate_ac: self.animate_ac = False; self.status = "Animation stopped."

                if e.button == 1: # Left Click
                    if mx > CANVAS_W: # Clicking in palette area
                        self.reset_action(keep_freq=True) # Reset most actions, keep freq prompt possibility
                        for t_val, r_pal in self.pal_rects.items():
                            if r_pal.collidepoint(mx, my):
                                self.drag_type = t_val; self.status = f"Dragging {LBL[t_val]}..."
                                break
                    else: # Clicking in canvas area
                        if self.mode_wire_del:
                            w_del = self.get_wire_at(mx, my)
                            if w_del: self.circ.delete_wire(w_del); self.status="Wire deleted."; self.circ.solve(self.freq)
                        elif self.mode_edit:
                            c_edit = self.get_comp_at(mx, my)
                            if c_edit and c_edit.ct!=CType.GND: # Cannot edit GND value
                                self.active_edit = c_edit
                                self.prompt = f"New val for {c_edit!r} ({c_edit.label()}): "; self.input_buffer=""
                                self.status = self.prompt
                        else: # Normal click/drag
                            p_clk = self.get_pin_at(mx, my)
                            c_clk = self.get_comp_at(mx, my) if not p_clk else None # Prioritize pins

                            if p_clk: # Clicked a pin (for wiring)
                                self.drag_type = None # Cancel any palette drag
                                if self.wire_start and p_clk is not self.wire_start and p_clk.c is not self.wire_start.c:
                                    # Complete wire
                                    self.circ.add_wire(self.wire_start, p_clk); self.status = "Wire added."
                                    self.wire_start = None # Reset wire start
                                    self.circ.solve(self.freq) # Re-solve after wiring
                                elif self.wire_start is p_clk: # Clicked the same pin again
                                    self.wire_start=None; self.status="Wiring cancelled."
                                else: # Start a new wire
                                    self.wire_start = p_clk; self.status = "Click another pin..."

                            elif c_clk: # Clicked a component (for selection/dragging)
                                self.drag_type=None; self.wire_start=None # Cancel other modes
                                # Handle selection: Ctrl+Click to toggle, plain click to select only this one
                                if not ctrl_pressed and not c_clk.sel:
                                    for c_oth in self.circ.comps: c_oth.sel=False # Deselect others
                                c_clk.sel = not c_clk.sel if ctrl_pressed else True # Toggle or select
                                self.status = f"{c_clk!r} selected."
                                self.drag_inst = c_clk # Start drag for this component

                            else: # Clicking empty canvas
                                self.drag_type=None; self.wire_start=None
                                # Deselect all unless Ctrl is held
                                if not ctrl_pressed:
                                    for c_oth in self.circ.comps: c_oth.sel=False
                                self.status = "" # Clear status

                elif e.button == 3: # Right Click (for deleting components)
                    if not self.prompt: # Don't delete while prompting
                        c_del = self.get_comp_at(mx, my)
                        if c_del:
                            self.circ.delete_comp(c_del); self.status=f"Deleted {c_del!r}."
                            self.circ.solve(self.freq) # Re-solve after deleting
                            self.animate_ac = False # Stop animation if deleting

            # Handle mouse button release
            if e.type == pg.MOUSEBUTTONUP and not self.show_lu_popup and not self.prompt:
                if e.button == 1:
                    # If dragging from palette and released on canvas
                    if self.drag_type is not None and mx < CANVAS_W:
                        self.place(mx,my)
                    # Reset drag state
                    self.drag_type = None
                    self.drag_inst = None

            # Handle mouse motion (for dragging, wire preview, tooltips)
            if e.type == pg.MOUSEMOTION and not self.show_lu_popup and not self.prompt:
                # Dragging selected components
                if self.drag_inst is not None and e.buttons[0]: # Check if left button is held
                    mdx = mx + self.offx - self.drag_inst.x; mdy = my + self.offy - self.drag_inst.y
                    moved_any = False
                    for c_mv in self.circ.comps:
                        if c_mv.sel: # Move all selected components
                            new_x, new_y = c_mv.x + mdx, c_mv.y + mdy
                            # Component dimensions calculation uses global SIZE_SCALE
                            hbw, hbh = (Comp.BW//2, Comp.BH//2) if c_mv.ct != CType.GND else (int(10*SIZE_SCALE), int(12*SIZE_SCALE))
                            # Bounds checking uses global CANVAS_W, CANVAS_H
                            c_mv.x = max(hbw, min(CANVAS_W - hbw, new_x))
                            c_mv.y = max(hbh, min(CANVAS_H - hbh, new_y))
                            moved_any = True
                    if moved_any:
                         self.circ.invalidate() # Invalidate simulation if any components moved
                         self.animate_ac = False # Stop animation if moving
                # Hovering over palette items for tooltips
                elif self.drag_type is None: # Only show tooltips if not currently dragging a part type
                    current_tooltip = ""
                    if mx > CANVAS_W: # Mouse is in the palette area
                        for t_val, r_pal in self.pal_rects.items():
                            if r_pal.collidepoint(mx, my):
                                current_tooltip = self.pal_tooltips.get(t_val, LBL.get(t_val,""))
                                break
                    # Update status only if a tooltip is found or if the current status is a tooltip
                    # Avoid overwriting animation status unless a specific tooltip is active
                    if not self.animate_ac or current_tooltip: # Only update if not animating OR there's a specific tooltip
                        if current_tooltip and current_tooltip != self.status:
                            self.status = current_tooltip
                        # Clear tooltip if mouse moves away from the palette area and current status is a tooltip
                        elif not current_tooltip and mx <= CANVAS_W and self.status.startswith(tuple(self.pal_tooltips.values())):
                            self.status = "" # Clear tooltip if mouse leaves palette


            # Handle mouse wheel (for adjusting component values)
            if e.type == pg.MOUSEWHEEL and not self.show_lu_popup and not self.prompt:
                 if not self.prompt: # Don't adjust values while prompting
                    # --- Stop animation on mouse wheel value change ---
                    if self.animate_ac: self.animate_ac = False; self.status = "Animation stopped."
                    c_adj = self.get_comp_at(mx, my)
                    # Only adjust values for specific component types
                    if c_adj and c_adj.ct in (CType.R,CType.C,CType.L,CType.VDC,CType.VAC):
                        factor = 10.0 if e.y > 0 else 0.1 # Base factor 10x or 0.1x
                        if shift_pressed: factor = (10**0.5) if e.y > 0 else (1/(10**0.5)) # Shift for ~3.16x fine tune
                        c_adj.val *= factor
                        # Ensure R, C, L values don't go to zero or negative
                        if c_adj.ct in (CType.R,CType.C,CType.L): c_adj.val = max(c_adj.val, 1e-15) # Minimum non-zero value

                        self.status = f"Set {c_adj!r} to {c_adj.label()}";
                        self.circ.invalidate(); # Invalidate simulation
                        self.circ.solve(self.freq) # Re-solve immediately


    def process_prompt_input(self):
        """Processes input from the status bar prompt."""
        input_str = self.input_buffer.strip()
        if not input_str:
            self.reset_action(); self.status = "Edit cancelled."; return # Cancel if input is empty

        if self.prompt.startswith("New val"):
            c = self.active_edit
            if c is None: self.reset_action(); return # Should not happen if active_edit is managed correctly

            new_val = parse_value(input_str) # Parse the number+prefix
            new_phase = None
            if c.ct == CType.VAC:
                 # Look for phase part (@ or ∠) after the numerical value
                 match_phase = re.search(r"(?:@|∠|\u2220)\s*([-+]?\s*\d+\.?\d*)", self.input_buffer)
                 if match_phase:
                     try: new_phase = float(match_phase.group(1).replace(" ","")) # Parse phase in degrees
                     except ValueError: new_phase = None # Invalid phase input

            changed = False
            # Update value if valid number was parsed
            if new_val is not None:
                if c.ct in (CType.R,CType.C,CType.L):
                    if new_val > 1e-15: # Prevent zero/negative R, C, L
                        c.val = new_val; changed = True
                else: # VDC, VAC can be zero or negative
                    c.val = new_val; changed = True

            # Update phase if valid phase was parsed for a VAC source
            if new_phase is not None and c.ct == CType.VAC:
                 # Normalize phase to -180 to 180 degrees
                 new_phase = (new_phase + 180) % 360 - 180
                 c.phase = new_phase
                 changed = True

            if changed:
                 self.status=f"Updated {c!r} to {c.label()}";
                 self.circ.invalidate(); # Invalidate simulation
                 self.circ.solve(self.freq) # Re-solve with new value
            else:
                 self.status = f"Invalid input: '{self.input_buffer}'" # Report invalid input

        elif self.prompt.startswith("Frequency"):
            new_freq = parse_value(input_str) # Parse frequency value
            if new_freq is not None and new_freq >= 0:
                self.freq = new_freq; self.status = f"Freq set to {value_to_str(self.freq)}Hz."
                self.circ.invalidate(); # Invalidate simulation
                self.circ.solve(self.freq) # Re-solve with new frequency
            else:
                self.status = f"Invalid freq: '{self.input_buffer}'" # Report invalid input

        # Always clear prompt state after processing
        self.prompt = ""; self.input_buffer = ""; self.active_edit = None
        # --- Stop animation after prompt input ---
        self.animate_ac = False


    def reset_action(self, keep_freq=False, keep_edit=False, keep_wire_del=False):
        """Resets various transient UI states like dragging, wiring, prompting."""
        self.drag_type = None; self.drag_inst = None; self.wire_start = None
        if not keep_edit: self.mode_edit = False
        if not keep_wire_del: self.mode_wire_del = False

        # Clear prompt state unless specifically keeping an edit prompt
        if not (keep_edit and self.active_edit is not None):
             self.prompt = ""; self.input_buffer = ""; self.active_edit = None

        # --- Stop animation on most resets ---
        self.animate_ac = False
        self.sim_time = 0.0 # Reset animation time

        # --- Clear selection unless specific modes or overlays are active ---
        # Keep selection if in Edit mode, Wire Delete mode, or viewing an overlay/popup
        if not (keep_edit or keep_wire_del or self.show_help or self.show_page > 0 or self.show_lu_popup):
             for c_ in self.circ.comps: c_.sel = False # Deselect all

        # Set status text based on current state
        # Prioritize popups/overlays, then modes, then specific actions
        if self.show_lu_popup:
             self.status = "LU Decomposition View"
        elif self.show_help:
             self.status = "Legend shown"
        elif self.show_page > 0:
             pages = ["Overlay", "Heatmap", "Y*x=b", "Circuit Info"]
             # Adjusted page index/count for status bar display (excluding page 0)
             visible_pages = [1, 2, 3]
             page_index_in_visible = visible_pages.index(self.show_page) + 1 if self.show_page in visible_pages else self.show_page # Handle potential off-by-one if page 0 ever got here
             total_visible_pages = len(visible_pages)
             self.status = f"Info: {pages[self.show_page]} (Pg {page_index_in_visible}/{total_visible_pages})"
        elif self.prompt: # Prompt is active
            self.status = self.prompt # Status bar shows the prompt text
        elif keep_edit:
             self.status = "Component Edit ON"
        elif keep_wire_del:
             self.status = "Wire Delete ON"
        elif keep_freq: # Frequency prompt was intended but maybe not shown yet
            self.status = f"Set Frequency (now={value_to_str(self.freq)}Hz)"
        else: # No specific state is kept, return to default or cancelled
             mx, my = pg.mouse.get_pos()
             # Check if mouse is currently hovering a palette item for a tooltip status
             if mx > CANVAS_W:
                 current_tooltip = ""
                 for t_val, r_pal in self.pal_rects.items():
                     if r_pal.collidepoint(mx, my): current_tooltip = self.pal_tooltips.get(t_val, LBL.get(t_val,"")); break
                 if current_tooltip:
                      self.status = current_tooltip
                 else:
                      self.status = "Cancelled." # Default status if no tooltip
             else:
                 self.status = "Cancelled."


    def delete_selection(self):
        """Deletes all currently selected components and their connected wires."""
        sel = [c for c in self.circ.comps if c.sel]
        if not sel: return False # Nothing selected

        # --- Stop animation before deleting ---
        self.animate_ac = False
        self.sim_time = 0.0 # Reset animation time

        for c_del in sel:
            self.circ.delete_comp(c_del) # delete_comp also removes wires

        # Re-solve the circuit after deletion
        self.circ.solve(self.freq)
        return True # Indicate that something was deleted


    def place(self, x: int, y: int):
        """Places a component of the currently dragged type at the given canvas coordinates."""
        if self.drag_type is None:
            self.status = "Error: No type selected to place."; return # Should not happen if drag_type is managed correctly

        # --- Stop animation before placing ---
        self.animate_ac = False
        self.sim_time = 0.0 # Reset animation time

        ctype = self.drag_type
        # Snap to grid
        gx, gy = round(x/GRID)*GRID, round(y/GRID)*GRID
        # Prevent placing components completely outside the canvas
        hbw,hbh = (Comp.BW//2, Comp.BH//2) if ctype!=CType.GND else (int(10*SIZE_SCALE),int(12*SIZE_SCALE)) # Half width/height for bounds check
        gx = max(hbw, min(CANVAS_W-hbw, gx))
        gy = max(hbh, min(CANVAS_H-hbh, gy))

        new_c = Comp(ctype, gx, gy) # Create new component instance
        self.circ.add_comp(new_c) # Add to circuit

        self.status = f'Placed {new_c!r}'

        # Optional: Auto-wire to nearby pins (within a small radius)
        SNAP_SQ = (PIN_R * 3.5)**2 # Squared distance threshold for snapping
        snapped_wire = False
        # Find all pins *except* those on the new component
        other_pins = [p for c_ in self.circ.comps if c_ is not new_c for p in c_.pins]

        for p1 in new_c.pins:
            for p2 in other_pins:
                # Check distance between pins
                if (p1.pos[0]-p2.pos[0])**2 + (p1.pos[1]-p2.pos[1])**2 <= SNAP_SQ:
                    self.circ.add_wire(p1, p2); snapped_wire = True # Add a wire

        if snapped_wire:
            self.status += ' (auto-wired)' # Update status if wires were added

        # Re-solve the circuit after adding components/wires
        self.circ.solve(self.freq)

    # --- Drawing Methods ---
    # All draw_* methods below implicitly use the global constants
    # WIN_W, WIN_H, CANVAS_W, CANVAS_H, PAL_W, INFO_H, GRID, PIN_R, SIZE_SCALE, FONT_SCALE
    # which are now set based on the config file before App() is initialized.
    # The logic within them remains the same, relying on these globals.

    def draw(self):
        """Main drawing function."""
        self.scr.fill((40, 45, 50)) # Background color

        self.draw_grid()
        self.draw_wires()
        self.draw_comps() # This internally calls draw_solution_overlay if needed
        self.draw_palette()
        self.draw_info_bar()

        # Draw overlays/popups on top
        if self.show_graph_popup:
            self.draw_live_graph_popup()
        elif self.show_lu_popup:
            self.draw_lu_popup_overlay()
        elif self.show_help:
            self.draw_legend()
        else:
            self.draw_mini_legend() # Mini legend is shown by default

        pg.display.flip() # Update the full screen


    def draw_mini_legend(self):
        """Draws a small, always-visible legend in the corner."""
        if self.show_help or self.show_lu_popup: return # Don't show if main overlays are up

        # Updated mini legend to include 'A' and 'U' (T is now Legend)
        lines = ["S: Solve", "M: Info", "U: LU View", "A: Animate AC", "T: Legend", "Esc: Cancel"]
        line_h = self.font_tiny.get_linesize() # Use the smallest font

        padding = int(5 * SIZE_SCALE) # Spacing based on scale
        # Calculate total height needed for the text block
        total_h = len(lines) * line_h + max(0, len(lines)-1) * int(2*SIZE_SCALE) # Add extra spacing between lines

        # Position the text block in the bottom right corner
        start_y = WIN_H - padding - total_h
        start_x = WIN_W - padding # Right-aligned

        for i, text in enumerate(lines):
            surf = self.font_tiny.render(text, True, (170, 170, 190)) # Light gray text
            y_pos_line = start_y + i * (line_h + int(2*SIZE_SCALE)) # Position each line
            rect = surf.get_rect(right=start_x, top=y_pos_line) # Anchor to top-right
            self.scr.blit(surf, rect)


    def draw_grid(self):
        """Draws the background grid on the canvas area."""
        grid_c = (60,65,70) # Dark gray color
        # Draw vertical lines
        for x_ in range(0, CANVAS_W, GRID):
            pg.draw.line(self.scr, grid_c, (x_,0), (x_,CANVAS_H))
        # Draw horizontal lines
        for y_ in range(0, CANVAS_H, GRID):
            pg.draw.line(self.scr, grid_c, (0,y_), (CANVAS_W,y_))

    def draw_wires(self):
        """Draws all circuit wires."""
        wire_c = (220,220,130) # Yellowish color
        active_c = (255,255,0) # Bright yellow (for wire preview)
        del_c = (255,100,100) # Reddish (for delete mode hover)

        thick = max(1,int(2*SIZE_SCALE)) # Wire thickness
        active_thick = max(1,int(1*SIZE_SCALE)) # Preview wire thickness

        mx,my = pg.mouse.get_pos() # Current mouse position
        hover_w = self.get_wire_at(mx,my) if self.mode_wire_del else None # Check for hover only in delete mode

        for w in self.circ.wires:
            col = del_c if w is hover_w else wire_c # Use delete color if hovered
            # Draw line between the two pins of the wire
            pg.draw.line(self.scr, col, w.a.pos, w.b.pos, thick)

        # Draw preview wire if currently starting a wire
        if self.wire_start:
            pg.draw.line(self.scr, active_c, self.wire_start.pos, (mx,my), active_thick)

    def draw_comps(self):
        """Draws all components."""
        # Draw existing components
        for c in self.circ.comps:
            self.draw_single_comp(c)

        # Draw ghost preview of component being dragged from palette
        if self.drag_type is not None and pg.mouse.get_pos()[0] < CANVAS_W: # Check if dragging AND mouse is over canvas
            mx,my = pg.mouse.get_pos()
            # Snap ghost to grid
            gx,gy = round(mx/GRID)*GRID, round(my/GRID)*GRID

            # Get dimensions for bounds checking (GND is smaller)
            hbw,hbh = (Comp.BW//2, Comp.BH//2) if self.drag_type!=CType.GND else (int(10*SIZE_SCALE),int(12*SIZE_SCALE))

            # Keep ghost within canvas bounds
            gx = max(hbw, min(CANVAS_W-hbw, gx))
            gy = max(hbh, min(CANVAS_H-hbh, gy))

            ghost_color_rgb = COL[self.drag_type] # Get base color
            rect_w = Comp.BW if self.drag_type != CType.GND else hbw*2
            rect_h = Comp.BH if self.drag_type != CType.GND else hbh*2

            # Create a surface with alpha for transparency
            ghost_s = pg.Surface((rect_w, rect_h), pg.SRCALPHA)
            ghost_rect_on_surf = pg.Rect(0,0,rect_w, rect_h) # Rect relative to the surface

            # Draw the ghost rectangle with transparency
            pg.draw.rect(ghost_s, (*ghost_color_rgb, GHOST_A), ghost_rect_on_surf, border_radius=int(6*SIZE_SCALE))

            # Draw the component label on the ghost surface
            label_s = self.font.render(LBL[self.drag_type], True, (255,255,255, GHOST_A+80)) # Label text color (white with more alpha)
            ghost_s.blit(label_s, label_s.get_rect(center=ghost_rect_on_surf.center)) # Center label on ghost rect

            # Blit the ghost surface onto the main screen
            self.scr.blit(ghost_s, (gx - rect_w//2, gy - rect_h//2)) # Position centered at snapped grid location


    def draw_single_comp(self, c: Comp):
        """Draws a single component instance."""
        r = c.rect(); # Get the component's drawing rectangle
        base_col = COL[c.ct]; # Component type color
        border_c=(60,60,60); # Dark gray border color
        text_c=(0,0,0); # Black text color for name/value

        pin_c=(240,240,120); # Yellowish pin color
        sel_c=(255,255,0); # Bright yellow for selected components
        edit_c=(255,100,0); # Orange for component hovered in edit mode

        br = int(6*SIZE_SCALE); # Border radius based on scale
        thick = max(1,int(SIZE_SCALE)); # Border thickness
        sel_inf = int(6*SIZE_SCALE); # Inflation for selection/edit highlight border

        # Special drawing for Ground component
        if c.ct == CType.GND:
            x,y = c.x,c.y; # Center position
            stem_h=int(12*SIZE_SCALE); # Height of the vertical stem
            bar_y0=int(4*SIZE_SCALE); # Y offset for the top bar
            bar_ystep=int(4*SIZE_SCALE); # Spacing between horizontal bars
            bar_ws=[int(w_val*SIZE_SCALE) for w_val in (14,10,6)]; # Widths of the bars (top to bottom)
            gnd_thick = max(1, int(2*SIZE_SCALE)); # Thickness for ground symbol lines

            # Draw stem
            pg.draw.line(self.scr,base_col,(x,y-stem_h),(x,y),gnd_thick)
            # Draw horizontal bars
            for i,w_bar in enumerate(bar_ws):
                 pg.draw.line(self.scr,base_col,(x-w_bar//2,y+bar_y0+i*bar_ystep),(x+w_bar//2,y+bar_y0+i*bar_ystep),gnd_thick)
        else:
            # Draw rectangular components (R, C, L, VDC, VAC)
            pg.draw.rect(self.scr,base_col,r,border_radius=br) # Filled rectangle
            pg.draw.rect(self.scr,border_c,r,thick,border_radius=br) # Border

            # Draw component name (e.g., "R1", "C2")
            name_s = self.font.render(c.name, True, text_c)
            self.scr.blit(name_s, name_s.get_rect(center=(c.x, r.top + int(10*SIZE_SCALE)))) # Position above center

            # Draw component value/label (e.g., "1kΩ", "10uF")
            val_s = self.font_symbol.render(c.label(), True, text_c)
            self.scr.blit(val_s, val_s.get_rect(center=(c.x, r.bottom - int(10*SIZE_SCALE)))) # Position below center

        # Draw pins for all components
        for p_pin in c.pins:
            pg.draw.circle(self.scr,pin_c,p_pin.pos,PIN_R) # Filled circle for pin
            pg.draw.circle(self.scr,border_c,p_pin.pos,PIN_R,thick) # Border for pin

        # Draw selection highlight
        if c.sel:
             pg.draw.rect(self.scr,sel_c,r.inflate(sel_inf,sel_inf),max(1,int(2*SIZE_SCALE)),border_radius=br+int(2*SIZE_SCALE))

        # Draw edit mode hover highlight
        if self.mode_edit and c.ct!=CType.GND and r.collidepoint(pg.mouse.get_pos()):
             pg.draw.rect(self.scr,edit_c,r.inflate(sel_inf,sel_inf),max(1,int(2*SIZE_SCALE)),border_radius=br+int(2*SIZE_SCALE))

        # Draw simulation solution overlay (voltage/current values)
        # Only draw if solution exists and on the primary overlay page (page 0)
        if self.circ.solution is not None and self.show_page == 0 and self.circ.solution.size > 0 :
            self.draw_solution_overlay(c)


    def draw_solution_overlay(self, c: Comp):
        """Draws voltage and current values on components based on the solution."""
        if self.circ.solution is None or not self.circ.solution.size: return # No solution to display

        # Define colors for static phasor vs. animated instantaneous display
        volt_col_static = (120,255,120) # Green (Phasor Voltage)
        curr_col_static = (255,200,120) # Orange (Phasor Current)
        volt_col_animated = (180, 255, 255) # Cyan (Instantaneous Voltage)
        curr_col_animated = (255, 180, 255) # Magenta (Instantaneous Current)

        font_to_use=self.font_symbol # Font for values

        # Split solution vector into node voltages and voltage source currents
        node_voltages = self.circ.solution[:self.circ.node_cnt]
        vsrc_currents = self.circ.solution[self.circ.node_cnt:]
        net_map = self.circ.nmap # Mapping from net ID to node index

        omega = 2 * math.pi * self.freq # Angular frequency
        current_time = self.sim_time # Use elapsed simulation time for animation

        # --- Draw Voltage for Pins ---
        for p_pin in c.pins:
            volt_phasor = 0j # Default voltage is 0 relative to ground

            # Get voltage phasor for the pin's net (if it's not ground, net 0)
            if p_pin.net > 0 and p_pin.net in net_map:
                node_idx = net_map[p_pin.net] # Get the MNA node index for this net
                if 0 <= node_idx < len(node_voltages):
                    volt_phasor = node_voltages[node_idx] # Get the solved voltage phasor

            v_str = "" # String to display
            display_animated = self.animate_ac and self.freq > 0 # Are we animating AC?

            if display_animated:
                # Calculate instantaneous voltage: |V| * cos(ωt + φ)
                mag, ph_rad = cmath.polar(volt_phasor)
                instantaneous_v = mag * math.cos(omega * current_time + ph_rad)
                # Format as a real number string (e.g., "5.12V")
                v_str = f"{instantaneous_v:.3g}V"
                v_surf = font_to_use.render(v_str, True, volt_col_animated) # Use animated color
                # Optional: Vary color intensity based on instantaneous value (more complex)
                # intensity = min(255, max(0, int(abs(instantaneous_v / (mag + 1e-9)) * 255)))
                # v_surf = font_to_use.render(v_str, True, (180, 180+intensity//2, 180+intensity)) # Example color variation
            else:
                # Display static phasor value
                mag, ph_rad = cmath.polar(volt_phasor); ph_deg = math.degrees(ph_rad)
                # Format as phasor string (e.g., "5V∠0°")
                v_str = "0V" if mag < 1e-9 else f"{value_to_str(mag)}V\u2220{ph_deg:.0f}\u00B0"
                v_surf = font_to_use.render(v_str, True, volt_col_static) # Use static color

            # Position the voltage text relative to the pin
            v_off_x = PIN_R + int(3*SIZE_SCALE) # Horizontal offset from pin center
            # Position left or right based on which side of the component the pin is
            text_anchor_x = p_pin.pos[0] + v_off_x if p_pin.dx >= 0 else p_pin.pos[0] - v_off_x
            anchor_point_key = 'midleft' if p_pin.dx >= 0 else 'midright'

            if p_pin.dx == 0 and p_pin.dy == 0: # Special case for a single pin at the center (like GND)
                # Position above the pin
                v_rect = v_surf.get_rect(center=(p_pin.pos[0], p_pin.pos[1] - PIN_R - int(5*SIZE_SCALE)))
            else:
                # Position left or right, vertically centered with the pin
                v_rect = v_surf.get_rect(**{anchor_point_key: (text_anchor_x, p_pin.pos[1])})

            self.scr.blit(v_surf, v_rect)


        # --- Draw Current for Components ---
        current_phasor = None # Default current is None

        # Calculate current for passive components (R, C, L) if they have 2 terminals
        if c.ct in (CType.R,CType.C,CType.L) and len(c.pins)==2:
            p1,p2=c.pins[0],c.pins[1]; # Get the two pins (order matters for current direction later)
            v_at_p1,v_at_p2=0j,0j # Default voltage at pins is 0j

            # Get voltage phasors at the pins (if connected to non-ground nets)
            if p1.net > 0 and p1.net in net_map and net_map[p1.net] < len(node_voltages):
                 v_at_p1 = node_voltages[net_map[p1.net]]
            if p2.net > 0 and p2.net in net_map and net_map[p2.net] < len(node_voltages):
                 v_at_p2 = node_voltages[net_map[p2.net]]

            # Calculate impedance based on type and frequency
            imp = float('inf') # Default to infinite impedance (open circuit)

            if c.ct==CType.R:
                 imp = max(c.val,1e-12) # Resistance (avoid division by zero)
            elif c.ct==CType.C:
                if c.val != 0 and omega != 0:
                    imp = 1/(1j*omega*c.val) # Capacitor impedance
                else:
                    imp = float('inf') # Open circuit (C=0 or DC)
            elif c.ct==CType.L:
                if c.val != 0 and omega != 0:
                    imp = (1j*omega*c.val) # Inductor impedance
                else:
                    imp = R_L_DC_SHORT # Short circuit (L=0 or DC) - use small resistance

            # Calculate current phasor I = (V_p1 - V_p2) / Z
            # Note: This assumes current flows from pin 1 to pin 2 internally for R,C,L.
            # MNA calculates current flow differently (pos to neg for sources),
            # but for passive components connected between two nodes, this is valid.
            # The direction displayed will be from pin 1 to pin 2.
            if imp!=float('inf') and abs(imp)>1e-12:
                 current_phasor = (v_at_p1 - v_at_p2) / imp
            else:
                 current_phasor=0j # No current if impedance is infinite or voltage difference is zero

        # Get current for voltage sources (VDC, VAC)
        elif c.ct in (CType.VDC,CType.VAC):
             # The current variable for a voltage source is a direct result from MNA
            vs_idx=self.circ.vsrc_idx.get(c) # Get the index of this source in the MNA current variables
            # The current variables are located after all node voltage variables in the solution vector
            if vs_idx is not None and self.circ.node_cnt + vs_idx < len(self.circ.solution):
                current_phasor = self.circ.solution[self.circ.node_cnt + vs_idx] # Get the solved current phasor


        # Draw Current text if current_phasor was calculated for this component type
        if current_phasor is not None:
            curr_str = "" # String to display
            display_animated = self.animate_ac and self.freq > 0 # Are we animating AC?

            if display_animated:
                 # Calculate instantaneous current: |I| * cos(ωt + φ)
                 mag, ph_rad = cmath.polar(current_phasor)
                 instantaneous_i = mag * math.cos(omega * current_time + ph_rad)
                 # Format as a real number string (e.g., "12.3mA")
                 curr_str = f"{instantaneous_i:.3g}A"
                 curr_surf=font_to_use.render(curr_str,True,curr_col_animated) # Use animated color
            else:
                # Display static phasor value (Magnitude∠Phase)
                mag,ph_rad=cmath.polar(current_phasor); ph_deg=math.degrees(ph_rad)
                # Format as phasor string (e.g., "10mA∠90°")
                curr_str="0A" if mag<1e-9 else f"{value_to_str(mag)}A\u2220{ph_deg:.0f}\u00B0"
                curr_surf=font_to_use.render(curr_str,True,curr_col_static) # Use static color

            # Position the current text below the component center
            self.scr.blit(curr_surf,curr_surf.get_rect(midtop=(c.x,c.rect().bottom+int(3*SIZE_SCALE)))) # Uses global SIZE_SCALE offset


    def draw_palette(self):
        """Draws the component palette on the right side."""
        # Draw background rectangle for the palette area
        pg.draw.rect(self.scr, (45,50,55), (CANVAS_W,0,PAL_W,WIN_H))

        br = int(6*SIZE_SCALE); # Border radius for palette items
        thick = max(1,int(SIZE_SCALE)); # Border thickness

        mx,my = pg.mouse.get_pos() # Current mouse position

        # Draw each palette item
        for t_val, r_pal in self.pal_rects.items():
            is_hover = r_pal.collidepoint(mx,my) and not self.drag_type # Hover if mouse is over and not dragging something else
            is_active_drag = self.drag_type == t_val # Active if dragging this type

            # Adjust color slightly on hover or active drag
            base_c = COL[t_val]
            draw_c = [min(255,x+30) for x in base_c] if is_hover or is_active_drag else base_c

            pg.draw.rect(self.scr, draw_c, r_pal, border_radius=br) # Fill rectangle
            pg.draw.rect(self.scr, (80,80,80), r_pal, thick, border_radius=br) # Draw border

            # Draw component label ("R", "C", etc.) in the center of the item
            lbl_s = self.font.render(LBL[t_val],True,(0,0,0)) # Black text
            self.scr.blit(lbl_s, lbl_s.get_rect(center=r_pal.center))


    def draw_info_bar(self):
        """Draws the status bar and info display area at the bottom."""
        # Rectangle for the entire info bar area
        bar_r = pg.Rect(0,CANVAS_H,WIN_W,INFO_H)
        pg.draw.rect(self.scr, (28,30,32), bar_r) # Dark background color
        # Top border line
        pg.draw.line(self.scr, (60,65,70), (0,CANVAS_H), (WIN_W,CANVAS_H), max(1,int(SIZE_SCALE)))

        disp_txt = self.status # Text to display initially is the current status
        txt_c = (220,220,255) # Default text color (light blueish-white)
        cur = "" # Cursor for input prompt

        # If a prompt is active, display the prompt and buffer
        if self.prompt:
            disp_txt = self.prompt + self.input_buffer
            txt_c = (255,255,180) # Yellowish color for input text
            # Blinking cursor
            if int(pg.time.get_ticks()/400)%2==0: cur="_"

        # Add simulation status (AC/DC and frequency/animation state) if solved
        elif self.circ.solution is not None and self.circ.solution.size > 0:
            sim_type = "AC" if self.freq > 0 else "DC"
            if sim_type == "AC":
                if self.animate_ac:
                     disp_txt += f"    |    {sim_type} (Animating @ {value_to_str(self.freq)}Hz, t={self.sim_time:.2f}s)"
                else:
                     disp_txt += f"    |    {sim_type} (Steady-State Phasors @ {value_to_str(self.freq)}Hz)"
            else: # DC analysis (freq == 0)
                 disp_txt += f"    |    {sim_type} Analysis"

        # Render and blit the status text
        status_s = self.font.render(disp_txt+cur, True, txt_c)
        # Position status text in the bottom left corner of the bar
        self.scr.blit(status_s, (int(15*SIZE_SCALE), CANVAS_H + int(5*SIZE_SCALE)))


        # Draw detailed info for the current info page (if not page 0)
        if self.show_page > 0 and self.circ.Y is not None and self.circ.Y.shape[0] > 0:
            # Mapping of page index to page name
            page_names = {1:"Heatmap", 2:"Y*x=b Text", 3:"Circuit Info"}
            if self.show_page in page_names:
                # Calculate the page number display (excluding page 0 from count)
                visible_pages = [1, 2, 3] # Pages that have content in the info bar
                page_index_in_visible = visible_pages.index(self.show_page) + 1 if self.show_page in visible_pages else self.show_page # +1 because humans start counting from 1
                total_visible_pages = len(visible_pages)

                page_title = f"Pg {page_index_in_visible}/{total_visible_pages}: {page_names[self.show_page]}"
                title_s = self.font_small.render(page_title, True, (150,150,150)) # Gray title text
                # Position page title in the bottom right corner of the bar, below status text
                self.scr.blit(title_s, title_s.get_rect(topright=(WIN_W-int(10*SIZE_SCALE), CANVAS_H + int(25*SIZE_SCALE))))

                # Call the specific draw function for the current info page
                draw_func = {1:self.draw_heat_map, 2:self.draw_matrix_text, 3:self.draw_stamps_info}.get(self.show_page)
                if draw_func:
                    draw_func() # Execute the drawing function


    def draw_legend(self):
        """Draws the full legend overlay."""
        # ── Section data ──────────────────────────────────────────────────────────
        sections = {
            "Mouse": [("Drag Palette", "Place component from palette"),
                      ("Pin\u2192Pin", "Wire between component pins"),
                      ("Drag Part",    "Move selected component(s)"),
                      ("R-Click",      "Delete component"),
                      ("Wheel",        "Adjust component value \u00d710 or \u00f710"),
                      ("Sh+Whl",       "Adjust component value \u00d7\u221a10 or \u00f7\u221a10"),
                      ("Ctrl+Click",   "Toggle selection for multi-select")],
            "Keyboard": [("S",             "Solve circuit"),
                         ("C",             "Toggle component value edit mode"),
                         ("F",             "Set simulation frequency"),
                         ("W",             "Toggle wire delete mode"),
                         ("A",             "Toggle AC steady-state animation"),
                         ("M",             "Cycle info pages (Overlay, Heatmap, Y\u00d7=b, Info)"),
                         ("E",             "Go directly to Circuit Info page"),
                         ("U",             "Toggle LU factorization popup view"),
                         ("Del / Backspace","Delete selected component(s)"),
                         ("Ctrl+A",        "Select all components"),
                         ("Esc",           "Cancel / close overlay"),
                         ("T / L / H / ?", "Toggle this Legend")],
            "Info Pages (M/E keys)": [("0", "Solution Overlay (default)"),
                                      ("1", "MNA Matrix Heatmap"),
                                      ("2", "Y\u00d7=b Matrix Text View"),
                                      ("3", "Circuit Component & Net Info")]
        }

        # ── Colours, fonts, geometry ─────────────────────────────────────────────
        legend_bg = (25, 30, 35, 235)
        hdr_c     = (180, 220, 255)
        txt_c     = (220, 220, 220)
        key_c     = (255, 255, 180)
        bord_c    = (100, 120, 140, 220)

        pad      = int(20 * SIZE_SCALE)
        txt_sp   = int(5  * SIZE_SCALE)
        sect_sp  = int(15 * SIZE_SCALE)
        br       = int(12 * SIZE_SCALE)
        bord_th  = max(1, int(2 * SIZE_SCALE))

        f_leg   = self.font_small
        f_hdr   = self.font
        f_title = pg.font.SysFont(FONT_PRIMARY_NAME, int(18 * FONT_SCALE), bold=True)

        # ── Pre-compute sizes ────────────────────────────────────────────────────
        max_kw   = 0            # key column width
        max_dw   = 0            # description column width
        sect_hs  = []           # per-section heights

        line_h_leg = f_leg.get_linesize() + txt_sp
        hdr_h_leg  = f_hdr.get_linesize() + txt_sp

        for sn, itms in sections.items():
            sh = hdr_h_leg                   # section height starts with header
            current_max_dw_section = 0

            for k, d_orig in itms:
                max_kw = max(max_kw, f_leg.render(k + ":", 1, key_c).get_width())

                # wrap description, track width & lines
                words = d_orig.split(' ')
                current_line_desc = ""
                desc_lines_count  = 1
                temp_max_desc_w_for_item = 0
                max_desc_render_width = int(300 * SIZE_SCALE)

                for word in words:
                    test_line_desc = current_line_desc + word + " "
                    if (f_leg.render(test_line_desc, 1, txt_c).get_width() >
                        max_desc_render_width and current_line_desc):
                        temp_max_desc_w_for_item = max(
                            temp_max_desc_w_for_item,
                            f_leg.render(current_line_desc.strip(), 1, txt_c).get_width()
                        )
                        current_line_desc = word + " "
                        desc_lines_count += 1
                    else:
                        current_line_desc = test_line_desc

                if current_line_desc.strip():
                    temp_max_desc_w_for_item = max(
                        temp_max_desc_w_for_item,
                        f_leg.render(current_line_desc.strip(), 1, txt_c).get_width()
                    )
                else:
                    temp_max_desc_w_for_item = 0  # empty description

                current_max_dw_section = max(current_max_dw_section,
                                             temp_max_desc_w_for_item)

                sh += line_h_leg * desc_lines_count

            max_dw = max(max_dw, current_max_dw_section)
            sect_hs.append(sh)

        total_th = sum(sect_hs) + sect_sp * (len(sect_hs) - 1) if sect_hs else 0

        max_kw = max(max_kw, int(50 * SIZE_SCALE))
        max_dw = max(max_dw, int(100 * SIZE_SCALE))

        itm_ind = int(10 * SIZE_SCALE)
        k_d_gap = int(10 * SIZE_SCALE)

        content_text_width = itm_ind + max_kw + k_d_gap + max_dw

        title_h_leg = f_title.get_linesize() + txt_sp
        leg_w = min(WIN_W - pad * 2, content_text_width + pad * 2)
        leg_h = total_th + title_h_leg + pad * 2

        leg_r = pg.Rect(
            (WIN_W - leg_w) // 2,
            max(int(10 * SIZE_SCALE), (WIN_H - leg_h) // 2),
            leg_w,
            min(WIN_H - int(20 * SIZE_SCALE), leg_h)
        )

        # ── Surface & frame ──────────────────────────────────────────────────────
        leg_s = pg.Surface(leg_r.size, pg.SRCALPHA)
        pg.draw.rect(leg_s, legend_bg, leg_s.get_rect(), border_radius=br)
        pg.draw.rect(leg_s, bord_c,  leg_s.get_rect(), bord_th, border_radius=br)

        # ── Title ────────────────────────────────────────────────────────────────
        title_s = f_title.render("CONTROLS & REFERENCE", 1, (255, 255, 255))
        leg_s.blit(title_s, title_s.get_rect(midtop=(leg_w // 2, pad)))

        # ── Sections & items ─────────────────────────────────────────────────────
        cur_y = title_s.get_rect(midtop=(leg_w // 2, pad)).bottom + pad // 2
        x_st  = pad

        for sn, itms in sections.items():
            hdr_s = f_hdr.render(sn, 1, hdr_c)
            leg_s.blit(hdr_s, (x_st, cur_y))
            cur_y += hdr_h_leg

            for k, d_orig in itms:
                key_s = f_leg.render(k + ":", 1, key_c)
                leg_s.blit(key_s, (x_st + itm_ind, cur_y))

                # wrapped description
                words = d_orig.split(' ')
                current_line_desc = ""
                line_start_x_desc = x_st + itm_ind + max_kw + k_d_gap
                temp_y_desc = cur_y
                available_desc_width = leg_w - line_start_x_desc - pad

                for word in words:
                    test_line_desc = current_line_desc + word + " "
                    if (f_leg.render(test_line_desc, 1, txt_c).get_width() >
                        available_desc_width and current_line_desc):
                        desc_s = f_leg.render(current_line_desc.strip(), 1, txt_c)
                        leg_s.blit(desc_s, (line_start_x_desc, temp_y_desc))
                        temp_y_desc += line_h_leg
                        current_line_desc = word + " "
                    else:
                        current_line_desc = test_line_desc

                if current_line_desc.strip():
                    desc_s = f_leg.render(current_line_desc.strip(), 1, txt_c)
                    leg_s.blit(desc_s, (line_start_x_desc, temp_y_desc))
                    cur_y = temp_y_desc + line_h_leg
                else:
                    cur_y += line_h_leg

            cur_y += sect_sp // 2

        # ── Blit & close hint ───────────────────────────────────────────────────
        self.scr.blit(leg_s, leg_r)

        if leg_r.bottom < WIN_H - int(10 * SIZE_SCALE):
            close_s = self.font_tiny.render(
                "Press T, L, H, ?, or Esc to close", 1, (180, 180, 180)
            )
            self.scr.blit(
                close_s,
                close_s.get_rect(midbottom=(WIN_W // 2, WIN_H - int(5 * SIZE_SCALE)))
            )


    def draw_heat_map(self):
        """Draws a heatmap representation of the MNA matrix and RHS vector."""
        Y = self.circ.Y; # MNA matrix
        b_vec = self.circ.b; # Right-hand side vector

        if Y is None or b_vec is None or Y.shape[0] == 0:
            # Display a message if no matrix data is available
            no_data_surf = self.font_small.render("No matrix data. Solve circuit (S).", True, (200,200,200))
            self.scr.blit(no_data_surf, no_data_surf.get_rect(centerx=WIN_W//2, top=CANVAS_H + int(50*SIZE_SCALE)))
            return

        dim = Y.shape[0] # Dimension of the matrix/vector

        # Calculate positioning and cell size based on available space and matrix dimension
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE) # Start position below status text
        available_h = WIN_H - info_bar_content_y_start - int(5 * SIZE_SCALE) # Available height in the info bar area
        available_w = WIN_W - int(40 * SIZE_SCALE) # Available width (with some padding)

        # Determine optimal cell size
        max_cell_s = int(12*SIZE_SCALE) # Maximum allowed cell size
        # Calculate cell size based on available space divided by dimension (+ a little extra for 'b' vector)
        # Ensure division by zero is handled
        cell_s = max(1, min(available_w // (dim + 1.5) if (dim+1.5)>0 else available_w,
                             available_h // dim if dim>0 else available_h,
                             max_cell_s)) # Take the minimum of width/height constraints and max size

        spacing = max(1, int(cell_s / 8)) # Small spacing between cells

        # Calculate total dimensions of the drawn matrix block (Y and b)
        total_h = dim * (cell_s + spacing) - spacing if dim > 0 else 0 # Height of the matrix
        total_w_y = dim * (cell_s + spacing) - spacing if dim > 0 else 0 # Width of Y matrix block
        total_w_b = cell_s # Width of b vector block
        b_gap = spacing * 2 # Gap between Y and b blocks

        # Calculate origin (top-left) position to center the heatmap block
        ox = (WIN_W - total_w_y - b_gap - total_w_b) // 2
        oy = info_bar_content_y_start + (available_h - total_h) // 2

        # Determine maximum absolute values for normalization (scaling colors)
        mY = np.max(np.abs(Y)) if np.any(Y) else 1.0; mY = max(mY, 1e-12) # Max value in Y (avoid division by zero)
        mb = np.max(np.abs(b_vec)) if np.any(b_vec) else 1.0; mb = max(mb, 1e-12) # Max value in b (avoid division by zero)

        border_th = max(1, int(0.5*SIZE_SCALE)) # Border thickness for cells

        # Draw Y matrix heatmap
        for r_iter in range(dim):
            for c_iter in range(dim):
                v = np.abs(Y[r_iter,c_iter]); # Get magnitude of the element
                br_val = min(1, v/mY if mY!=0 else 0) # Normalize to 0-1 range
                col_val = int(220*br_val) # Scale color intensity (from 0 to 220)

                # Calculate rectangle for this cell
                cell_r = pg.Rect(ox+c_iter*(cell_s+spacing), oy+r_iter*(cell_s+spacing), cell_s, cell_s)
                # Draw filled rectangle (grayscale based on magnitude)
                pg.draw.rect(self.scr, (col_val,col_val,col_val), cell_r)
                # Draw cell border if cells are large enough
                if cell_s > 3:
                     pg.draw.rect(self.scr, (80,80,80), cell_r, border_th) # Dark gray border


        # Draw b vector heatmap
        bx_pos = ox + total_w_y + b_gap # X position for the b vector block
        for r_idx in range(dim):
             v_val_b = b_vec[r_idx]; # Get magnitude of the element
             v_abs = np.abs(v_val_b);
             br_val_b = min(1, v_abs/mb if mb!=0 else 0) # Normalize to 0-1 range
             col_val_b = int(220*br_val_b) # Scale color intensity

             # Calculate rectangle for this cell
             cell_r_b = pg.Rect(bx_pos, oy+r_idx*(cell_s+spacing), cell_s, cell_s)
             # Draw filled rectangle (greenish color based on magnitude)
             pg.draw.rect(self.scr, (int(col_val_b*0.6),col_val_b,int(col_val_b*0.6)), cell_r_b)
             # Draw cell border if cells are large enough
             if cell_s > 3:
                  pg.draw.rect(self.scr, (80,80,80), cell_r_b, border_th)


    def draw_matrix_text(self):
        """Draws a text representation of the MNA system Y*x=b."""
        Y, b_vec = self.circ.Y, self.circ.b
        if Y is None or b_vec is None or Y.shape[0] == 0:
            # Display a message if no matrix data
            no_data_surf = self.font_small.render("No matrix data. Solve circuit (S).", True, (200,200,200))
            self.scr.blit(no_data_surf, no_data_surf.get_rect(centerx=WIN_W//2, top=CANVAS_H + int(50*SIZE_SCALE)))
            return

        dim = Y.shape[0] # Dimension of the matrix/vector

        # Helper function to format complex numbers for text display
        def fmt_c(val):
            mag = abs(val)
            # Use a small threshold for considering a value zero
            if mag < 1e-12: return "0"
            # If imaginary part is much smaller than real, show only real
            if abs(val.imag) < 1e-9 * max(abs(val.real), 1e-12): return f"{val.real:.2g}"
            # If real part is much smaller than imaginary, show only imaginary
            if abs(val.real) < 1e-9 * max(abs(val.imag), 1e-12): return f"{val.imag:.2g}j"
            # Otherwise, show in polar form (Magnitude∠Phase)
            return f"{mag:.2g}\u2220{math.degrees(cmath.phase(val)):.0f}\u00B0" # Using unicode angle symbol

        # Get the labels for the unknown variables (voltages and currents)
        unknown_labels = [self.circ.idx_to_nodename.get(i, f"x{i}") for i in range(dim)]

        # Format all elements of Y and b as strings
        Y_str_elements = [[fmt_c(Y[r, c_idx]) for c_idx in range(dim)] for r in range(dim)]
        b_str_elements = [fmt_c(b_val) for b_val in b_vec]

        # Calculate required width for each column based on the widest string in that column (including headers)
        # Initialize max_row_label_w for the row labels column
        max_row_label_w = 0
        if unknown_labels:
             max_row_label_w = max(self.font_tiny.render(lbl + ": ", True, (0,0,0)).get_width() for lbl in unknown_labels)

        Y_col_widths = [0] * dim
        for c_idx in range(dim):
            # Start with the width of the column header (unknown variable label)
            header_w = self.font_tiny.render(unknown_labels[c_idx], True, (0,0,0)).get_width()
            Y_col_widths[c_idx] = header_w
            # Check all elements in the column
            for r_idx in range(dim):
                Y_col_widths[c_idx] = max(Y_col_widths[c_idx], self.font_tiny.render(Y_str_elements[r_idx][c_idx], True, (0,0,0)).get_width())

        # Calculate required width for the 'x' vector and 'b' vector columns
        x_vec_col_width = 0
        if unknown_labels:
             x_vec_col_width = max(self.font_tiny.render(lbl, True, (0,0,0)).get_width() for lbl in unknown_labels)
        b_vec_col_width = 0
        if b_str_elements:
             b_vec_col_width = max(self.font_tiny.render(s, True, (0,0,0)).get_width() for s in b_str_elements)

        # Define spacing between columns and operators
        col_spacing = int(10 * SIZE_SCALE); # Space between matrix columns
        op_spacing = int(8 * SIZE_SCALE); # Space around operators (*, =)
        bracket_spacing = int(3*SIZE_SCALE); # Space between matrix/vector content and brackets

        line_h = self.font_tiny.get_linesize() + int(2 * SIZE_SCALE); # Height of a line of text plus spacing

        # Starting position for the text display
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE);
        start_x = int(10 * SIZE_SCALE); # Left padding
        current_y = info_bar_content_y_start; # Current drawing Y position


        # Draw column headers (unknown variable labels above Y and x vectors)
        current_x_col_header = start_x + max_row_label_w + bracket_spacing # Start position for the first Y column header
        for c_idx in range(dim):
            header_surf = self.font_tiny.render(unknown_labels[c_idx], True, (180, 180, 220)); # Header text color
            # Position header centered over its column
            header_rect = header_surf.get_rect(centerx=current_x_col_header + Y_col_widths[c_idx] // 2, top=current_y);
            self.scr.blit(header_surf, header_rect);
            current_x_col_header += Y_col_widths[c_idx] + col_spacing; # Move x position to the next column header

        # Position header for the x vector column
        current_x_x_header = current_x_col_header - col_spacing + bracket_spacing + op_spacing + self.font_tiny.render("*", True, (0,0,0)).get_width() + op_spacing + bracket_spacing
        if dim > 0: # Only draw if there's at least one unknown
             x_header_surf = self.font_tiny.render("x", True, (180, 180, 220)) # Generic label for x vector
             x_header_rect = x_header_surf.get_rect(centerx=current_x_x_header + x_vec_col_width // 2, top=current_y)
             self.scr.blit(x_header_surf, x_header_rect)


        current_y += line_h; # Move Y down for the matrix content


        # Draw each row of the matrix equation Y * x = b
        for r_idx in range(dim):
            # Stop drawing if we run out of space in the info bar
            if current_y + line_h > WIN_H - int(5*SIZE_SCALE): break

            current_x_row = start_x; # Start position for the current row

            # Draw row label (unknown variable label before the matrix)
            row_label_surf = self.font_tiny.render(unknown_labels[r_idx] + ":", True, (180, 180, 220));
            row_label_rect = row_label_surf.get_rect(topleft=(current_x_row, current_y));
            self.scr.blit(row_label_surf, row_label_rect);
            current_x_row += max_row_label_w; # Move x position past the row label

            # Draw opening bracket for Y matrix
            bracket_surf = self.font_tiny.render("|", True, (200,210,200)); # Bracket color
            self.scr.blit(bracket_surf, (current_x_row, current_y));
            current_x_row += bracket_surf.get_width() + bracket_spacing; # Move x past bracket and add spacing

            # Draw Y matrix elements for this row
            for c_idx in range(dim):
                # Stop drawing columns if we run out of horizontal space (unlikely given calculation, but safe)
                if current_x_row + Y_col_widths[c_idx] + col_spacing > WIN_W - int(10*SIZE_SCALE):
                     # Indicate truncation if not all columns are shown
                     if c_idx < dim -1:
                          dot_surf = self.font_tiny.render("...", True, (200,210,200))
                          self.scr.blit(dot_surf, (current_x_row, current_y + line_h // 2 - dot_surf.get_height() // 2 -1))
                     break # Stop drawing columns in this row

                elem_surf = self.font_tiny.render(Y_str_elements[r_idx][c_idx], True, (200,210,200)); # Element text color
                # Position element text right-aligned within its column width, vertically centered
                elem_rect = elem_surf.get_rect(right=current_x_row + Y_col_widths[c_idx], centery=current_y + line_h // 2 -1); # Adjust centery slightly for visual alignment
                self.scr.blit(elem_surf, elem_rect);
                current_x_row += Y_col_widths[c_idx] + col_spacing; # Move x position to the next element

            # Draw closing bracket for Y matrix
            # Need to step back by the last col_spacing before drawing the bracket
            current_x_row -= col_spacing;
            self.scr.blit(bracket_surf, (current_x_row, current_y));
            current_x_row += bracket_surf.get_width() + op_spacing; # Move x past bracket and add operator spacing

            # Draw multiplication operator "*"
            mul_surf = self.font_tiny.render("*", True, (220,220,100)); # Operator color
            self.scr.blit(mul_surf, (current_x_row, current_y));
            current_x_row += mul_surf.get_width() + op_spacing; # Move x past operator and add spacing

            # Draw opening bracket for x vector
            self.scr.blit(bracket_surf, (current_x_row, current_y));
            current_x_row += bracket_surf.get_width() + bracket_spacing; # Move x past bracket and add spacing

            # Draw x vector element for this row (the unknown variable label)
            x_elem_surf = self.font_tiny.render(unknown_labels[r_idx], True, (200,210,200));
            # Position element text centered within its column width, vertically centered
            x_elem_rect = x_elem_surf.get_rect(centerx=current_x_row + x_vec_col_width // 2, centery=current_y + line_h // 2 -1);
            self.scr.blit(x_elem_surf, x_elem_rect);
            current_x_row += x_vec_col_width; # Move x position past the x element

            # Draw closing bracket for x vector
            self.scr.blit(bracket_surf, (current_x_row, current_y));
            current_x_row += bracket_spacing + op_spacing; # Move x past bracket and add operator spacing


            # Draw equals operator "="
            eq_surf = self.font_tiny.render("=", True, (220,220,100)); # Operator color
            self.scr.blit(eq_surf, (current_x_row, current_y));
            current_x_row += eq_surf.get_width() + op_spacing; # Move x past operator and add spacing

            # Draw opening bracket for b vector
            self.scr.blit(bracket_surf, (current_x_row, current_y));
            current_x_row += bracket_spacing; # Move x past bracket and add spacing

            # Draw b vector element for this row
            b_elem_surf = self.font_tiny.render(b_str_elements[r_idx], True, (180,220,180)); # b element text color (greenish)
            # Position element text centered within its column width, vertically centered
            b_elem_rect = b_elem_surf.get_rect(centerx=current_x_row + b_vec_col_width // 2, centery=current_y + line_h // 2 -1);
            self.scr.blit(b_elem_surf, b_elem_rect);
            current_x_row += b_vec_col_width; # Move x position past the b element

            # Draw closing bracket for b vector
            self.scr.blit(bracket_surf, (current_x_row, current_y));

            current_y += line_h; # Move Y down for the next row


    def draw_stamps_info(self):
        """Draws detailed circuit information (nodes, sources, component details)."""
        if self.circ.Y is None or self.circ.Y.shape[0] == 0: return # No solution data

        lines = [] # List to hold lines of text to display

        # Add summary info
        lines.append(f"Nodes (Excl.GND): {self.circ.node_cnt}, V-Sources: {len(self.circ.vsrc_idx)}, Matrix Dim: {self.circ.Y.shape[0]}")
        lines.append(f"Frequency: {value_to_str(self.freq)}Hz")
        # Display the loaded simulation constants
        lines.append(f"Gmin: {GMIN_DEFAULT:.1e}, L=0/DC R: {R_L_DC_SHORT:.1e}")
        lines.append("-" * 40) # Separator line

        # Add solved variable values
        lines.append("Solved Variables (MNA Unknowns):")
        if self.circ.solution is not None:
            for i in range(self.circ.Y.shape[0]):
                var_name = self.circ.idx_to_nodename.get(i, f"Var{i}") # Get the variable name (Voltage or Current)
                val_str = ""
                if i < len(self.circ.solution): # Ensure index is within bounds
                    phasor = self.circ.solution[i]
                    mag, ph_rad = cmath.polar(phasor); ph_deg = math.degrees(ph_rad) # Get magnitude and phase
                    unit = "(V)" if var_name.startswith("V(") else "(A)" # Determine unit based on name format
                    # Format the value string
                    val_str = " = 0" if mag < 1e-12 else f" = {value_to_str(mag)}\u2220{ph_deg:.0f}\u00B0 {unit}" # Use value_to_str for magnitude
                lines.append(f"  x[{i}]: {var_name}{val_str}") # Add line for each unknown
        else:
             lines.append("  (Solution not available - solver error?)")


        lines.append("-" * 40) # Separator line

        # Add component details and net connections
        lines.append("Component Details & Connections:")
        omega = 2 * math.pi * self.freq # Angular frequency for impedance calcs

        for c_comp in self.circ.comps:
            pin_info_list = []
            # Get connection info for each pin
            for p_pin in c_comp.pins:
                net_display = "GND" if p_pin.net == 0 else f"N{p_pin.net}" # Display net ID (use "GND" for net 0)
                pin_info_list.append(f"{p_pin.name}:{net_display}") # Format "PinName:NetID"
            nets_str = ", ".join(pin_info_list) # Combine pin info strings

            details_str = "" # String for component specific details
            try:
                if c_comp.ct == CType.GND: details_str = "Ground Reference (Net 0)"
                elif c_comp.ct == CType.R:
                    details_str = f"R={c_comp.label()}, G={1/max(1e-12,c_comp.val):.2g}S" # Display Resistance and Conductance
                elif c_comp.ct == CType.C:
                    if c_comp.val == 0: details_str=f"C={c_comp.label()}, Yc=0S (Open)"
                    elif omega==0: details_str=f"C={c_comp.label()}, Yc=0S (DC Open)"
                    else: details_str=f"C={c_comp.label()}, Yc={1j*omega*c_comp.val:.2g}S" # Display Capacitance and Admittance
                elif c_comp.ct == CType.L:
                    if c_comp.val == 0: details_str=f"L={c_comp.label()}, Yl={1/R_L_DC_SHORT:.2g}S (0H Short)" # Use loaded R_L_DC_SHORT
                    elif omega==0: details_str=f"L={c_comp.label()}, Yl={1/R_L_DC_SHORT:.2g}S (DC Short)" # Use loaded R_L_DC_SHORT
                    else: details_str=f"L={c_comp.label()}, Yl={1/(1j*omega*c_comp.val):.2g}S" # Display Inductance and Admittance
                elif c_comp.ct in (CType.VDC,CType.VAC):
                    details_str=f"V={c_comp.label()}" # Display Source Voltage
            except (ZeroDivisionError, OverflowError):
                 details_str="Math Error in calc" # Catch potential math errors (e.g., division by zero)

            lines.append(f"  {c_comp.name} [{nets_str}]: {details_str}") # Add line for each component

        # Calculate line height and available drawing area
        line_h_info = self.font_tiny.get_linesize()
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE)
        start_y_info = info_bar_content_y_start + int(2*SIZE_SCALE) # Start below the page title
        left_pad_info = int(10*SIZE_SCALE) # Left padding

        # Determine how many lines can fit in the available height
        max_fit_lines = (WIN_H - start_y_info - int(5*SIZE_SCALE)) // line_h_info if line_h_info > 0 else 0

        # Draw the text lines, truncating if necessary
        for i, line_txt in enumerate(lines[:max_fit_lines]):
            surf = self.font_tiny.render(line_txt, True, (180, 220, 180)) # Greenish text color
            self.scr.blit(surf, (left_pad_info, start_y_info + i * line_h_info))


    def draw_lu_popup_overlay(self):
        """Draws the LU decomposition matrices in a popup window."""
        # Check if LU factors are available
        if self.circ.Y is None or self.circ.P_matrix is None or \
           self.circ.L_matrix is None or self.circ.U_matrix is None:
            no_data_surf = self.font_small.render("LU factors not available. Solve circuit first (S).", True, (200,50,50)) # Error message color
            self.scr.blit(no_data_surf, no_data_surf.get_rect(center=(WIN_W//2, WIN_H//2))) # Position centered
            return

        # Calculate popup dimensions and position
        popup_w = int(WIN_W * 0.9); # 90% of window width
        popup_h = int(WIN_H * 0.85); # 85% of window height
        popup_x = (WIN_W - popup_w) // 2; # Centered horizontally
        popup_y = (WIN_H - popup_h) // 2; # Centered vertically

        # Create a surface for the popup with alpha channel
        popup_surf = pg.Surface((popup_w, popup_h), pg.SRCALPHA);

        # Define colors for the popup elements
        popup_bg_col = (35, 40, 45, 245); # Dark background, slightly transparent
        border_col = (100, 120, 140, 230); # Border color
        title_col = (200, 200, 255); # Light blueish title
        matrix_title_col = (190, 190, 230); # Slightly darker matrix titles
        expl_text_col = (170, 170, 170); # Gray explanation text
        text_col = (210, 210, 210); # Light gray matrix element text
        header_col = (180, 180, 220); # Node/Column header text

        # Draw background and border for the popup surface
        pg.draw.rect(popup_surf, popup_bg_col, popup_surf.get_rect(), border_radius=10);
        pg.draw.rect(popup_surf, border_col, popup_surf.get_rect(), 2, border_radius=10); # 2px border

        # Define vertical spacing
        padding = int(15 * SIZE_SCALE); # Padding around content
        line_h_matrix_elem = self.font_tiny.get_linesize() + int(1*SIZE_SCALE); # Line height for matrix elements
        line_h_explanation = self.font_tiny.get_linesize() + int(1*SIZE_SCALE); # Line height for explanation text
        line_h_title = self.font_small.get_linesize() + int(2*SIZE_SCALE); # Line height for matrix titles

        current_y = padding; # Starting Y position for content within the popup surface

        # Draw main title
        main_title_surf = self.font.render("LU Decomposition: P Y = L U", True, title_col);
        popup_surf.blit(main_title_surf, main_title_surf.get_rect(centerx=popup_w//2, top=current_y));
        current_y += self.font.get_linesize() + padding // 2; # Move Y down after title

        # Explanation text about LU decomposition
        explanation_intro = [
            "The MNA matrix Y is decomposed for solving Yx=b.",
            "P is a Permutation matrix that reorders rows of Y for numerical stability (pivoting).",
            "L is Lower triangular (ones on diagonal), U is Upper triangular.",
            "The factorization achieved is: P @ Y = L @ U.",
            "Then, L@(U@x) = P@b is solved via: 1) L@z = P@b (for z), 2) U@x = z (for x)."
        ]
        # Draw explanation lines
        for line_text in explanation_intro:
            # Check if there's enough space before drawing the next line
            if current_y + line_h_explanation > popup_h - padding * 3: break
            expl_surf = self.font_tiny.render(line_text, True, expl_text_col);
            popup_surf.blit(expl_surf, (padding, current_y));
            current_y += line_h_explanation; # Move Y down

        current_y += padding // 2; # Add extra space after explanation

        # Data for the matrices to be displayed (Title, Matrix Data, Explanation, Matrix Type String for fmt_val)
        mats_to_display_data = [
            ("Original Y Matrix (MNA System Matrix):", self.circ.Y, "This is the matrix to be solved.", 'Y'),
            ("Permutation Matrix P (from Pivoting):", self.circ.P_matrix, "P reorders rows of Y: P@Y is factored.", 'P'),
            ("Lower Triangular Factor L (Diagonal = 1):", self.circ.L_matrix, "L is a factor of the permuted Y matrix.", 'L'),
            ("Upper Triangular Factor U:", self.circ.U_matrix, "U is the other factor: P@Y = L@U.", 'U')
        ]

        dim = self.circ.Y.shape[0] # Dimension of the matrices (assuming they are square and same size)
        unknown_labels = [self.circ.idx_to_nodename.get(i, f"x{i}") for i in range(dim)] # Labels for rows/columns

        # Helper function to format matrix elements for the popup display
        # Added matrix_type argument to determine formatting (e.g., for P)
        def fmt_val_popup(val, matrix_type): # Renamed to avoid conflict if needed later, clearer scope
            if val is None: return "---"
            # Use matrix_type to decide formatting specifics
            if matrix_type == 'P': return f"{int(round(val.real))}" # P matrix elements are integers (0 or 1)
            # Generic complex formatting for Y, L, U
            mag = abs(val)
            if mag < 1e-11: return "0" # Treat very small values as zero
            # Use 1 significant figure for real/imaginary parts if one is dominant
            if abs(val.imag) < 1e-9 * max(abs(val.real),1e-11): return f"{val.real:.1g}"
            if abs(val.real) < 1e-9 * max(abs(val.imag),1e-11): return f"{val.imag:.1g}j"
            # Otherwise, show in polar form with 1 significant figure for magnitude
            return f"{mag:.1g}\u2220{math.degrees(cmath.phase(val)):.0f}" # Unicode angle, 0 decimal places for phase


        # Estimate how many columns can fit based on a typical element width
        test_char_width = self.font_tiny.render("-0.0X",0,(0,0,0)).get_width(); # Width of a sample string
        # Available width for matrices = popup width - padding - a bit extra for column headers/labels/spacing
        max_cols_to_show_estimate = (popup_w - padding * 3 - int(40*SIZE_SCALE)) // (test_char_width if test_char_width > 0 else 10);
        max_cols_to_show = max(1, min(dim, max_cols_to_show_estimate)); # Cap by dimension and minimum 1


        # Estimate how many rows can fit per matrix section based on available height
        remaining_h_for_matrices = popup_h - current_y - padding - self.font_small.get_linesize(); # Remaining height
        # Average overhead per matrix block (title + explanation + one line of header)
        avg_overhead_per_matrix = line_h_title + line_h_explanation + line_h_matrix_elem;
        if len(mats_to_display_data) > 0 and line_h_matrix_elem > 0:
            max_rows_per_matrix_est = (remaining_h_for_matrices // len(mats_to_display_data) - avg_overhead_per_matrix) // line_h_matrix_elem;
        else:
            max_rows_per_matrix_est = 1;
        max_rows_to_show_per_matrix = max(1, min(dim, max_rows_per_matrix_est if max_rows_per_matrix_est > 0 else 1)); # Cap by dimension and minimum 1


        # Draw each matrix section
        # Added matrix_type to the loop unpacking
        for title_str, matrix_data, expl_str, matrix_type in mats_to_display_data:
            # Check if there's enough space for this matrix block
            if current_y + line_h_title + line_h_explanation + line_h_matrix_elem * 2 > popup_h - padding - self.font_small.get_linesize():
                 break # Stop drawing matrices if no space

            # Draw matrix title
            mat_title_surf = self.font_small.render(title_str, True, matrix_title_col);
            popup_surf.blit(mat_title_surf, (padding, current_y));
            current_y += line_h_title; # Move Y down

            # Draw matrix explanation
            expl_surf = self.font_tiny.render(expl_str, True, expl_text_col);
            popup_surf.blit(expl_surf, (padding + int(10*SIZE_SCALE), current_y)); # Indent explanation
            current_y += line_h_explanation + int(2*SIZE_SCALE); # Move Y down, add slight gap

            # If matrix data is not available (e.g., solve failed), display N/A
            if matrix_data is None:
                popup_surf.blit(self.font_tiny.render("N/A", True, text_col), (padding + int(10*SIZE_SCALE), current_y));
                current_y += line_h_matrix_elem * 2; # Add space for roughly two lines
                continue # Skip to next matrix

            # Determine actual number of rows/columns to show based on available space and matrix size
            current_max_rows = min(dim, max_rows_to_show_per_matrix)
            current_max_cols = min(dim, max_cols_to_show)

            # Format the elements that will be displayed as strings
            # Pass the matrix_type string to the formatting function
            str_elements = [[fmt_val_popup(matrix_data[r,c], matrix_type) for c in range(current_max_cols)] for r in range(current_max_rows)]

            # Calculate column widths for the *displayed* part of the matrix
            col_widths = [0] * current_max_cols
            row_label_max_w = 0

            # Calculate width for row labels and column headers
            for i_calc in range(max(current_max_rows, current_max_cols)): # Iterate up to the max of rows/cols being shown
                 if i_calc < current_max_rows:
                    # Row labels are unknown names for Y, generic r# for L/U/P
                    row_label_text_calc = unknown_labels[i_calc] if title_str.startswith("Original Y") else f"r{i_calc}"
                    row_label_max_w = max(row_label_max_w, self.font_tiny.render(row_label_text_calc + ":", True, (0,0,0)).get_width())
                 if i_calc < current_max_cols:
                    # Column headers are unknown names for Y, generic c# for L/U/P
                    col_header_text_calc = unknown_labels[i_calc] if title_str.startswith("Original Y") else f"c{i_calc}"
                    col_widths[i_calc] = self.font_tiny.render(col_header_text_calc, True, (0,0,0)).get_width()
                    # Update column width based on element strings
                    for r_idx_calc in range(current_max_rows):
                        col_widths[i_calc] = max(col_widths[i_calc], self.font_tiny.render(str_elements[r_idx_calc][i_calc],True,(0,0,0)).get_width())

            x_offset_matrix = padding + row_label_max_w + int(5*SIZE_SCALE); # Start position for the matrix content (after row labels)

            # Draw column headers above the matrix
            temp_x = x_offset_matrix
            for c_idx_draw in range(current_max_cols):
                col_header_text_draw = unknown_labels[c_idx_draw] if title_str.startswith("Original Y") else f"c{c_idx_draw}"
                h_surf = self.font_tiny.render(col_header_text_draw, True, header_col);
                # Position header centered over its calculated column width
                popup_surf.blit(h_surf, h_surf.get_rect(centerx=temp_x + col_widths[c_idx_draw]//2, top=current_y));
                temp_x += col_widths[c_idx_draw] + int(3*SIZE_SCALE); # Move x to the next column header

            # Indicate if columns are truncated
            if dim > current_max_cols:
                popup_surf.blit(self.font_tiny.render("...", True, header_col), (temp_x, current_y)); # Draw ellipsis

            current_y += line_h_matrix_elem; # Move Y down for matrix elements

            # Draw matrix elements row by row
            for r_idx_draw in range(current_max_rows):
                # Check space before drawing row
                if current_y + line_h_matrix_elem > popup_h - padding - self.font_small.get_linesize():
                     # Indicate truncation if not all rows are shown
                     if r_idx_draw < dim - 1:
                          popup_surf.blit(self.font_tiny.render("... (more rows)", True, text_col), (x_offset_matrix, current_y))
                          current_y += line_h_matrix_elem # Add space for the ellipsis line
                     break # Stop drawing rows

                # Draw row label
                row_label_text_draw = unknown_labels[r_idx_draw] if title_str.startswith("Original Y") else f"r{r_idx_draw}"
                r_label_s = self.font_tiny.render(row_label_text_draw + ":", True, header_col);
                # Position row label right-aligned before the matrix content
                popup_surf.blit(r_label_s, r_label_s.get_rect(right=x_offset_matrix - int(3*SIZE_SCALE) , top=current_y)); # Add slight gap

                temp_x_elems = x_offset_matrix # Start X for elements in this row

                # Draw elements in the current row
                for c_idx_draw_elem in range(current_max_cols):
                    elem_s = self.font_tiny.render(str_elements[r_idx_draw][c_idx_draw_elem], True, text_col);
                    # Position element text right-aligned within its column width
                    elem_rect = elem_s.get_rect(right=temp_x_elems + col_widths[c_idx_draw_elem], centery=current_y + line_h_matrix_elem//2-1);
                    popup_surf.blit(elem_s, elem_rect);
                    temp_x_elems += col_widths[c_idx_draw_elem] + int(3*SIZE_SCALE); # Move x to the next element

                # Indicate if columns are truncated in this row
                if dim > current_max_cols:
                    popup_surf.blit(self.font_tiny.render("...", True, text_col), (temp_x_elems, current_y));

                current_y += line_h_matrix_elem; # Move Y down for the next row

            # Indicate if rows are truncated
            if dim > current_max_rows and current_y + line_h_matrix_elem <= popup_h - padding - self.font_small.get_linesize():
                 popup_surf.blit(self.font_tiny.render("... (more rows)", True, text_col), (x_offset_matrix, current_y));
                 current_y += line_h_matrix_elem;

            current_y += int(8*SIZE_SCALE); # Add spacing after the matrix block


        # Draw close instruction at the bottom of the popup
        close_msg = self.font_small.render("Press U or Esc to close LU View", True, (200,200,200)); # Light gray text, updated key
        popup_surf.blit(close_msg, close_msg.get_rect(centerx=popup_w//2, bottom=popup_h-padding//2)); # Position centered at the bottom

        # Blit the popup surface onto the main screen
        self.scr.blit(popup_surf, (popup_x, popup_y));


    # ───────── Live graph popup ──────────
    def _maybe_sample_live_waveforms(self):
        if self.sim_time < self._next_sample_t:
            return
        self._next_sample_t = self.sim_time + 1.0 / self.graph_samples_hz

        omega = 2*math.pi*self.freq
        max_src_v = max((c.val for c in self.circ.comps
                         if c.ct == CType.VAC), default=1.0)

        for c in self.circ.comps:
            v_i = self._instant_v_i(c, omega, self.sim_time)
            if v_i is None:
                continue
            v, i = v_i
            # normalise to max source magnitude (avoid ÷0)
            nv = v / max_src_v
            ni = i / max_src_v
            self.history[(c, 'V')].append((self.sim_time, nv))
            self.history[(c, 'I')].append((self.sim_time, ni))

    def _instant_v_i(self, c:Comp, omega:float, t:float):
        """Return instantaneous (V, I) for component *c*
           or None if we cannot compute it yet."""
        if self.circ.solution is None or self.circ.solution.size == 0:
            return None
        # === voltage exactly like overlay ===
        net_map = self.circ.nmap
        node_volt = self.circ.solution[:self.circ.node_cnt]

        def pin_v(p):
            return (0j if p.net == 0 else node_volt[net_map[p.net]]
                   ) if p.net in net_map else 0j

        # Instantaneous voltage at *first* pin
        vp = pin_v(c.pins[0])
        mag, ph = cmath.polar(vp)
        v_inst_1 = mag * math.cos(omega*t + ph)

        # Same for the second pin (if any) to get ΔV
        if len(c.pins) > 1:
            vn = pin_v(c.pins[1])
            mag2, ph2 = cmath.polar(vn)
            v_inst_2 = mag2 * math.cos(omega*t + ph2)
            v_inst = v_inst_1 - v_inst_2
        else:
            v_inst = v_inst_1

        # === current exactly like overlay ===
        cur_phasor = None
        if c.ct in (CType.R,CType.C,CType.L) and len(c.pins)==2:
            # reuse overlay logic
            p1,p2 = c.pins
            v1, v2 = vp, vn
            Z = {CType.R: lambda: max(c.val,1e-12),
                 CType.C: lambda: (1e30 if omega==0 or c.val==0 else 1/(1j*omega*c.val)),
                 CType.L: lambda: (R_L_DC_SHORT if omega==0 or c.val==0 else 1j*omega*c.val)}[c.ct]()
            cur_phasor = (v1 - v2) / Z
        elif c.ct in (CType.VDC, CType.VAC):
            idx = self.circ.vsrc_idx.get(c)
            if idx is not None:
                cur_phasor = self.circ.solution[self.circ.node_cnt+idx]

        if cur_phasor is None:
            return v_inst, 0.0

        magI, phI = cmath.polar(cur_phasor)
        i_inst = magI * math.cos(omega*t + phI)
        return v_inst, i_inst

    def draw_live_graph_popup(self):
        if not self.history:
            return

        popup_w, popup_h = int(WIN_W*0.9), int(WIN_H*0.6)
        px, py = (WIN_W-popup_w)//2, (WIN_H-popup_h)//2
        surf = pg.Surface((popup_w, popup_h), pg.SRCALPHA)
        pg.draw.rect(surf,(30,35,40,240),surf.get_rect(),border_radius=8)
        pg.draw.rect(surf,(120,140,160,230),surf.get_rect(),2,border_radius=8)

        title = self.font.render("Live waveforms – last 5 s (normalised)",True,(220,220,255))
        surf.blit(title,title.get_rect(midtop=(popup_w//2,8)))

        # graph area
        pad = int(40*SIZE_SCALE)
        gx0, gy0 = pad,  title.get_rect().bottom + pad//2
        gx1, gy1 = popup_w-pad, popup_h-pad
        w,  h   = gx1-gx0, gy1-gy0

        # axes
        pg.draw.line(surf,(90,90,90),(gx0,gy1),(gx1,gy1))
        pg.draw.line(surf,(90,90,90),(gx0,gy0),(gx0,gy1))
        # zero-line
        pg.draw.line(surf,(70,70,70),(gx0,(gy0+gy1)//2),(gx1,(gy0+gy1)//2),1)

        # time scale
        t_now = self.sim_time
        t_min = max(0.0, t_now-self.graph_history_s)
        def x_of(t): return gx0 + (t-t_min)/(self.graph_history_s)*w
        def y_of(v): return (gy0+gy1)//2 - v* (h/2-4)

        colours = {}
        colour_cycle = [(255,100,100),(100,200,255),(255,220,120),
                        (180,255,180),(255,180,255),(160,140,255)]
        for k in self.history.keys():
            if k[0] not in colours:
                colours[k[0]] = colour_cycle[len(colours)%len(colour_cycle)]

        # draw each trace
        for (comp, kind), dq in self.history.items():
            col = colours[comp]
            if kind == 'I':
                col = tuple(min(255,x+60) for x in col)  # lighter for current
            pts = [(x_of(t), y_of(v)) for t,v in dq if t>=t_min]
            if len(pts) > 1:
                pg.draw.lines(surf,col,False,pts,1)

        # legend
        lx, ly = gx0, gy0- int(22*SIZE_SCALE)
        for comp, col in colours.items():
            name_s = self.font_tiny.render(f"{comp.name}  (V)",True,col)
            surf.blit(name_s,(lx,ly)); ly-=name_s.get_height()+2
            name_s = self.font_tiny.render(f"{comp.name}  (I)",True,
                                           tuple(min(255,x+60) for x in col))
            surf.blit(name_s,(lx,ly)); ly-=name_s.get_height()+4

        # close hint
        hint = self.font_small.render("Press G or Esc to close",True,(200,200,200))
        surf.blit(hint,hint.get_rect(midbottom=(popup_w//2,popup_h-6)))

        self.scr.blit(surf,(px,py))


# ───────── main ─────────
# Ensure the main execution block runs the App
if __name__=="__main__":
    try:
        App().run()
    except KeyboardInterrupt:
        # Handle Ctrl+C gracefully in terminal
        print("\nCircuit Canvas stopped.")
        pass
    finally:
        # Ensure Pygame is quit properly
        pg.quit()
        # sys.exit() # sys.exit() might be needed depending on environment