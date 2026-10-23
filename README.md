# Molecular Dynamics: Lennard-Jones Potential & Velocity Verlet

A 2D molecular dynamics simulation of \(N\) interacting particles in a region with periodic boundary conditions. The interactions between particles are described by the Lennard-Jones potential, and the equations of motion are numerically integrated using the Velocity Verlet algorithm.

## 1. Project Structure

- `src/` — C source code (`main_verle.c`, `mt19937-64.c`, `mt64.h`)
- `scripts/` — Helper scripts (`animate.py` for visualization)
- `run_all.sh` — Main script for compilation and execution
- `lab1.md` — Original laboratory report text

## 2. System Requirements

The program runs on Linux (Ubuntu, Debian, Fedora, Arch). The following are required:
- **GCC** compiler
- **Math library** (linked with `-lm`)
- **Python 3** with `numpy` and `matplotlib`

Install Python packages:
```bash
pip3 install numpy matplotlib
```

## 3. Getting Started

**Preparation:**
Before the first run, make the script executable:
```bash
chmod +x run_all.sh
```

**Execution:**
Run the program:
```bash
./run_all.sh
```

After launching, a menu will appear to select the initial state:
```text
1 - GAS (random positions)
2 - GAS (regular lattice)
3 - LIQUID (small)
4 - LIQUID (large)
5 - SOLID (square lattice)
6 - SOLID (triangular lattice)
```
Enter a number (1-6) and press `Enter`.

## 4. What Happens After Selection

The script automatically performs the following:
1. Compiles the C code with the required macro (e.g., `-DSTATE_GAS_REGULAR`).
2. Runs the simulation.
3. Prints the kinetic, potential, and total energy at each saved step to the console.
4. Opens a window with particle animation (`animate.py`).
5. Waits for `Enter`, then closes the window and terminates.

## 5. Simulation Parameters

| Mode | State | N | L | Initial Velocity | Lattice |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | Gas (random) | 10 | 30.0 | 0.15 | None |
| 2 | Gas (regular) | 10 | 12.0 | 0.15 | Square |
| 3 | Liquid (small) | 36 | 12.0 | 0.06 | Square |
| 4 | Liquid (large) | 100 | 20.0 | 0.06 | Square |
| 5 | Solid (square) | 25 | 12.0 | 0.0008 | Square |
| 6 | Solid (triangular) | 25 | 12.0 | 0.0008 | Triangular |

**General Parameters:**
- Time step: `h = 0.0005`
- Total steps: `STEPS = 10000`
- Data saving interval: `SAVE_EVERY = 100` steps
- Potential cutoff radius: `rcut = 3.5 * sigma`

## 6. Program Output

The console outputs lines in the following format:
```text
Step X: Ek = value, Ep = value, E = value
```
- **Ek** — Kinetic energy (always positive)
- **Ep** — Potential energy (negative in liquid and solid states)
- **E = Ek + Ep** — Total energy (must be conserved)

## 7. Animation

The animation window displays particle motion:
- **Blue dots** represent particles.
- Field boundaries are fixed based on the first frame.
- The window title displays the current step number and the number of particles.
- Frame delay is **50 milliseconds** (20 FPS).

## 8. Algorithm & Physics

### Velocity Verlet Algorithm
1. Compute forces based on current coordinates: `a(t) = F(t) / m`
2. Update coordinates: `r(t+h) = r(t) + v(t)·h + 0.5·a(t)·h²`
3. Compute forces based on new coordinates: `a(t+h) = F(t+h) / m`
4. Update velocities: `v(t+h) = v(t) + 0.5·(a(t) + a(t+h))·h`

*Advantage:* Good energy conservation with a sufficiently small time step `h`.

### Periodic Boundary Conditions (PBC)
Implemented via the `min_image(dx)` function, which maps the coordinate difference between two particles into the range `[-L/2, L/2]`. This emulates an infinite system. 
- If `dx > L/2` → `dx = dx - L`
- If `dx < -L/2` → `dx = dx + L`
*(Same for dy)*. Coordinates may go outside `[0, L]`, but distances are always calculated correctly.

## 9. Correctness Criteria

The simulation is considered correct if:
- Total energy `E` remains approximately constant (drift ≤ 0.03 over 10,000 steps).
- Kinetic and potential energies oscillate in anti-phase (when `Ek` rises, `Ep` falls).
- No sudden energy jumps occur.
- The animation shows smooth particle motion.

## 10. Troubleshooting

| Problem | Solution |
| :--- | :--- |
| Compilation error | Ensure `main_verle.c`, `mt19937-64.c`, `mt64.h` are in `src/`, and GCC is installed. |
| `No module named numpy/matplotlib` | Run `pip3 install numpy matplotlib`. |
| Animation does not open | Ensure CSV files exist in `data/` and `animate.py` is in `scripts/`. |
| Animation window does not close | Press `Ctrl+C` in the terminal to abort the script, then close the window manually. |

## 11. Windows Compatibility

The program is written for Linux because it uses:
- Bash scripts (`run_all.sh`)
- POSIX paths (`/`)
- System calls (`system()`, `mkdir`)

To run it on Windows, you must use **WSL** (Windows Subsystem for Linux).

---

<details>
<summary><b>Click to expand: Theoretical Background (P.S.)</b></summary>

### 1. Lennard-Jones Potential
\[
U(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]
\]
- \(r\) — distance between particles
- \(\varepsilon = 1\) — depth of the potential well
- \(\sigma = 1\) — distance at which \(U(r) = 0\)
- \(r_{min} = 1.12\sigma\) — distance where energy is minimal
- \(U(r_{min}) = -\varepsilon = -1\)
- As \(r \to 0\): repulsion (\(U \to +\infty\))
- As \(r \to \infty\): attraction (\(U \to 0\))

### 2. Kinetic Energy
\[
E_k = \sum \frac{1}{2} m (v_x^2 + v_y^2)
\]
- \(m = 1\) — particle mass
- Summed over all \(N\) particles

### 3. Total Energy
\[
E = E_k + E_p
\]
\[
E_p = \sum U(r_{ij})
\]
In a closed system, \(E\) must be conserved.

### 4. Why Total Energy Changes Slightly
Due to numerical errors:
- Algorithm error is \(O(h^2)\).
- With \(h = 0.0005\), the error is ~0.03 over 10,000 steps.
- If there is no long-term drift, the simulation is correct.

### 5. Three States of Matter
- **Gas:** \(N = 10\), \(L = 30\), velocity \(0.15\), \(E \approx 0\). Particles fly chaotically and fill the volume.
- **Liquid:** \(N = 36\), \(L = 12\), velocity \(0.06\), \(E \approx -3\) to \(-5\). Particles move but stay together.
- **Solid:** \(N = 25\), \(L = 12\), velocity \(0.0008\), \(E \approx -20\) to \(-30\). Particles oscillate around lattice nodes.

### 6. Square vs. Triangular Lattices
- **Square lattice:** Each particle has 4 neighbors.
- **Triangular lattice:** Each particle has 6 neighbors (denser packing).
- The triangular lattice is energetically more favorable (lower \(E_p\)).

</details>
