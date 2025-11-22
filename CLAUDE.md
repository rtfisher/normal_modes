# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a physics education simulation that visualizes the normal modes of a 1D coupled oscillator system. The system consists of N masses connected by springs between two fixed walls, demonstrating classical mechanics concepts for PHY313.

## Running the Simulation

```bash
# Display animation interactively
python normal_modes.py

# Save animation to file (requires ffmpeg)
python normal_modes.py --outfile normal_modes.mp4

# Customize output settings
python normal_modes.py --outfile output.mp4 --fps 30 --dpi 150
```

**Command line options:**
- `--outfile FILENAME`: Save animation to file instead of displaying (requires ffmpeg)
- `--fps FPS`: Frames per second for output video (default: calculated from dt as 1/dt = 50)
- `--dpi DPI`: Resolution for output video (default: 100)

The script creates an animated matplotlib visualization showing each normal mode sequentially, with each mode displayed for 4 seconds before transitioning to the next.

## Physics Model

The system simulates N equal masses (default: 9) connected by equal springs with fixed boundary conditions:

- **Mode frequencies**: `ω_n = 2ω₀ sin(nπ/(2(N+1)))` for n = 1, 2, ..., N
- **Mode shapes**: The displacement of mass j in mode n is proportional to `sin(n·j·π/(N+1))`
- Eigenvectors are normalized to have maximum displacement of 1

## Key Parameters (normal_modes.py:7-21)

- `N`: Number of masses (line 8)
- `m`, `k`: Mass and spring constant (lines 9-10)
- `L`: Total system length (line 14)
- `dt`, `t_per_mode`: Animation time step and duration per mode (lines 19-20)
- `amplitude`: Oscillation amplitude (line 39)

## Code Structure

The simulation has three main sections:

1. **Physical setup** (lines 7-36): Defines system parameters, calculates mode frequencies analytically, and constructs normalized eigenvector matrix
2. **Visualization setup** (lines 41-77): Creates matplotlib figure with Circle artists for masses, line artists for springs, and Rectangle patches for walls
3. **Animation loop** (lines 98-134): Updates mass positions and spring geometry each frame using `displacements = amplitude * modes[:, current_mode] * sin(ω*t)`

The `draw_spring()` function (lines 79-88) renders springs as zigzag patterns that compress/extend as masses move.

## Modifying the Simulation

- To change the number of masses, modify `N` at line 8
- To adjust animation speed, change `dt` (line 19) or `interval` parameter in FuncAnimation (line 139)
- Mass and spring visual sizing automatically scales with `mass_spacing` (lines 50-51)
