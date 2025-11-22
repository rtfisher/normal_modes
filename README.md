# Normal Modes Visualization

[![CI](https://github.com/rtfisher/normal_modes/workflows/CI/badge.svg)](https://github.com/rtfisher/normal_modes/actions)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An interactive physics simulation that visualizes the normal modes of a 1D coupled oscillator system. This educational tool demonstrates classical mechanics concepts by animating N masses connected by springs between two fixed walls.

## Features

- **Accurate Physics**: Analytical calculation of normal mode frequencies and eigenvectors for fixed-boundary conditions
- **Velocity-Based Visualization**: Color-coded masses using divergent colormap (red = moving right, blue = moving left)
- **Adaptive Rendering**: Spring coil count and width scale dynamically with compression/extension
- **Flexible Configuration**: Customize number of masses, amplitude, duration, and which modes to display
- **Video Export**: Save high-quality animations to MP4 format with configurable FPS and DPI
- **Progress Feedback**: Real-time rendering progress display during video export
- **Input Validation**: Comprehensive parameter validation to prevent invalid configurations

## Installation

### Prerequisites

- Python 3.8 or higher
- ffmpeg (required only for saving animations to video files)

### Install Python Dependencies

```bash
pip install -r requirements.txt
```

### Install FFmpeg (Optional)

**macOS:**
```bash
brew install ffmpeg
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install ffmpeg
```

**Windows:**
Download from [ffmpeg.org](https://ffmpeg.org/download.html)

## Quick Start

### Display Animation Interactively

```bash
python normal_modes.py
```

This will display all 9 normal modes (default) in sequence, with each mode playing for 4 seconds.

### Save Animation to File

```bash
python normal_modes.py --outfile normal_modes.mp4
```

### Show Specific Mode

```bash
# Display only the third normal mode
python normal_modes.py --mode 3
```

## Usage

```bash
python normal_modes.py [OPTIONS]
```

### Command-Line Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--outfile FILENAME` | str | None | Save animation to file (requires ffmpeg) |
| `--fps FPS` | int | 50 | Frames per second for video output |
| `--dpi DPI` | int | 100 | Resolution (dots per inch) for video |
| `--num-masses N` | int | 9 | Number of masses in the system (1-100) |
| `--amplitude A` | float | 0.35 | Oscillation amplitude |
| `--mode N` | int | None | Show only specific mode (1 to N) |
| `--duration D` | float | 4.0 | Time per mode in seconds |
| `--pause P` | float | 0.5 | Pause between mode transitions (seconds) |

### Examples

**Create high-quality video with 15 masses:**
```bash
python normal_modes.py --num-masses 15 --outfile demo.mp4 --fps 60 --dpi 150
```

**Quick preview of first three modes:**
```bash
python normal_modes.py --num-masses 9 --duration 2 --pause 0.2
```

**Study the fundamental mode:**
```bash
python normal_modes.py --mode 1 --duration 10 --amplitude 0.4
```

**Large amplitude oscillation (5 masses):**
```bash
python normal_modes.py --num-masses 5 --amplitude 0.6 --outfile big_amplitude.mp4
```

## Physics Model

### System Description

The simulation models N equal masses connected by equal springs with fixed boundaries:

```
[WALL]--spring--[mass₁]--spring--[mass₂]--...--[massₙ]--spring--[WALL]
```

- **Masses**: All equal (m = 1.0)
- **Springs**: All have equal spring constant (k = 1.0)
- **Boundary Conditions**: Fixed at both ends (Dirichlet)

### Normal Mode Frequencies

For N masses with fixed boundaries, the angular frequencies are:

```
ωₙ = 2ω₀ sin(nπ / (2(N+1)))    for n = 1, 2, ..., N
```

where `ω₀ = √(k/m)` is the natural frequency of a single oscillator.

### Mode Shapes (Eigenvectors)

The displacement of mass j in mode n is:

```
uⱼ⁽ⁿ⁾ = sin(n·j·π / (N+1))    for j = 1, 2, ..., N
```

Eigenvectors are normalized such that the maximum displacement is 1.

### Time Evolution

For a pure mode n, the position of mass j at time t is:

```
xⱼ(t) = x₀ⱼ + A·uⱼ⁽ⁿ⁾·sin(ωₙt)
```

where:
- `x₀ⱼ` is the equilibrium position
- `A` is the amplitude
- `uⱼ⁽ⁿ⁾` is the mode shape
- `ωₙ` is the mode frequency

The velocity is:

```
vⱼ(t) = A·ωₙ·uⱼ⁽ⁿ⁾·cos(ωₙt)
```

## Visualization

### Color Coding

Masses are colored based on their instantaneous velocity using the RdBu_r (Red-Blue reversed) divergent colormap:

- **Red**: Moving to the right (positive velocity)
- **White**: At rest (zero velocity, turning points)
- **Blue**: Moving to the left (negative velocity)

This helps identify:
- **Nodes**: Masses that remain white (zero displacement and velocity)
- **Phase relationships**: Masses with the same color are moving in the same direction
- **Maximum velocity**: Brightest colors indicate highest speed

### Spring Rendering

Springs are drawn as zigzag patterns with:
- **Adaptive coil count**: Number of coils scales with spring length
- **Dynamic width**: Coil width scales with system size (L-aware)
- **Realistic compression**: Visually shows when springs are compressed or extended

## Code Structure

The simulation consists of three main components:

### 1. Physical Setup (lines 62-90)
- Define system parameters (N, m, k, L)
- Calculate analytical mode frequencies
- Construct normalized eigenvector matrix
- Validate input parameters

### 2. Visualization Setup (lines 92-133)
- Create matplotlib figure and axes
- Initialize Circle artists for masses
- Create spring line artists
- Setup colormap for velocity-based coloring
- Add wall rectangles

### 3. Animation Loop (lines 159-228)
- Update mass positions based on mode and time
- Calculate velocities for color coding
- Update spring geometry dynamically
- Handle mode transitions with pauses
- Render frames for display or video export

## Development

### Running Tests

```bash
# Run all tests
pytest test_normal_modes.py -v

# Run specific test category
pytest test_normal_modes.py::TestPhysics -v
pytest test_normal_modes.py::TestPlotting -v
pytest test_normal_modes.py::TestUI -v

# Run with coverage
pytest test_normal_modes.py --cov=normal_modes --cov-report=html
```

### Continuous Integration

The project uses GitHub Actions for automated testing. On each push or pull request, the CI workflow:

1. Tests on multiple Python versions (3.8, 3.9, 3.10, 3.11)
2. Runs the complete test suite
3. Generates coverage reports
4. Validates code quality

See `.github/workflows/ci.yml` for details.

## Educational Use

This simulation is designed for junior-level physics major classical mechanics classes and demonstrates:

- **Normal mode decomposition** of coupled systems
- **Eigenvalue problems** in classical mechanics
- **Orthogonality** of normal modes
- **Standing waves** in discrete systems
- **Frequency spectrum** of coupled oscillators

### Key Concepts to Explore

1. **Fundamental Mode (n=1)**: All masses move in phase, lowest frequency
2. **Highest Mode (n=N)**: Adjacent masses move out of phase, highest frequency
3. **Nodes**: Higher modes have stationary points (masses that don't move)
4. **Frequency Ordering**: Mode frequencies increase with mode number
5. **Symmetry**: Mode shapes are symmetric about the center

## Troubleshooting

### "ffmpeg not found" Error

If you get this error when trying to save a video:

1. Install ffmpeg (see Installation section)
2. Verify installation: `ffmpeg -version`
3. Make sure ffmpeg is in your system PATH

### "Amplitude too large" Error

The validation prevents amplitudes that would cause mass overlap:

```
Maximum amplitude = 0.4 × (L / (N+1))
```

For N=9 masses with L=10, max amplitude ≈ 0.4

### Interactive Display Not Working

If `plt.show()` doesn't display anything:

1. Check that you have a GUI backend for matplotlib
2. Try setting a backend explicitly: `export MPLBACKEND=TkAgg`
3. Alternatively, always use `--outfile` to save to video

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Make your changes and add tests
4. Run the test suite (`pytest`)
5. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
6. Push to the branch (`git push origin feature/AmazingFeature`)
7. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this simulation in your research or teaching, please cite:

```bibtex
@software{normal_modes_2025,
  title = {Normal Modes Visualization},
  author = {Robert Fisher},
  year = {2025},
  url = {https://github.com/rtfisher/normal-modes}
}
```

## Acknowledgments

- Physics model based on classical mechanics textbook treatment of coupled oscillators
- Developed for junior-level phyiscs major classical mechanics class
- Visualization inspired by interactive physics demonstrations

## References

1. Thornton, S. T., & Marion, J. B. (2003). *Classical Dynamics of Particles and Systems*. Brooks/Cole.
2. Taylor, J. R. (2005). *Classical Mechanics*. University Science Books.
3. Goldstein, H., Poole, C., & Safko, J. (2002). *Classical Mechanics*. Addison-Wesley.
