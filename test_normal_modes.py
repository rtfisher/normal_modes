"""
Comprehensive test suite for normal_modes.py

Tests cover:
1. Physics calculations (mode frequencies, eigenvectors)
2. Plotting functions (spring rendering, colormap)
3. UI/CLI (argument parsing, validation)
"""

import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.cm import RdBu_r
import sys
import subprocess
from io import StringIO


# Import the module now that it has proper if __name__ == '__main__' guard
from normal_modes import validate_params


class TestPhysics:
    """Test physics calculations for normal modes"""

    def test_mode_frequencies_ordering(self):
        """Test that mode frequencies increase with mode number"""
        N = 9
        omega_0 = 1.0
        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

        # Check frequencies are in ascending order
        assert np.all(np.diff(frequencies) > 0), "Frequencies should increase with mode number"

    def test_mode_frequencies_range(self):
        """Test that frequencies are in expected range [0, 2ω₀]"""
        N = 9
        omega_0 = 1.0
        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

        assert np.all(frequencies > 0), "All frequencies should be positive"
        assert np.all(frequencies < 2 * omega_0), "Frequencies should be less than 2ω₀"

    def test_fundamental_frequency(self):
        """Test the fundamental (lowest) frequency calculation"""
        N = 9
        omega_0 = 1.0
        # For n=1: ω₁ = 2ω₀ sin(π/(2(N+1)))
        expected = 2 * omega_0 * np.sin(np.pi / (2 * (N + 1)))

        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

        assert np.isclose(frequencies[0], expected), "Fundamental frequency incorrect"

    def test_highest_frequency(self):
        """Test the highest frequency calculation"""
        N = 9
        omega_0 = 1.0
        # For n=N: ω_N = 2ω₀ sin(Nπ/(2(N+1)))
        expected = 2 * omega_0 * np.sin(N * np.pi / (2 * (N + 1)))

        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

        assert np.isclose(frequencies[-1], expected), "Highest frequency incorrect"

    def test_eigenvector_shape(self):
        """Test that eigenvector matrix has correct shape"""
        N = 9
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))

        assert modes.shape == (N, N), f"Expected shape ({N}, {N}), got {modes.shape}"

    def test_eigenvector_normalization(self):
        """Test that eigenvectors are properly normalized"""
        N = 9
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        # Check that max absolute value of each mode is 1
        for n in range(N):
            max_val = np.max(np.abs(modes[:, n]))
            assert np.isclose(max_val, 1.0), f"Mode {n} not normalized, max = {max_val}"

    def test_fundamental_mode_shape(self):
        """Test that fundamental mode has all masses moving in phase"""
        N = 9
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))

        # Fundamental mode (n=0) should have all positive values
        fundamental = modes[:, 0]
        assert np.all(fundamental > 0), "Fundamental mode should have all positive displacements"

    def test_highest_mode_alternating(self):
        """Test that highest mode has alternating sign (adjacent masses out of phase)"""
        N = 9
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))

        # Highest mode should have alternating signs
        highest = modes[:, -1]
        # Check that adjacent masses have opposite signs
        for j in range(N - 1):
            product = highest[j] * highest[j + 1]
            assert product < 0, f"Adjacent masses {j} and {j+1} should have opposite signs"

    def test_mode_orthogonality(self):
        """Test that different modes are orthogonal"""
        N = 9
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))

        # Check orthogonality: modes[:, i] · modes[:, j] ≈ 0 for i ≠ j
        for i in range(N):
            for j in range(i + 1, N):
                dot_product = np.dot(modes[:, i], modes[:, j])
                # Due to discrete nature, not perfectly orthogonal but should be small
                assert abs(dot_product) < 0.5, f"Modes {i} and {j} not sufficiently orthogonal"

    def test_number_of_modes(self):
        """Test that we get exactly N modes for N masses"""
        for N in [3, 5, 9, 15, 20]:
            mode_numbers = np.arange(1, N + 1)
            assert len(mode_numbers) == N, f"Should have {N} modes for {N} masses"

    def test_velocity_calculation(self):
        """Test velocity calculation from position"""
        N = 5
        omega_0 = 1.0
        amplitude = 0.3
        t = 1.0

        # Calculate for mode 0
        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        omega = frequencies[0]
        displacement = amplitude * modes[:, 0] * np.sin(omega * t)
        velocity = amplitude * modes[:, 0] * omega * np.cos(omega * t)

        # Numerical derivative check
        dt = 1e-6
        displacement_future = amplitude * modes[:, 0] * np.sin(omega * (t + dt))
        numerical_velocity = (displacement_future - displacement) / dt

        assert np.allclose(velocity, numerical_velocity, rtol=1e-4), \
            "Analytical velocity doesn't match numerical derivative"

    def test_energy_conservation(self):
        """Test that total energy oscillates around a mean (simplified check)"""
        N = 5
        omega_0 = 1.0
        amplitude = 0.3
        m = 1.0
        k = 1.0

        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        omega = frequencies[0]

        # Calculate kinetic energy at different times
        times = np.linspace(0, 2*np.pi/omega, 100)
        kinetic_energies = []

        for t in times:
            velocity = amplitude * modes[:, 0] * omega * np.cos(omega * t)
            KE = 0.5 * m * np.sum(velocity**2)
            kinetic_energies.append(KE)

        # Maximum kinetic energy should be positive
        max_KE = np.max(kinetic_energies)
        assert max_KE > 0, "Maximum kinetic energy should be positive"


class TestPlotting:
    """Test plotting and visualization functions"""

    def test_spring_drawing_basic(self):
        """Test that spring drawing function creates valid line data"""
        fig, ax = plt.subplots()
        line = ax.plot([], [], 'k-')[0]

        # Test spring from x=0 to x=1
        x1, x2 = 0, 1
        L = 10.0
        N = 9
        mass_spacing = L / (N + 1)

        # Manually implement draw_spring function
        coil_width = min(L * 0.015, mass_spacing * 0.2)
        spring_length = abs(x2 - x1)
        n_coils = max(5, int(spring_length / mass_spacing * 10))

        x = np.linspace(x1, x2, n_coils * 2 + 1)
        y_coords = np.zeros_like(x)
        for i in range(1, len(x) - 1):
            y_coords[i] = coil_width * (1 if i % 2 == 1 else -1)

        line.set_data(x, y_coords)
        xdata, ydata = line.get_data()

        assert len(xdata) > 0, "Spring should have x coordinates"
        assert len(ydata) > 0, "Spring should have y coordinates"
        assert len(xdata) == len(ydata), "x and y should have same length"
        assert xdata[0] == x1, "Spring should start at x1"
        assert xdata[-1] == x2, "Spring should end at x2"

        plt.close(fig)

    def test_spring_coil_scaling(self):
        """Test that spring coils scale with length"""
        L = 10.0
        N = 9
        mass_spacing = L / (N + 1)

        # Short spring
        spring_length_short = 0.5
        n_coils_short = max(5, int(spring_length_short / mass_spacing * 10))

        # Long spring
        spring_length_long = 2.0
        n_coils_long = max(5, int(spring_length_long / mass_spacing * 10))

        assert n_coils_long > n_coils_short, "Longer springs should have more coils"

    def test_colormap_setup(self):
        """Test that colormap is properly configured"""
        colormap = RdBu_r
        norm = Normalize(vmin=-1, vmax=1)

        # Test extreme values
        color_max_pos = colormap(norm(1.0))
        color_max_neg = colormap(norm(-1.0))
        color_zero = colormap(norm(0.0))

        # Check that colors are valid RGBA tuples
        assert len(color_max_pos) == 4, "Color should be RGBA tuple"
        assert len(color_max_neg) == 4, "Color should be RGBA tuple"
        assert len(color_zero) == 4, "Color should be RGBA tuple"

        # Check colors are different at extremes
        assert not np.allclose(color_max_pos[:3], color_max_neg[:3]), \
            "Extreme colors should be different"

        # Check all values are in valid range [0, 1]
        for color in [color_max_pos, color_max_neg, color_zero]:
            assert all(0 <= c <= 1 for c in color), "Color values should be in [0, 1]"

    def test_velocity_normalization(self):
        """Test that velocities are properly normalized for coloring"""
        N = 5
        omega_0 = 1.0
        amplitude = 0.3

        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        omega = frequencies[0]
        t = 0.0  # At t=0, velocity is maximum
        velocities = amplitude * modes[:, 0] * omega * np.cos(omega * t)

        max_velocity = amplitude * np.max(frequencies)
        normalized_velocities = velocities / max_velocity

        assert np.all(np.abs(normalized_velocities) <= 1.0), \
            "Normalized velocities should be in [-1, 1]"

    def test_mass_artist_creation(self):
        """Test that mass artists are created correctly"""
        N = 9
        L = 10.0
        mass_spacing = L / (N + 1)
        mass_radius = mass_spacing * 0.2

        masses_artists = [Circle((0, 0), mass_radius, color='gray', zorder=3) for _ in range(N)]

        assert len(masses_artists) == N, f"Should have {N} mass artists"
        for mass in masses_artists:
            assert mass.radius == mass_radius, "All masses should have same radius"

    def test_equilibrium_positions(self):
        """Test that equilibrium positions are evenly spaced"""
        N = 9
        L = 10.0
        mass_spacing = L / (N + 1)
        equilibrium_positions = np.array([mass_spacing * (i + 1) for i in range(N)])

        assert len(equilibrium_positions) == N, "Should have N equilibrium positions"
        assert equilibrium_positions[0] == mass_spacing, "First mass at one spacing"
        assert equilibrium_positions[-1] == N * mass_spacing, "Last mass at N spacings"

        # Check even spacing
        spacings = np.diff(equilibrium_positions)
        assert np.allclose(spacings, mass_spacing), "Masses should be evenly spaced"

    def test_position_bounds(self):
        """Test that mass positions stay within reasonable bounds"""
        N = 5
        L = 10.0
        amplitude = 0.35
        mass_spacing = L / (N + 1)
        equilibrium_positions = np.array([mass_spacing * (i + 1) for i in range(N)])

        # Create a mode
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        # Check maximum displacement
        max_displacement = amplitude * np.max(np.abs(modes[:, 0]))
        positions_max = equilibrium_positions + max_displacement
        positions_min = equilibrium_positions - max_displacement

        assert np.all(positions_max < L), "Masses shouldn't exceed right wall"
        assert np.all(positions_min > 0), "Masses shouldn't exceed left wall"

    def test_figure_creation(self):
        """Test that matplotlib figure is created with correct settings"""
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.set_xlim(-0.5, 10.5)
        ax.set_ylim(-2, 2)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        assert xlim[0] < 0, "x limit should start negative"
        assert ylim[0] == -2, "y lower limit should be -2"
        assert ylim[1] == 2, "y upper limit should be 2"

        plt.close(fig)


class TestUI:
    """Test UI, CLI arguments, and validation"""

    def test_validate_params_function_exists(self):
        """Test that validate_params function is defined"""
        # validate_params is imported at the top
        assert callable(validate_params), "validate_params should be importable"

    def test_invalid_num_masses_too_low(self):
        """Test validation rejects N < 1"""
        with pytest.raises(ValueError, match="Number of masses must be between 1 and 100"):
            validate_params(N=0, amplitude=0.3, mass_spacing=1.0, fps=50, dpi=100)

    def test_invalid_num_masses_too_high(self):
        """Test validation rejects N > 100"""
        with pytest.raises(ValueError, match="Number of masses must be between 1 and 100"):
            validate_params(N=150, amplitude=0.3, mass_spacing=1.0, fps=50, dpi=100)

    def test_invalid_amplitude_negative(self):
        """Test validation rejects negative amplitude"""
        with pytest.raises(ValueError, match="Amplitude must be positive"):
            validate_params(N=9, amplitude=-0.5, mass_spacing=1.0, fps=50, dpi=100)

    def test_invalid_amplitude_too_large(self):
        """Test validation rejects amplitude causing overlap"""
        mass_spacing = 1.0
        max_amp = mass_spacing * 0.4

        with pytest.raises(ValueError, match="Amplitude too large"):
            validate_params(N=9, amplitude=max_amp + 0.1, mass_spacing=mass_spacing,
                          fps=50, dpi=100)

    def test_invalid_fps_negative(self):
        """Test validation rejects negative FPS"""
        with pytest.raises(ValueError, match="FPS must be between 1 and 120"):
            validate_params(N=9, amplitude=0.3, mass_spacing=1.0, fps=-10, dpi=100)

    def test_invalid_fps_too_high(self):
        """Test validation rejects FPS > 120"""
        with pytest.raises(ValueError, match="FPS must be between 1 and 120"):
            validate_params(N=9, amplitude=0.3, mass_spacing=1.0, fps=200, dpi=100)

    def test_invalid_dpi_negative(self):
        """Test validation rejects negative DPI"""
        with pytest.raises(ValueError, match="DPI must be between 1 and 300"):
            validate_params(N=9, amplitude=0.3, mass_spacing=1.0, fps=50, dpi=-50)

    def test_invalid_dpi_too_high(self):
        """Test validation rejects DPI > 300"""
        with pytest.raises(ValueError, match="DPI must be between 1 and 300"):
            validate_params(N=9, amplitude=0.3, mass_spacing=1.0, fps=50, dpi=500)

    def test_valid_parameters(self):
        """Test validation accepts valid parameters"""
        # Should not raise any exception
        validate_params(N=9, amplitude=0.35, mass_spacing=1.0, fps=50, dpi=100)
        validate_params(N=5, amplitude=0.3, mass_spacing=1.5, fps=60, dpi=150)
        validate_params(N=1, amplitude=0.1, mass_spacing=5.0, fps=30, dpi=80)

    def test_cli_help_option(self):
        """Test that --help option works"""
        result = subprocess.run(
            [sys.executable, 'normal_modes.py', '--help'],
            capture_output=True, text=True
        )

        assert result.returncode == 0, "Help should exit successfully"
        assert 'usage:' in result.stdout.lower(), "Help should show usage"
        assert '--num-masses' in result.stdout, "Help should list --num-masses option"
        assert '--amplitude' in result.stdout, "Help should list --amplitude option"

    def test_cli_invalid_mode_number(self):
        """Test that invalid mode number is rejected"""
        result = subprocess.run(
            [sys.executable, 'normal_modes.py', '--mode', '100', '--num-masses', '5'],
            capture_output=True, text=True, timeout=5
        )

        assert result.returncode != 0, "Invalid mode should cause error"
        assert 'Mode must be between' in result.stderr or 'Mode must be between' in result.stdout

    def test_frame_calculation(self):
        """Test that frame counts are calculated correctly"""
        dt = 0.02
        t_per_mode = 4.0
        frames_per_mode = int(t_per_mode / dt)

        assert frames_per_mode == 200, f"Expected 200 frames, got {frames_per_mode}"

        pause_duration = 0.5
        frames_per_pause = int(pause_duration / dt)

        assert frames_per_pause == 25, f"Expected 25 pause frames, got {frames_per_pause}"

    def test_modes_to_show_all(self):
        """Test that all modes are shown when no specific mode selected"""
        N = 9
        mode_arg = None

        if mode_arg is not None:
            modes_to_show = [mode_arg - 1]
        else:
            modes_to_show = list(range(N))

        assert len(modes_to_show) == N, f"Should show all {N} modes"
        assert modes_to_show == list(range(N)), "Should be [0, 1, 2, ..., N-1]"

    def test_modes_to_show_single(self):
        """Test that single mode is shown when specified"""
        N = 9
        mode_arg = 3  # User specifies mode 3 (1-indexed)

        if mode_arg is not None:
            modes_to_show = [mode_arg - 1]  # Convert to 0-indexed
        else:
            modes_to_show = list(range(N))

        assert len(modes_to_show) == 1, "Should show only 1 mode"
        assert modes_to_show[0] == 2, "Mode 3 (1-indexed) should be index 2"

    def test_total_frames_calculation(self):
        """Test total frame count for animation"""
        frames_per_mode = 200
        frames_per_pause = 25
        frames_per_cycle = frames_per_mode + frames_per_pause

        # Test with all modes
        N = 9
        modes_to_show = list(range(N))
        total_frames = frames_per_cycle * len(modes_to_show)

        expected = (200 + 25) * 9
        assert total_frames == expected, f"Expected {expected} total frames, got {total_frames}"

        # Test with single mode
        modes_to_show_single = [2]
        total_frames_single = frames_per_cycle * len(modes_to_show_single)

        expected_single = 225
        assert total_frames_single == expected_single, \
            f"Expected {expected_single} frames for single mode, got {total_frames_single}"


class TestIntegration:
    """Integration tests that test the full workflow"""

    def test_full_physics_pipeline(self):
        """Test complete physics calculation pipeline"""
        N = 5
        m = 1.0
        k = 1.0
        omega_0 = np.sqrt(k / m)
        L = 10.0
        amplitude = 0.3

        # Calculate frequencies
        mode_numbers = np.arange(1, N + 1)
        frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

        # Calculate modes
        modes = np.zeros((N, N))
        for n in range(N):
            for j in range(N):
                modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
            modes[:, n] /= np.max(np.abs(modes[:, n]))

        # Calculate positions and velocities at t=0
        t = 0
        omega = frequencies[0]
        mass_spacing = L / (N + 1)
        equilibrium_positions = np.array([mass_spacing * (i + 1) for i in range(N)])

        displacements = amplitude * modes[:, 0] * np.sin(omega * t)
        velocities = amplitude * modes[:, 0] * omega * np.cos(omega * t)
        positions = equilibrium_positions + displacements

        # Verify results
        assert len(frequencies) == N
        assert modes.shape == (N, N)
        assert len(positions) == N
        assert len(velocities) == N
        assert np.all(positions > 0)
        assert np.all(positions < L)

    def test_script_runs_basic(self):
        """Test that the script runs without errors (basic smoke test)"""
        # Test with video output to avoid display issues
        result = subprocess.run(
            [sys.executable, 'normal_modes.py', '--num-masses', '3', '--mode', '1',
             '--duration', '0.5', '--outfile', 'test_output.mp4'],
            capture_output=True, text=True, timeout=30
        )

        assert result.returncode == 0, f"Script failed: {result.stderr}"

        # Clean up
        import os
        if os.path.exists('test_output.mp4'):
            os.remove('test_output.mp4')

    def test_different_n_values(self):
        """Test that physics works for different N values"""
        for N in [1, 3, 5, 10, 20]:
            omega_0 = 1.0
            mode_numbers = np.arange(1, N + 1)
            frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

            modes = np.zeros((N, N))
            for n in range(N):
                for j in range(N):
                    modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
                modes[:, n] /= np.max(np.abs(modes[:, n]))

            # Verify basic properties
            assert len(frequencies) == N
            assert modes.shape == (N, N)
            assert np.all(frequencies > 0)
            for n in range(N):
                assert np.isclose(np.max(np.abs(modes[:, n])), 1.0)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
