import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, Rectangle
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import RdBu_r
import argparse
import sys


def validate_params(N, amplitude, mass_spacing, fps, dpi):
    """Validate input parameters"""
    if N < 1 or N > 100:
        raise ValueError("Number of masses must be between 1 and 100")
    if amplitude <= 0:
        raise ValueError("Amplitude must be positive")
    # Allow up to 80% of spacing to account for mass radius (20% each side = 40% total)
    max_amplitude = mass_spacing * 0.4
    if amplitude > max_amplitude:
        raise ValueError(f"Amplitude too large ({amplitude:.2f}) - masses may overlap. Max recommended: {max_amplitude:.2f}")
    if fps is not None and (fps <= 0 or fps > 120):
        raise ValueError("FPS must be between 1 and 120")
    if dpi <= 0 or dpi > 300:
        raise ValueError("DPI must be between 1 and 300")


def main():
    """Main function to run the simulation"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Visualize normal modes of coupled oscillators')
    parser.add_argument('--outfile', type=str, default=None,
                        help='Output filename for animation (e.g., normal_modes.mp4). Requires ffmpeg.')
    parser.add_argument('--fps', type=int, default=None,
                        help='Frames per second for output video (default: 1/dt)')
    parser.add_argument('--dpi', type=int, default=100,
                        help='DPI (resolution) for output video (default: 100)')
    parser.add_argument('--num-masses', type=int, default=9,
                        help='Number of masses in the system (default: 9)')
    parser.add_argument('--amplitude', type=float, default=0.35,
                        help='Oscillation amplitude (default: 0.35)')
    parser.add_argument('--mode', type=int, default=None,
                        help='Show only specific mode (1 to N, default: show all)')
    parser.add_argument('--duration', type=float, default=4.0,
                        help='Time per mode in seconds (default: 4.0)')
    parser.add_argument('--pause', type=float, default=0.5,
                        help='Pause duration between modes in seconds (default: 0.5)')
    args = parser.parse_args()

    # Parameters
    N = args.num_masses
    m = 1.0  # Mass (all equal)
    k = 1.0  # Spring constant (all equal)
    omega_0 = np.sqrt(k / m)  # Natural frequency

    # System geometry
    L = 10.0  # Total length
    mass_spacing = L / (N + 1)
    equilibrium_positions = np.array([mass_spacing * (i + 1) for i in range(N)])

    # Time parameters
    dt = 0.02
    t_per_mode = args.duration
    pause_duration = args.pause
    frames_per_mode = int(t_per_mode / dt)
    frames_per_pause = int(pause_duration / dt)

    # Amplitude of oscillation
    amplitude = args.amplitude

    # Validate parameters
    fps_to_validate = args.fps if args.fps is not None else int(1/dt)
    validate_params(N, amplitude, mass_spacing, fps_to_validate, args.dpi)

    # Validate and handle mode selection
    if args.mode is not None:
        if args.mode < 1 or args.mode > N:
            raise ValueError(f"Mode must be between 1 and {N}")
        modes_to_show = [args.mode - 1]  # Convert to 0-indexed
    else:
        modes_to_show = list(range(N))  # Show all modes

    # Calculate normal mode frequencies
    # For N masses with fixed boundaries, the frequencies are:
    # omega_n = 2 * omega_0 * sin(n * pi / (2 * (N + 1))) for n = 1, 2, ..., N
    mode_numbers = np.arange(1, N + 1)
    frequencies = 2 * omega_0 * np.sin(mode_numbers * np.pi / (2 * (N + 1)))

    # Normal mode eigenvectors (displacement patterns)
    # For mode n, the displacement of mass j is proportional to sin(n * j * pi / (N + 1))
    modes = np.zeros((N, N))
    for n in range(N):
        for j in range(N):
            modes[j, n] = np.sin((n + 1) * (j + 1) * np.pi / (N + 1))
        # Normalize
        modes[:, n] /= np.max(np.abs(modes[:, n]))

    # Setup figure
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_xlim(-0.5, L + 0.5)
    ax.set_ylim(-2, 2)
    ax.set_aspect('equal')
    ax.axis('off')

    # Setup colormap for velocity-based coloring
    # RdBu_r: Red for positive velocity (moving right), Blue for negative (moving left)
    colormap = RdBu_r
    norm = Normalize(vmin=-1, vmax=1)  # Will normalize velocities to [-1, 1]

    # Create artists for masses and springs
    # Scale mass radius based on spacing - ensure no overlap even with displacement
    mass_radius = mass_spacing * 0.2  # 20% of spacing
    masses_artists = [Circle((0, 0), mass_radius, color='gray', zorder=3) for _ in range(N)]
    for mass in masses_artists:
        ax.add_patch(mass)

    # Springs (lines connecting masses and to walls)
    spring_lines = []
    # Spring from left wall to first mass
    spring_lines.append(ax.plot([], [], 'k-', linewidth=2, zorder=1)[0])
    # Springs between masses
    for _ in range(N - 1):
        spring_lines.append(ax.plot([], [], 'k-', linewidth=2, zorder=1)[0])
    # Spring from last mass to right wall
    spring_lines.append(ax.plot([], [], 'k-', linewidth=2, zorder=1)[0])

    # Walls
    wall_width = 0.2
    wall_height = 1.5
    left_wall = Rectangle((-wall_width, -wall_height/2), wall_width, wall_height,
                          color='gray', zorder=2)
    right_wall = Rectangle((L, -wall_height/2), wall_width, wall_height,
                           color='gray', zorder=2)
    ax.add_patch(left_wall)
    ax.add_patch(right_wall)

    # Title - position dynamically based on ylim
    title_y = ax.get_ylim()[1] * 0.85
    title_text = ax.text(L/2, title_y, '', fontsize=16, ha='center', weight='bold')

    def draw_spring(line, x1, x2, y=0):
        """Draw a spring as a zigzag line with adaptive coil count"""
        # Scale coil width with system length (L-aware)
        coil_width = min(L * 0.015, mass_spacing * 0.2)

        # Scale number of coils with spring length for realistic appearance
        spring_length = abs(x2 - x1)
        n_coils = max(5, int(spring_length / mass_spacing * 10))

        x = np.linspace(x1, x2, n_coils * 2 + 1)
        y_coords = np.zeros_like(x)
        for i in range(1, len(x) - 1):
            y_coords[i] = coil_width * (1 if i % 2 == 1 else -1)
        y_coords += y
        line.set_data(x, y_coords)

    def init():
        for mass in masses_artists:
            mass.center = (0, 0)
        for line in spring_lines:
            line.set_data([], [])
        title_text.set_text('')
        return masses_artists + spring_lines + [title_text]

    def animate(frame):
        # Calculate which cycle (mode + pause) we're in
        frames_per_cycle = frames_per_mode + frames_per_pause
        cycle_num = frame // frames_per_cycle
        frame_in_cycle = frame % frames_per_cycle

        # Determine if we're in animation or pause
        if frame_in_cycle >= frames_per_mode:
            # We're in pause - freeze at equilibrium
            positions = equilibrium_positions.copy()
            velocities = np.zeros(N)

            # Get the mode that just finished
            if cycle_num < len(modes_to_show):
                current_mode_idx = modes_to_show[cycle_num]
                n = current_mode_idx + 1
                freq_value = frequencies[current_mode_idx] / omega_0
                title_text.set_text(f'Mode {n}: ω = 2ω₀ sin({n}π/{2*(N+1)}) = {freq_value:.4f}ω₀')
            else:
                title_text.set_text('')
        else:
            # We're animating a mode
            if cycle_num >= len(modes_to_show):
                # Animation complete, stay at equilibrium
                positions = equilibrium_positions.copy()
                velocities = np.zeros(N)
                title_text.set_text('')
            else:
                current_mode_idx = modes_to_show[cycle_num]

                # Time within this mode
                t = frame_in_cycle * dt

                # Calculate positions and velocities for current mode
                omega = frequencies[current_mode_idx]
                displacements = amplitude * modes[:, current_mode_idx] * np.sin(omega * t)
                velocities = amplitude * modes[:, current_mode_idx] * omega * np.cos(omega * t)
                positions = equilibrium_positions + displacements

                # Update title with frequency
                n = current_mode_idx + 1
                freq_value = frequencies[current_mode_idx] / omega_0
                title_text.set_text(f'Mode {n}: ω = 2ω₀ sin({n}π/{2*(N+1)}) = {freq_value:.4f}ω₀')

        # Normalize velocities for colormap (relative to max possible velocity)
        max_velocity = amplitude * np.max(frequencies)
        if max_velocity > 0:
            normalized_velocities = velocities / max_velocity
        else:
            normalized_velocities = np.zeros(N)

        # Update masses with position and color
        for i, mass in enumerate(masses_artists):
            mass.center = (positions[i], 0)
            # Color based on velocity
            color = colormap(norm(normalized_velocities[i]))
            mass.set_facecolor(color)

        # Update springs
        # Left wall to first mass
        draw_spring(spring_lines[0], 0, positions[0])

        # Between masses
        for i in range(N - 1):
            draw_spring(spring_lines[i + 1], positions[i], positions[i + 1])

        # Last mass to right wall
        draw_spring(spring_lines[N], positions[N - 1], L)

        return masses_artists + spring_lines + [title_text]

    # Create animation
    # Total frames = (mode animation + pause) for each mode to show
    frames_per_cycle = frames_per_mode + frames_per_pause
    total_frames = frames_per_cycle * len(modes_to_show)

    # Progress callback for video rendering
    def progress_callback(current_frame, total_frames):
        """Print progress during video rendering"""
        if current_frame % 10 == 0 or current_frame == total_frames - 1:
            percent = (current_frame + 1) / total_frames * 100
            print(f"\rRendering: {current_frame + 1}/{total_frames} frames ({percent:.1f}%)", end='', flush=True)
        if current_frame == total_frames - 1:
            print()  # New line at end

    anim = FuncAnimation(fig, animate, init_func=init, frames=total_frames,
                        interval=dt*1000, blit=True, repeat=True)

    plt.tight_layout()

    # Save animation if output file specified
    if args.outfile:
        # Check for ffmpeg availability
        import subprocess
        try:
            subprocess.run(['ffmpeg', '-version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Error: ffmpeg not found. Please install ffmpeg to save animations.")
            print("  macOS: brew install ffmpeg")
            print("  Linux: sudo apt-get install ffmpeg")
            print("  Windows: Download from https://ffmpeg.org/download.html")
            sys.exit(1)

        fps = args.fps if args.fps is not None else int(1/dt)
        print(f"Saving animation to {args.outfile}")
        print(f"Settings: {len(modes_to_show)} mode(s), fps={fps}, dpi={args.dpi}")
        print(f"Total frames: {total_frames}")

        try:
            anim.save(args.outfile, writer='ffmpeg', fps=fps, dpi=args.dpi,
                      progress_callback=progress_callback)
            print(f"\nAnimation saved successfully to {args.outfile}")
        except Exception as e:
            print(f"\nError saving animation: {e}")
            sys.exit(1)
    else:
        plt.show()


if __name__ == '__main__':
    main()
