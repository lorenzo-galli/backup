import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os

# ==============================================================================
# CONFIGURAZIONE PRINCIPALE PER LA VISUALIZZAZIONE
# ==============================================================================
DENSITY = 0.1                         # Densità da visualizzare (es. 0.07, 0.1, 0.2)
TRAJECTORY_FILE = f"run_rho_{DENSITY}_traj.xyz"  # File di traiettoria generato
# ==============================================================================

def read_trajectory(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Impossibile trovare il file di traiettoria '{filename}'. Esegui prima la simulazione con md_simulation.py!")

    with open(filename, "r") as f:
        lines = f.readlines()
        
    N = int(lines[0].strip())
    frames = len(lines) // (N + 2)
    
    trajectory = np.zeros((frames, N, 2))
    
    for frame in range(frames):
        start = frame * (N + 2) + 2
        for i in range(N):
            parts = lines[start + i].split()
            trajectory[frame, i, 0] = float(parts[1])
            trajectory[frame, i, 1] = float(parts[2])
            
    return N, frames, trajectory

def run_visualizer(filename, density):
    print(f"Caricamento traiettoria da '{filename}' (rho = {density})...")
    N, frames, trajectory = read_trajectory(filename)
    
    L = np.sqrt(N / density)
    
    # Wrap delle coordinate nel box primario per la visualizzazione
    trajectory = trajectory % L
    
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_title(f"MD Simulation Live (rho={density})")
    
    # Scatter plot iniziale
    scatter = ax.scatter(trajectory[0, :, 0], trajectory[0, :, 1], c='orange', edgecolors='red', s=50)
    
    def update(frame):
        scatter.set_offsets(trajectory[frame])
        ax.set_title(f"MD Live (rho={density}) - Frame {frame}/{frames}")
        return scatter,
        
    ani = animation.FuncAnimation(fig, update, frames=frames, interval=20, blit=True)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            density = float(sys.argv[1])
            filename = f"run_rho_{density}_traj.xyz"
        except ValueError:
            filename = sys.argv[1]
            density = DENSITY
    else:
        density = DENSITY
        filename = TRAJECTORY_FILE
        
    print(f"Avvio del visualizzatore per densità: {density}")
    print("Suggerimento: puoi anche avviarlo da riga di comando con: python3 md_visualizer.py <densita>\n")
    run_visualizer(filename, density)
