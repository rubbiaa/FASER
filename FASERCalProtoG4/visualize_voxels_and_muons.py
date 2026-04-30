#!/usr/bin/env python3
"""
3D Visualization of Voxels and Muon Tracks with Active/Inactive Regions
Parses Geant4 muon tracking data file and visualizes the optically isolated voxels with muon hits
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import argparse

class VoxelVisualizer:
    def __init__(self, data_file, voxel_size=10.0, refl_thickness=0.5, event_id=None):
        """
        Initialize visualizer
        
        Parameters:
        -----------
        data_file : str
            Path to muon tracking data file
        voxel_size : float
            Total voxel size in mm (default: 10 mm)
        refl_thickness : float
            Reflective layer thickness in mm (default: 0.5 mm)
        event_id : int or None
            Specific event ID to visualize, or None for all events
        """
        self.data_file = data_file
        self.voxel_size = voxel_size
        self.refl_thickness = refl_thickness
        self.active_size = voxel_size - 2 * refl_thickness  # 9 mm for 10mm voxel
        self.event_id = event_id
        
        # Detector geometry: 11×48×48 voxels
        self.n_voxels_x = 11
        self.n_voxels_y = 48
        self.n_voxels_z = 48
        
        self.muon_steps = []
        self.electron_steps = []
        self.active_hits = []
        self.inactive_hits = []
        self.electron_active_hits = []
        self.electron_inactive_hits = []
        self.events_found = set()
        
        # Track-by-track data: key = (event_id, track_id, pdg), value = list of steps
        self.tracks = {}
        
    def parse_data_file(self):
        """Parse the muon tracking data file"""
        print(f"Parsing data file: {self.data_file}")
        
        with open(self.data_file, 'r') as f:
            for line in f:
                # Skip comments
                if line.startswith('#'):
                    continue
                    
                # Parse data line
                parts = line.strip().split()
                if len(parts) != 12:
                    continue
                
                try:
                    event_id = int(parts[0])
                    pdg = int(parts[1])
                    track_id = int(parts[2])
                    edep = float(parts[3])
                    gx, gy, gz = float(parts[4]), float(parts[5]), float(parts[6])
                    lx, ly, lz = float(parts[7]), float(parts[8]), float(parts[9])
                    layer = int(parts[10])
                    is_active = bool(int(parts[11]))
                    
                    # Filter by event ID if specified
                    if self.event_id is not None and event_id != self.event_id:
                        continue
                    
                    self.events_found.add(event_id)
                    
                    step_data = {
                        'event_id': event_id,
                        'pdg': pdg,
                        'track_id': track_id,
                        'edep': edep,
                        'global_pos': np.array([gx, gy, gz]),
                        'local_pos': np.array([lx, ly, lz]),
                        'layer': layer,
                        'is_active': is_active
                    }
                    
                    # Store by track
                    track_key = (event_id, track_id, pdg)
                    if track_key not in self.tracks:
                        self.tracks[track_key] = []
                    self.tracks[track_key].append(step_data)
                    
                    # Separate muons and electrons by PDG code
                    # PDG: 13/-13 = muon/antimuon, 11/-11 = electron/positron
                    if abs(pdg) == 13:  # Muon
                        self.muon_steps.append(step_data)
                        if step_data['is_active']:
                            self.active_hits.append(step_data)
                        else:
                            self.inactive_hits.append(step_data)
                    elif abs(pdg) == 11:  # Electron
                        self.electron_steps.append(step_data)
                        if step_data['is_active']:
                            self.electron_active_hits.append(step_data)
                        else:
                            self.electron_inactive_hits.append(step_data)
                    else:
                        # Other particles - treat as muon for now
                        self.muon_steps.append(step_data)
                        if step_data['is_active']:
                            self.active_hits.append(step_data)
                        else:
                            self.inactive_hits.append(step_data)
                        
                except (ValueError, IndexError) as e:
                    continue
        
        print(f"Found {len(self.events_found)} event(s): {sorted(self.events_found)}")
        print(f"Total muon steps: {len(self.muon_steps)}")
        print(f"  - {len(self.active_hits)} in active regions")
        print(f"Total electron steps: {len(self.electron_steps)}")
        print(f"  - {len(self.electron_active_hits)} in active regions")
        print(f"  - {len(self.inactive_hits)} in inactive (reflective) regions")
        
    def draw_voxel_box(self, ax, center, size, color='blue', alpha=0.1, edge_color='black', linewidth=0.5):
        """Draw a 3D box representing a voxel or detector outline"""
        cx, cy, cz = center
        # Handle both scalar and array sizes
        if np.isscalar(size):
            dx = dy = dz = size / 2
        else:
            dx, dy, dz = np.array(size) / 2
        
        # Define the 8 vertices of the box
        vertices = np.array([
            [cx-dx, cy-dy, cz-dz],
            [cx+dx, cy-dy, cz-dz],
            [cx+dx, cy+dy, cz-dz],
            [cx-dx, cy+dy, cz-dz],
            [cx-dx, cy-dy, cz+dz],
            [cx+dx, cy-dy, cz+dz],
            [cx+dx, cy+dy, cz+dz],
            [cx-dx, cy+dy, cz+dz]
        ])
        
        # Define the 6 faces
        faces = [
            [vertices[0], vertices[1], vertices[5], vertices[4]],  # front
            [vertices[2], vertices[3], vertices[7], vertices[6]],  # back
            [vertices[0], vertices[3], vertices[7], vertices[4]],  # left
            [vertices[1], vertices[2], vertices[6], vertices[5]],  # right
            [vertices[0], vertices[1], vertices[2], vertices[3]],  # bottom
            [vertices[4], vertices[5], vertices[6], vertices[7]]   # top
        ]
        
        # Only create face collection if we have a valid color
        if color != 'none' and alpha > 0:
            face_collection = Poly3DCollection(faces, alpha=alpha, facecolor=color, 
                                              edgecolor=edge_color, linewidth=linewidth)
            ax.add_collection3d(face_collection)
        else:
            # Just draw the edges as lines for outline
            edges = [
                [vertices[0], vertices[1]], [vertices[1], vertices[2]],
                [vertices[2], vertices[3]], [vertices[3], vertices[0]],
                [vertices[4], vertices[5]], [vertices[5], vertices[6]],
                [vertices[6], vertices[7]], [vertices[7], vertices[4]],
                [vertices[0], vertices[4]], [vertices[1], vertices[5]],
                [vertices[2], vertices[6]], [vertices[3], vertices[7]]
            ]
            for edge in edges:
                edge_array = np.array(edge)
                ax.plot3D(edge_array[:, 0], edge_array[:, 1], edge_array[:, 2],
                         color=edge_color, linewidth=linewidth)
        
    def draw_voxel_with_reflective_layers(self, ax, voxel_idx):
        """Draw a single voxel showing active and reflective regions"""
        ix, iy, iz = voxel_idx
        
        # Center of the voxel
        center = np.array([
            ix * self.voxel_size + self.voxel_size / 2,
            iy * self.voxel_size + self.voxel_size / 2,
            iz * self.voxel_size + self.voxel_size / 2
        ])
        
        # Draw outer box (total voxel with reflective layer)
        self.draw_voxel_box(ax, center, self.voxel_size, 
                           color='gray', alpha=0.05, edge_color='gray', linewidth=0.3)
        
        # Draw inner box (active scintillator region)
        self.draw_voxel_box(ax, center, self.active_size, 
                           color='cyan', alpha=0.15, edge_color='blue', linewidth=0.8)
        
    def plot_3d_visualization(self):
        """Create 3D visualization of voxels and muon hits"""
        if not self.muon_steps:
            print("No muon steps found. Please run the simulation first.")
            return
        
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111, projection='3d')
        
        # Calculate detector dimensions for reference
        det_size_x = self.n_voxels_x * self.voxel_size
        det_size_y = self.n_voxels_y * self.voxel_size
        det_size_z = self.n_voxels_z * self.voxel_size
        
        print(f"\nFull detector geometry:")
        print(f"  Voxels: {self.n_voxels_x} × {self.n_voxels_y} × {self.n_voxels_z}")
        print(f"  Size: {det_size_x} × {det_size_y} × {det_size_z} mm³")
        
        # Get range of global positions from hits
        all_global_pos = np.array([s['global_pos'] for s in self.muon_steps])
        min_global = np.min(all_global_pos, axis=0)
        max_global = np.max(all_global_pos, axis=0)
        
        print(f"\nHit region (mm):")
        print(f"  X: [{min_global[0]:.2f}, {max_global[0]:.2f}]")
        print(f"  Y: [{min_global[1]:.2f}, {max_global[1]:.2f}]")
        print(f"  Z: [{min_global[2]:.2f}, {max_global[2]:.2f}]")
        
        # Calculate which voxels are hit (in global coordinates)
        voxels_hit = set()
        for step in self.muon_steps + self.electron_steps:
            gpos = step['global_pos']
            # Calculate voxel indices from global position
            ix = int((gpos[0] - min_global[0]) // self.voxel_size)
            iy = int((gpos[1] - min_global[1]) // self.voxel_size)
            iz = int((gpos[2] - min_global[2]) // self.voxel_size)
            voxels_hit.add((ix, iy, iz))
        
        print(f"\nVoxels with hits: {len(voxels_hit)}")
        
        # Add neighboring voxels with multiple layers for better context
        # Use distance 1 to show 3x3x3 region around each hit for better XY plane visibility
        neighbor_distance = 0  # How many voxels out to show neighbors
        voxels_neighbors = set()
        for ix, iy, iz in voxels_hit:
            # Add all neighbors within distance
            for dx in range(-neighbor_distance, neighbor_distance + 1):
                for dy in range(-neighbor_distance, neighbor_distance + 1):
                    for dz in range(-neighbor_distance, neighbor_distance + 1):
                        if dx == 0 and dy == 0 and dz == 0:
                            continue  # Skip the center voxel itself
                        neighbor = (ix+dx, iy+dy, iz+dz)
                        # Check if neighbor is within reasonable bounds
                        nx, ny, nz = neighbor
                        if nx >= 0 and ny >= 0 and nz >= 0:  # Basic bounds check
                            voxels_neighbors.add(neighbor)
        
        # Remove voxels that already have hits from the neighbors set
        voxels_neighbors -= voxels_hit
        
        print(f"Neighboring voxels (within distance {neighbor_distance}): {len(voxels_neighbors)}")
        print(f"Drawing {len(voxels_hit)} hit voxels + {len(voxels_neighbors)} neighboring voxels")
        
        # Draw neighboring voxels first (so they appear behind hit voxels)
        for voxel_idx in sorted(voxels_neighbors):
            ix, iy, iz = voxel_idx
            # Calculate center in global coordinates
            center = np.array([
                min_global[0] + ix * self.voxel_size + self.voxel_size / 2,
                min_global[1] + iy * self.voxel_size + self.voxel_size / 2,
                min_global[2] + iz * self.voxel_size + self.voxel_size / 2
            ])
            
            # Draw outer box (total voxel) - very subtle wireframe
            self.draw_voxel_box(ax, center, self.voxel_size, 
                               color='lightgray', alpha=0.01, 
                               edge_color='lightgray', linewidth=0.3)
            
            # Draw inner box (active region) - cyan wireframe for neighbors
            self.draw_voxel_box(ax, center, self.active_size, 
                               color='cyan', alpha=0.1, 
                               edge_color='cyan', linewidth=0.6)
        
        # Draw voxels with hits (more prominent)
        voxels_to_draw = sorted(voxels_hit)
        
        # Draw voxels
        for voxel_idx in voxels_to_draw:
            ix, iy, iz = voxel_idx
            # Calculate center in global coordinates
            center = np.array([
                min_global[0] + ix * self.voxel_size + self.voxel_size / 2,
                min_global[1] + iy * self.voxel_size + self.voxel_size / 2,
                min_global[2] + iz * self.voxel_size + self.voxel_size / 2
            ])
            
            # Draw outer box (total voxel volume including reflective layer)
            self.draw_voxel_box(ax, center, self.voxel_size, 
                               color='lightgray', alpha=0.03, 
                               edge_color='gray', linewidth=0.4)
            
            # Draw inner box (active scintillator region) - prominent blue
            self.draw_voxel_box(ax, center, self.active_size, 
                               color='cyan', alpha=0.3, 
                               edge_color='blue', linewidth=1.0)
        
        # Plot muon active hits (use global positions)
        if self.active_hits:
            active_pos = np.array([s['global_pos'] for s in self.active_hits])
            ax.scatter(active_pos[:, 0], active_pos[:, 1], active_pos[:, 2],
                      c='green', marker='o', s=50, alpha=0.8, label='Muon active hits', edgecolors='darkgreen')
        
        # Plot muon inactive hits (use global positions)
        if self.inactive_hits:
            inactive_pos = np.array([s['global_pos'] for s in self.inactive_hits])
            ax.scatter(inactive_pos[:, 0], inactive_pos[:, 1], inactive_pos[:, 2],
                      c='red', marker='x', s=100, alpha=0.9, label='Muon inactive hits', linewidths=2)
        
        # Plot electron active hits in violet (use global positions)
        if self.electron_active_hits:
            e_active_pos = np.array([s['global_pos'] for s in self.electron_active_hits])
            ax.scatter(e_active_pos[:, 0], e_active_pos[:, 1], e_active_pos[:, 2],
                      c='violet', marker='o', s=50, alpha=0.8, label='Electron active hits', edgecolors='purple')
        
        # Plot electron inactive hits in violet (use global positions)
        if self.electron_inactive_hits:
            e_inactive_pos = np.array([s['global_pos'] for s in self.electron_inactive_hits])
            ax.scatter(e_inactive_pos[:, 0], e_inactive_pos[:, 1], e_inactive_pos[:, 2],
                      c='darkviolet', marker='x', s=100, alpha=0.9, label='Electron inactive hits', linewidths=2)
        
        # Draw muon track (use global positions)
        if len(self.muon_steps) > 1:
            track_pos = np.array([s['global_pos'] for s in self.muon_steps])
            ax.plot(track_pos[:, 0], track_pos[:, 1], track_pos[:, 2],
                   color='cyan', linewidth=2, alpha=0.7, label='Muon track')
        
        # Draw electron track in violet (use global positions)
        if len(self.electron_steps) > 1:
            e_track_pos = np.array([s['global_pos'] for s in self.electron_steps])
            ax.plot(e_track_pos[:, 0], e_track_pos[:, 1], e_track_pos[:, 2],
                   color='violet', linestyle='-', linewidth=2, alpha=0.7, label='Electron track')
        
        # Labels and styling
        ax.set_xlabel('X (mm)', fontsize=12)
        ax.set_ylabel('Y (mm)', fontsize=12)
        ax.set_zlabel('Z (mm)', fontsize=12)
        
        # Create title based on whether we're showing specific event or all events
        if self.event_id is not None:
            title = f'3D Voxel Structure with Particle Hits - Event {self.event_id}\n'
        elif len(self.events_found) == 1:
            title = f'3D Voxel Structure with Particle Hits - Event {list(self.events_found)[0]}\n'
        else:
            title = f'3D Voxel Structure with Particle Hits - {len(self.events_found)} Events\n'
        
        title += (f'Voxel size: {self.voxel_size} mm (Active: {self.active_size} mm, '
                 f'Reflective: {self.refl_thickness} mm on each face)')
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Set view limits to focus on hit region with some margin
        margin = 0.2  # 20% margin around hit region
        range_x = (max_global[0] - min_global[0]) * (1 + margin)
        range_y = (max_global[1] - min_global[1]) * (1 + margin)
        range_z = (max_global[2] - min_global[2]) * (1 + margin)
        
        # Center of hit region
        mid = (min_global + max_global) / 2
        
        # Use the largest range for all axes to maintain aspect ratio
        max_range = max(range_x, range_y, range_z) / 2
        
        ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
        ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
        ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
        
        plt.tight_layout()
        
        # Save figure
        if self.event_id is not None:
            output_file = f'voxel_visualization_event{self.event_id}_3d.png'
        else:
            output_file = 'voxel_visualization_all_events_3d.png'
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nVisualization saved to: {output_file}")
        
        plt.show()
        
    def print_statistics(self):
        """Print statistics about particle passage"""
        if not self.muon_steps and not self.electron_steps:
            return
        
        # Muon statistics
        muon_total_edep = sum(s['edep'] for s in self.muon_steps)
        muon_active_edep = sum(s['edep'] for s in self.active_hits)
        muon_inactive_edep = sum(s['edep'] for s in self.inactive_hits)
        
        # Electron statistics
        electron_total_edep = sum(s['edep'] for s in self.electron_steps)
        electron_active_edep = sum(s['edep'] for s in self.electron_active_hits)
        electron_inactive_edep = sum(s['edep'] for s in self.electron_inactive_hits)
        
        total_edep = muon_total_edep + electron_total_edep
        active_edep = muon_active_edep + electron_active_edep
        inactive_edep = muon_inactive_edep + electron_inactive_edep
        
        print("\n" + "="*60)
        print("PARTICLE TRACKING STATISTICS")
        print("="*60)
        
        if self.event_id is not None:
            print(f"Event ID: {self.event_id}")
        else:
            print(f"Events analyzed: {len(self.events_found)}")
            if len(self.events_found) <= 10:
                print(f"Event IDs: {sorted(self.events_found)}")
        
        # Muon statistics
        if self.muon_steps:
            print(f"\n--- MUONS ---")
            print(f"Total muon steps: {len(self.muon_steps)}")
            print(f"  Active region hits: {len(self.active_hits)} ({100*len(self.active_hits)/len(self.muon_steps):.1f}%)")
            print(f"  Inactive region hits: {len(self.inactive_hits)} ({100*len(self.inactive_hits)/len(self.muon_steps):.1f}%)")
            print(f"  Energy deposition: {muon_total_edep:.3f} MeV")
            print(f"    In active: {muon_active_edep:.3f} MeV ({100*muon_active_edep/muon_total_edep:.1f}%)")
            print(f"    In inactive: {muon_inactive_edep:.3f} MeV ({100*muon_inactive_edep/muon_total_edep:.1f}%)")
        
        # Electron statistics
        if self.electron_steps:
            print(f"\n--- ELECTRONS ---")
            print(f"Total electron steps: {len(self.electron_steps)}")
            print(f"  Active region hits: {len(self.electron_active_hits)} ({100*len(self.electron_active_hits)/len(self.electron_steps):.1f}%)")
            print(f"  Inactive region hits: {len(self.electron_inactive_hits)} ({100*len(self.electron_inactive_hits)/len(self.electron_steps):.1f}%)")
            print(f"  Energy deposition: {electron_total_edep:.3f} MeV")
            print(f"    In active: {electron_active_edep:.3f} MeV ({100*electron_active_edep/electron_total_edep:.1f}%)")
            print(f"    In inactive: {electron_inactive_edep:.3f} MeV ({100*electron_inactive_edep/electron_total_edep:.1f}%)")
        
        # Combined statistics
        if self.muon_steps and self.electron_steps:
            print(f"\n--- COMBINED ---")
            print(f"Total energy deposition: {total_edep:.3f} MeV")
            print(f"  In active regions: {active_edep:.3f} MeV ({100*active_edep/total_edep:.1f}%)")
            print(f"  In inactive regions: {inactive_edep:.3f} MeV ({100*inactive_edep/total_edep:.1f}%)")
        print(f"\nExpected active volume fraction: {(self.active_size/self.voxel_size)**3:.1%}")
        print("="*60)
        
    def analyze_tracks(self):
        """Analyze each track separately for active/inactive hit fractions"""
        if not self.tracks:
            return
        
        print("\n" + "="*80)
        print("TRACK-BY-TRACK ANALYSIS")
        print("="*80)
        print(f"Total tracks analyzed: {len(self.tracks)}")
        print()
        
        track_stats = []
        
        for track_key, steps in self.tracks.items():
            event_id, track_id, pdg = track_key
            
            if len(steps) < 2:
                continue  # Need at least 2 points for track length
            
            # Calculate track length
            positions = np.array([s['global_pos'] for s in steps])
            segments = np.diff(positions, axis=0)
            segment_lengths = np.linalg.norm(segments, axis=1)
            track_length = np.sum(segment_lengths)
            
            # Count active and inactive hits
            n_active = sum(1 for s in steps if s['is_active'])
            n_inactive = sum(1 for s in steps if not s['is_active'])
            n_total = len(steps)
            
            # Calculate fractions
            frac_active = n_active / n_total if n_total > 0 else 0
            frac_inactive = n_inactive / n_total if n_total > 0 else 0
            
            # Hits per unit length (hits per mm)
            hits_per_mm_active = n_active / track_length if track_length > 0 else 0
            hits_per_mm_inactive = n_inactive / track_length if track_length > 0 else 0
            hits_per_mm_total = n_total / track_length if track_length > 0 else 0
            
            # Energy statistics
            total_edep = sum(s['edep'] for s in steps)
            active_edep = sum(s['edep'] for s in steps if s['is_active'])
            inactive_edep = sum(s['edep'] for s in steps if not s['is_active'])
            
            # Calculate maximum gap between active/inactive transitions
            max_gap = 0.0
            if len(steps) > 1:
                for i in range(len(steps) - 1):
                    if steps[i]['is_active'] != steps[i+1]['is_active']:
                        # Transition detected - calculate distance
                        gap = np.linalg.norm(steps[i+1]['global_pos'] - steps[i]['global_pos'])
                        max_gap = max(max_gap, gap)
            
            particle_name = 'muon' if abs(pdg) == 13 else ('electron' if abs(pdg) == 11 else f'PDG{pdg}')
            
            track_stats.append({
                'event_id': event_id,
                'track_id': track_id,
                'pdg': pdg,
                'particle': particle_name,
                'n_steps': n_total,
                'n_active': n_active,
                'n_inactive': n_inactive,
                'frac_active': frac_active,
                'frac_inactive': frac_inactive,
                'track_length': track_length,
                'hits_per_mm_total': hits_per_mm_total,
                'hits_per_mm_active': hits_per_mm_active,
                'hits_per_mm_inactive': hits_per_mm_inactive,
                'total_edep': total_edep,
                'active_edep': active_edep,
                'inactive_edep': inactive_edep,
                'max_gap': max_gap
            })
        
        # Sort by track length (longest first)
        track_stats.sort(key=lambda x: x['track_length'], reverse=True)
        
        # Print summary statistics
        all_frac_active = [t['frac_active'] for t in track_stats]
        all_frac_inactive = [t['frac_inactive'] for t in track_stats]
        all_track_lengths = [t['track_length'] for t in track_stats]
        
        print(f"Track length statistics:")
        print(f"  Mean: {np.mean(all_track_lengths):.2f} mm")
        print(f"  Median: {np.median(all_track_lengths):.2f} mm")
        print(f"  Min: {np.min(all_track_lengths):.2f} mm")
        print(f"  Max: {np.max(all_track_lengths):.2f} mm")
        print()
        print(f"Fraction of hits in ACTIVE regions:")
        print(f"  Mean: {np.mean(all_frac_active):.1%}")
        print(f"  Median: {np.median(all_frac_active):.1%}")
        print(f"  Std: {np.std(all_frac_active):.1%}")
        print()
        print(f"Fraction of hits in INACTIVE regions:")
        print(f"  Mean: {np.mean(all_frac_inactive):.1%}")
        print(f"  Median: {np.median(all_frac_inactive):.1%}")
        print(f"  Std: {np.std(all_frac_inactive):.1%}")
        print()
        
        # Print top 10 longest tracks
        print(f"Top 10 longest tracks:")
        print(f"{'Event':>6} {'Track':>6} {'Particle':>10} {'Length(mm)':>12} {'Steps':>6} "
              f"{'Active':>6} {'Inactive':>8} {'Frac.Act':>9} {'Frac.Inact':>11} {'Edep(MeV)':>11}")
        print("-" * 100)
        
        for i, track in enumerate(track_stats[:10]):
            print(f"{track['event_id']:>6} {track['track_id']:>6} {track['particle']:>10} "
                  f"{track['track_length']:>12.2f} {track['n_steps']:>6} "
                  f"{track['n_active']:>6} {track['n_inactive']:>8} "
                  f"{track['frac_active']:>9.1%} {track['frac_inactive']:>11.1%} "
                  f"{track['total_edep']:>11.3f}")
        
        print("\n" + "="*80)
        
        # Print event-by-event summary
        print("\n" + "="*80)
        print("EVENT-BY-EVENT SUMMARY")
        print("="*80)
        print(f"{'EventID':>8} {'PDG':>6} {'Particle':>10} {'Steps':>6} {'Active':>7} {'Inactive':>9} "
              f"{'Frac.Act':>9} {'MaxGap(mm)':>11} {'Edep(MeV)':>11}")
        print("-" * 80)
        
        # Sort by event ID for event-by-event listing
        track_stats_by_event = sorted(track_stats, key=lambda x: x['event_id'])
        
        for track in track_stats_by_event:
            print(f"{track['event_id']:>8} {track['pdg']:>6} {track['particle']:>10} "
                  f"{track['n_steps']:>6} {track['n_active']:>7} {track['n_inactive']:>9} "
                  f"{track['frac_active']:>9.1%} {track['max_gap']:>11.2f} "
                  f"{track['total_edep']:>11.3f}")
        
        print("="*80)
        
        return track_stats

def main():
    parser = argparse.ArgumentParser(description='Visualize voxel structure and muon tracks')
    parser.add_argument('data_file', type=str, help='Path to muon tracking data file')
    parser.add_argument('--voxel-size', type=float, default=10.0, 
                       help='Total voxel size in mm (default: 10.0)')
    parser.add_argument('--refl-thickness', type=float, default=0.5,
                       help='Reflective layer thickness in mm (default: 0.5)')
    parser.add_argument('--event-id', type=int, default=None,
                       help='Specific event ID to visualize (default: all events)')
    
    args = parser.parse_args()
    
    # Create visualizer
    viz = VoxelVisualizer(args.data_file, args.voxel_size, args.refl_thickness, args.event_id)
    
    # Parse data file
    viz.parse_data_file()
    
    # Print statistics
    viz.print_statistics()
    
    # Analyze tracks
    viz.analyze_tracks()
    
    # Create 3D visualization (focuses on hit region only)
    viz.plot_3d_visualization()

if __name__ == '__main__':
    main()
