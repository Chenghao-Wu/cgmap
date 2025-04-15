# CGMap

CGMap is a Python tool for coarse-grain mapping of molecular dynamics trajectories. It provides efficient mapping of atomistic trajectories to coarse-grained representations with support for various output formats.

## Features

- Read LAMMPS dump files
- Configurable coarse-graining mapping through YAML files
- Multiple output formats supported:
  - XYZ trajectory files
  - LAMMPS data files
  - Compressed NPZ files
- Coordinate wrapping for periodic boundary conditions
- Detailed logging system
- Command-line interface

## Installation

### Requirements
- Python >= 3.8
- NumPy
- PyYAML

### Installing from source

```bash
git clone https://github.com/yourusername/cgmap.git
cd cgmap
pip install .
```

## Usage

### Basic Command

```bash
cgmap --dump trajectory.dump --system system.yaml --output cg_system
```

### Command-line Options
usage: cgmap [-h] --dump DUMP --system SYSTEM [--output OUTPUT]
[--format {xyz,data,npz}] [--wrap WRAP]
[--separate-xyz] [--frame FRAME]
[--log-dir LOG_DIR] [--debug]
optional arguments:
-h, --help show this help message and exit
--dump DUMP Input LAMMPS dump file
--system SYSTEM System configuration YAML file
--output OUTPUT Output file prefix (default: cg_system)
--wrap WRAP output coordinates in wrapped form (default: True)
--format {xyz,data,npz}
Output format: xyz, data, or npz (default: npz)
--separate-xyz Write separate XYZ files for each frame
--frame FRAME Frame index to write for LAMMPS data file (default: 0)
--log-dir LOG_DIR Directory for log file (default: no file logging)
--debug Enable debug logging

### Example System Configuration (system.yaml)

```yaml
system:
  names:
    - "mapping_wat.yaml"
  numbers:
    - 256
```

### Example Mapping Configuration (mapping_wat.yaml)

```yaml
site-types:
  WAT:
    index: [0, 1, 2]
    x-weight: [16.0, 1.0, 1.0]
    f-weight: [1.0, 1.0, 1.0]

config:
  - anchor: 0
    offset: 3
    repeat: 256
    sites:
      - ["WAT", 0]
```

## Output Formats

### XYZ Format
- Single trajectory file or separate files per frame
- Simple format for visualization
```bash
cgmap --dump traj.dump --system system.yaml --format xyz
```

### LAMMPS Data Format
- Full atom style
- Includes box dimensions and atom types
```bash
cgmap --dump traj.dump --system system.yaml --format data --frame 0
```

### NPZ Format
- Compressed NumPy format
- Contains coordinates, forces, site types, and box dimensions
```bash
cgmap --dump traj.dump --system system.yaml --format npz
```

## Logging

Enable detailed logging with the `--debug` flag and specify a log directory:
```bash
cgmap --dump traj.dump --system system.yaml --log-dir ./logs --debug
```

## Examples

1. Basic coarse-graining with NPZ output:
```bash
cgmap --dump trajectory.dump --system system.yaml
```

2. Generate XYZ trajectory with wrapped coordinates:
```bash
cgmap --dump trajectory.dump --system system.yaml --format xyz --wrap True
```

3. Create LAMMPS data file for a specific frame:
```bash
cgmap --dump trajectory.dump --system system.yaml --format data --frame 10
```

4. Generate separate XYZ files for each frame:
```bash
cgmap --dump trajectory.dump --system system.yaml --format xyz --separate-xyz
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use CGMap in your research, please cite:

@software{cgmap2025,
  author = {Wu, Zhenghao},
  title = {CGMap: A Tool for Coarse-Grain Mapping of Molecular Dynamics Trajectories},
  year = {2024},
  url = {https://github.com/yourusername/cgmap}
}
