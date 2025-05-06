import yaml
from cgmap.read_lammps import *
import sys
import numpy as np
import argparse
import logging
import os
from datetime import datetime

"""
usage: cgmap.py [-h] --dump DUMP --system SYSTEM [--output OUTPUT]
                [--format {xyz,data,both}] [--separate-xyz] [--frame FRAME]

Coarse-grain mapping tool for molecular dynamics trajectories

optional arguments:
  -h, --help            show this help message and exit
  --dump DUMP           Input LAMMPS dump file
  --system SYSTEM       System configuration YAML file
  --output OUTPUT       Output file prefix (default: cg_system)
  --format {xyz,data,npz}
                        Output format: xyz, data, or npz (default: npz)
  --separate-xyz        Write separate XYZ files for each frame
  --frame FRAME         Frame index to write for LAMMPS data file (default: 0)
  --target all
"""

def setup_logger(output_dir=None, debug=False):
    """
    Set up logger with file and console handlers.
    
    Args:
        output_dir (str): Directory for log file (default: None)
        debug (bool): If True, set log level to DEBUG
    """
    # Create logger
    logger = logging.getLogger('cgmap')
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    
    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_formatter = logging.Formatter(
        '%(levelname)s - %(message)s'
    )
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if output directory is specified
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = os.path.join(output_dir, f'cgmap_{timestamp}.log')
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger

def read_mapping_file(filepath):
    """
    Read and parse a YAML mapping file.
    
    Args:
        filepath (str): Path to the YAML file
        
    Returns:
        dict: Parsed YAML content
    """
    logger = logging.getLogger('cgmap')
    try:
        logger.info(f"Reading mapping file: {filepath}")
        with open(filepath, 'r') as file:
            mapping_data = yaml.safe_load(file)
        logger.debug(f"Successfully parsed mapping file")
        return mapping_data
    except FileNotFoundError:
        logger.error(f"Error: File '{filepath}' not found")
        return None
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML file: {e}")
        return None

def cg_map(frame_data,system_data,target='all'):
    logger = logging.getLogger('cgmap')
    logger.info("Starting coarse-graining mapping")
    
    xyz = get_frame(frame_data,'xyz')
    if target == 'coord':
        forces=np.array([])
    elif target == 'all':
        forces = get_frame(frame_data,'force')
    else:
        logger.error(f"Error selecting target: {target}")

    ixiyiz = get_frame(frame_data,'ixiyiz')
    box = get_frame(frame_data,'box')
    atom_id=0

    cg_group_site = []
    cg_group_coord = []
    cg_group_force = []
    cg_box = box

    logger.debug(f"Input data shapes - xyz: {xyz.shape}, forces: {forces.shape}")

    for index, group_i in enumerate(system_data['system']['numbers']):
        mapping_data = read_mapping_file(system_data['system']['names'][index])
        system_config = mapping_data['config'][0]
        repeat = system_config['repeat']  # 256
        offset = system_config['offset']
        total_atoms_group = repeat * offset

        for repeat_i in range(repeat):
            _id_start = repeat_i*offset + atom_id
            _id_end = (repeat_i+1)*offset + atom_id
            _xyz = xyz[_id_start:_id_end]
            _forces = forces[_id_start:_id_end]
            _ixiyiz = ixiyiz[_id_start:_id_end]
            #print(_id_start,_id_end,xyz.shape,_xyz.shape)
            for site_i in system_config['sites']:
                mol_name = site_i[0]
                anchor = site_i[1]
                mol_data = mapping_data['site-types'][mol_name]
                index_values = mol_data['index']  

                _id_site_start = anchor
                _id_site_end = anchor+len(index_values)

                _site_xyz = _xyz[_id_site_start:_id_site_end]
                _site_forces = _forces[_id_site_start:_id_site_end]
                _site_ixiyiz = _ixiyiz[_id_site_start:_id_site_end]
                
                x_weights = np.array(mol_data['x-weight'] )[:,np.newaxis].transpose()
                f_weights = np.array(mol_data['f-weight'] )[:,np.newaxis].transpose()

                _unwrap_coord = _site_ixiyiz * box + _site_xyz 
                #print(_unwrap_coord)
                if target == 'coord':
                    _x_weighted_sum = np.einsum("fd,cf->cd",_unwrap_coord,x_weights)
                    _x_weighted = _x_weighted_sum/x_weights.sum()
                    _cg_coord = _x_weighted
                    _cg_forces = np.array([])
                elif target == 'all':
                    _x_weighted_sum = np.einsum("fd,cf->cd",_unwrap_coord,x_weights)
                    _x_weighted = _x_weighted_sum/x_weights.sum()
                    _cg_coord = _x_weighted
                    _cg_forces = np.einsum("fd,cf->cd",_site_forces,f_weights)
                else:
                    logger.error(f"Error selecting target: {target}")
        
                cg_group_site.append(mol_name)
                cg_group_coord.append(_cg_coord.reshape(-1))
                cg_group_force.append(_cg_forces.reshape(-1))
        

        atom_id = atom_id + total_atoms_group
    logger.info(f"Mapped {len(cg_group_coord)} CG beads")
    return cg_group_site, cg_group_coord, cg_group_force, cg_box


# cg_data={'R':{},'z':{},'F':{},'cell':{}}

# for index,time_step in enumerate(list(dump_data.keys())):
#     frame_data = dump_data[time_step]
#     cg_group_site, cg_group_coord, cg_group_force, cg_box = cg_map(frame_data,system_data)
#     cg_data['R'][index] = cg_group_coord
#     cg_data['F'][index] = cg_group_site
#     cg_data['z'][index] = cg_group_force
#     cg_data['cell'][index] = cg_box

def write_xyz_trajectory(out_file, cg_data):
    """
    Write all frames of coarse-grained beads to a single XYZ trajectory file.
    
    Args:
        out_file (str): Path to the output XYZ file
        cg_data (dict): Dictionary containing CG data with keys 'R' for coordinates
                       and 'F' for bead types for each frame
    """
    try:
        with open(out_file, 'w') as f:
            n_frames = len(cg_data['R'])
            for frame in range(n_frames):
                coords = cg_data['R'][frame]
                sites = cg_data['F'][frame]
                
                # Write number of atoms
                n_beads = len(sites)
                f.write(f"{n_beads}\n")
                
                # Write comment line with frame number
                f.write(f"Frame {frame}\n")
                
                # Write coordinates
                for site, coord in zip(sites, coords):
                    x, y, z = coord
                    f.write(f"{site:4s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
                
    except Exception as e:
        print(f"Error writing XYZ trajectory file: {e}")
        return False
        
    return True

# Alternative function to write separate files for each frame
def write_xyz_separate_frames(prefix, cg_data):
    """
    Write each frame to a separate XYZ file.
    
    Args:
        prefix (str): Prefix for output files (will be appended with frame number)
        cg_data (dict): Dictionary containing CG data with keys 'R' for coordinates
                       and 'F' for bead types for each frame
    """
    try:
        n_frames = len(cg_data['R'])
        for frame in range(n_frames):
            filename = f"{prefix}_{frame:04d}.xyz"
            coords = cg_data['R'][frame]
            sites = cg_data['F'][frame]
            
            with open(filename, 'w') as f:
                # Write number of atoms
                n_beads = len(sites)
                f.write(f"{n_beads}\n")
                
                # Write comment line
                f.write(f"Frame {frame}\n")
                
                # Write coordinates
                for site, coord in zip(sites, coords):
                    x, y, z = coord
                    f.write(f"{site:4s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
                
    except Exception as e:
        print(f"Error writing XYZ files: {e}")
        return False
        
    return True

def write_lammps_data(out_file, cg_data, frame_idx=0, mol_id=None):
    """
    Write coarse-grained beads to LAMMPS data file in 'full' atom style.
    
    Args:
        out_file (str): Path to the output LAMMPS data file
        cg_data (dict): Dictionary containing CG data with keys:
                       'R' for coordinates
                       'F' for bead types
                       'cell' for box dimensions
        frame_idx (int): Frame index to write (default: 0)
        mol_id (array-like): Molecule IDs for each bead (default: all 1)
    """
    logger = logging.getLogger('cgmap')
    try:
        logger.info(f"Writing LAMMPS data file: {out_file}")
        coords = cg_data['R'][frame_idx]
        sites = cg_data['F'][frame_idx]
        box = cg_data['cell'][frame_idx]
        n_atoms = len(coords)
        
        logger.debug(f"Number of atoms: {n_atoms}")
        
        # Create molecule IDs if not provided
        if mol_id is None:
            mol_id = np.ones(n_atoms, dtype=int)
            logger.debug("Using default molecule IDs (all 1)")
            
        # Get unique atom types and create type mapping
        unique_types = sorted(set(sites))
        type_dict = {site: idx+1 for idx, site in enumerate(unique_types)}
        
        with open(out_file, 'w') as f:
            # Write header
            f.write("LAMMPS data file for coarse-grained system\n\n")
            
            # Write system dimensions
            f.write(f"{n_atoms} atoms\n")
            f.write(f"{len(unique_types)} atom types\n\n")
            
            # Write box dimensions
            f.write(f"0.0 {box[0]} xlo xhi\n")
            f.write(f"0.0 {box[1]} ylo yhi\n")
            f.write(f"0.0 {box[2]} zlo zhi\n\n")
            
            # Write masses (placeholder masses of 1.0)
            f.write("Masses\n\n")
            for i, site in enumerate(unique_types, 1):
                f.write(f"{i} 1.0  # {site}\n")
            f.write("\n")
            
            # Write atoms section
            # Format: atom-ID molecule-ID atom-type charge x y z
            f.write("Atoms  # full\n\n")
            for i, (site, coord, mol) in enumerate(zip(sites, coords, mol_id), 1):
                atom_type = type_dict[site]
                x, y, z = coord
                # Using 0.0 as placeholder for charge
                f.write(f"{i} {mol} {atom_type} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")
            
        logger.debug("Successfully wrote LAMMPS data file")
        return True
        
    except Exception as e:
        logger.error(f"Error writing LAMMPS data file: {e}")
        return False

# Add a new function to save npz file:
def save_cg_data_npz(filename, cg_data):
    """
    Save coarse-grained data to NPZ file format.
    
    Args:
        filename (str): Output NPZ filename
        cg_data (dict): Dictionary containing CG data with keys:
                       'R' for coordinates
                       'F' for bead types
                       'z' for forces
                       'cell' for box dimensions
    """
    try:
        # Convert site types (F) from list to array for better storage
        n_frames = len(cg_data['R'])
        site_types = np.array(list(cg_data['F'].values()), dtype='U4')  # U4 for 4-character strings
        
        # Save to npz file
        np.savez_compressed(filename,
                          coordinates=np.array(list(cg_data['R'].values())),
                          forces=np.array(list(cg_data['z'].values())),
                          site_types=site_types,
                          cell=np.array(list(cg_data['cell'].values())),
                          n_frames=n_frames)
        return True
    except Exception as e:
        print(f"Error saving NPZ file: {e}")
        return False

def wrap_coordinates(coords, box):
    """
    Wrap coordinates back into primary simulation box.
    
    Args:
        coords (list/numpy.array): List of coordinates to wrap, shape (n_atoms, 3)
        box (numpy.array): Box dimensions [Lx, Ly, Lz]
        
    Returns:
        numpy.array: Wrapped coordinates with same shape as input
    """
    coords = np.array(coords)  # Ensure numpy array
    wrapped_coords = coords.copy()
    
    # Wrap coordinates for each dimension
    for i in range(3):  # x, y, z dimensions
        # Subtract box length until coordinate is less than box length
        wrapped_coords[:, i] = wrapped_coords[:, i] - box[i] * np.floor(wrapped_coords[:, i] / box[i])
        
    return wrapped_coords

def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Coarse-grain mapping tool for molecular dynamics trajectories')
    
    # Input files
    parser.add_argument('--dump', type=str, required=True,
                      help='Input LAMMPS dump file')
    parser.add_argument('--system', type=str, required=True,
                      help='System configuration YAML file')
    
    # Output options
    parser.add_argument('--output', type=str, default='cg_system',
                      help='Output file prefix (default: cg_system)')
    parser.add_argument('--wrap', type=bool, default=True,
                      help='output coordinates in wrapped form (default: True)')
    parser.add_argument('--format', choices=['xyz', 'data', 'npz'], default='npz',
                      help='Output format: xyz, data, or npz (default: npz)')
    
    # XYZ output options
    parser.add_argument('--separate-xyz', action='store_true',
                      help='Write separate XYZ files for each frame')
    
    parser.add_argument('--target', type=str, default='all',
                      help='select the target property for coarse-graining, default: all (coordinates and forces)')
    
    # LAMMPS data output options
    parser.add_argument('--frame', type=int, default=0,
                      help='Frame index to write for LAMMPS data file (default: 0)')
    
    # Add logging-related arguments
    parser.add_argument('--log-dir', type=str, default=None,
                      help='Directory for log file (default: no file logging)')
    parser.add_argument('--debug', action='store_true',
                      help='Enable debug logging')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Setup logger
    logger = setup_logger(args.log_dir, args.debug)
    logger.info("Starting CGMap")
    logger.debug(f"Arguments: {args}")
    
    # Read system configuration
    system_data = read_mapping_file(args.system)
    if not system_data:
        logger.error("Failed to read system configuration file")
        sys.exit(1)
    
    # Read dump file
    logger.info(f"Reading dump file: {args.dump}")
    dump_data = read_dump_file(args.dump)
    if not dump_data:
        logger.error("Failed to read dump file")
        sys.exit(1)
    
    # Perform CG mapping
    logger.info("Starting coarse-graining process")
    cg_data = {'R':{},'z':{},'F':{},'cell':{}}
    n_frames = len(dump_data.keys())
    
    for index, time_step in enumerate(list(dump_data.keys())):
        logger.info(f"Processing frame {index+1}/{n_frames}")
        frame_data = dump_data[time_step]
        cg_group_site, cg_group_coord, cg_group_force, cg_box = cg_map(frame_data, system_data,target=args.target)
        
        if args.wrap:
            logger.debug("Wrapping coordinates")
            cg_group_coord = wrap_coordinates(cg_group_coord, cg_box)
        
        cg_data['R'][index] = cg_group_coord
        cg_data['F'][index] = cg_group_site
        cg_data['z'][index] = cg_group_force
        cg_data['cell'][index] = cg_box
    
    # Write output files
    logger.info(f"Writing output in {args.format} format")
    if args.format in ['xyz']:
        if args.separate_xyz:
            write_xyz_separate_frames(f'{args.output}', cg_data)
            logger.info(f"Written separate XYZ files: {args.output}_*.xyz")
        else:
            write_xyz_trajectory(f'{args.output}.xyz', cg_data)
            logger.info(f"Written XYZ trajectory file: {args.output}.xyz")
    
    elif args.format in ['data']:
        write_lammps_data(f'{args.output}.data', cg_data, frame_idx=args.frame)
        logger.info(f"Written LAMMPS data file: {args.output}.data (frame {args.frame})")
    
    elif args.format in ['npz']:
        npz_filename = f'{args.output}.npz'
        if save_cg_data_npz(npz_filename, cg_data):
            logger.info(f"Written NPZ file: {npz_filename}")
        else:
            logger.error("Failed to write NPZ file")
    
    logger.info("CGMap completed successfully")

if __name__ == "__main__":
    main()