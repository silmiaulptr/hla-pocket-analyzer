import pymol
from pymol import cmd
import csv
import os

# Initialize PyMOL
pymol.finish_launching()

# HLA Pocket mapping residues (Optimized to avoid overlap)
pocket_map = {
    "A": [5, 59, 159, 163, 167, 171],
    "B": [7, 9, 24, 34, 45, 63, 66, 67, 70, 99],
    "C": [9, 70, 73, 74],
    "D": [99, 155, 156, 160],
    "E": [97, 114, 147, 152],
    "F": [77, 80, 81, 84, 95, 123, 143, 146]
}


def identify_pocket(chain, resi):
    try:
        resi = int(resi)
        if chain != "A":
            return ""
        for pocket, resi_list in pocket_map.items():
            if resi in resi_list:
                return pocket
        return ""
    except:
        return ""

def get_full_atom_path(selection, index):
    """Get the full PyMOL path for an atom by index"""
    stored.atom_path = ""
    cmd.iterate(f"{selection} and index {index}", 
                "stored.atom_path = f'/{model}//{chain}/{resn}\\\''{resi}/{name}'")
    return stored.atom_path

# BASIC VISUALIZATION
cmd.hide("everything")
cmd.show("cartoon", "receptor")
cmd.color("lightblue", "receptor")
cmd.set("cartoon_transparency", 0.3, "receptor")
cmd.show("sticks", "peptide")
cmd.color("yellow", "peptide")

cmd.bg_color("white")
cmd.set("ray_opaque_background", "off")

# HYDROGEN BOND ANALYSIS
# Select H-bond candidates (within 3.5 Å distance)
cmd.select("hbonds", "(receptor within 3.5 of peptide) and (peptide within 3.5 of receptor)")

# Visualize H-bonds as dashed lines
try:
    # Delete existing hbond_lines if they exist
    cmd.delete('hbond_lines')
    # Create new hydrogen bond visualization
    cmd.distance("hbond_lines", "receptor", "peptide", mode=2, cutoff=3.5)
    cmd.set("dash_width", 2)
    cmd.set("dash_gap", 0.3)
    cmd.set("dash_color", "magenta")
    # Show distance labels with 2 decimal places
    cmd.set("label_digits", 2, "hbond_lines")
    cmd.set("label_distance", 1, "hbond_lines")
    
    # Only show hydrogen bonds with meaningful values (between 2.0 and 3.5 Å)
    cmd.hide("labels", "hbond_lines and distance > 3.5")
    cmd.hide("dashes", "hbond_lines and distance > 3.5")
except pymol.CmdException as e:
    print("Warning: Could not create hydrogen bond visualization:", str(e))

# HLA POCKET DEFINITIONS
cmd.select("pocketA", "chain A and resi 5+59+159+163+167+171")
cmd.select("pocketB", "chain A and resi 7+9+24+34+45+63+66+67+70+99")
cmd.select("pocketC", "chain A and resi 9+70+73+74")
cmd.select("pocketD", "chain A and resi 99+155+156+160")
cmd.select("pocketE", "chain A and resi 97+114+147+152")
cmd.select("pocketF", "chain A and resi 77+80+81+84+95+123+143+146")

# Color and show pocket sticks for qualitative analysis
cmd.color("red", "pocketA")
cmd.color("orange", "pocketB")
cmd.color("green", "pocketC")
cmd.color("blue", "pocketD")
cmd.color("cyan", "pocketE")
cmd.color("purple", "pocketF")
cmd.show("sticks", "pocketA or pocketB or pocketC or pocketD or pocketE or pocketF")

# DATA EXPORT & POCKET IDENTIFICATION
from pymol import cmd
import csv
import os

# Create directory if it doesn't exist
try:
    documents_folder = os.path.expanduser("~\\Documents")
    output_folder = os.path.join(documents_folder, "PyMOL_Results")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_path = os.path.join(output_folder, "hbonds_info.csv")
    print(f"CSV will be saved to: {output_path}")
except Exception as e:
    # Fallback to the script's directory if Documents folder is not accessible
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "hbonds_info.csv")
    print(f"Could not access Documents folder. CSV will be saved to: {output_path}")

# Pocket mapping residues
pocket_map = {
    "A": [5, 59, 159, 163, 167, 171],
    "B": [7, 9, 24, 34, 45, 63, 66, 67, 70, 99],
    "C": [9, 70, 73, 74],
    "D": [99, 155, 156, 160],
    "E": [97, 114, 147, 152],
    "F": [77, 80, 81, 84, 95, 123, 143, 146]
}

def identify_pocket(chain, resi):
    try:
        resi = int(resi)
    except:
        return ""
    if chain != "A":
        return ""
    for pocket, resi_list in pocket_map.items():
        if resi in resi_list:
            return pocket
    return ""

# Check if required selections exist
receptor_atoms = cmd.count_atoms("receptor")
peptide_atoms = cmd.count_atoms("peptide")

print(f"DEBUG: Found {receptor_atoms} atoms in 'receptor' selection")
print(f"DEBUG: Found {peptide_atoms} atoms in 'peptide' selection")

# Check minimum distance between receptor and peptide
min_dist = cmd.distance("temp_dist", "receptor", "peptide", mode=0)
cmd.delete("temp_dist")
print(f"DEBUG: Minimum distance between receptor and peptide: {min_dist:.2f} Å")

if not receptor_atoms:
    print("Error: 'receptor' selection does not exist. Please load and select the receptor first.")
    import sys
    sys.exit(1)

if not peptide_atoms:
    print("Error: 'peptide' selection does not exist. Please load and select the peptide first.")
    import sys
    sys.exit(1)

# Only proceed with H-bond analysis if selections exist
print("Found receptor and peptide selections, proceeding with H-bond analysis...")

# Create a minimal test CSV to verify we can write files
try:
    test_csv = os.path.join(os.path.dirname(output_path), "test_minimal.csv")
    with open(test_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Test", "CSV", "File"])
        writer.writerow(["Receptor", str(receptor_atoms), "atoms"])
        writer.writerow(["Peptide", str(peptide_atoms), "atoms"])
    print(f"DEBUG: Successfully created test CSV: {test_csv}")
except Exception as e:
    print(f"DEBUG: Error creating test CSV: {str(e)}")

# Initialize stored variables
cmd.do("stored.atom1 = []")
cmd.do("stored.atom2 = []")

hbonds_data = []

# Use a basic approach to visualize hydrogen bonds and extract their data
print("Creating hydrogen bond visualization and extracting data...")

try:
    # Step 1: Delete any existing hbond objects to start fresh
    try:
        cmd.delete("hbond_lines")
    except:
        pass
        
    # Step 2: Create a temporary selection of candidates for hydrogen bonds
    print("Finding hydrogen bond candidates...")
    try:
        # Select close atoms
        cmd.select("_close_atoms", "(receptor within 3.5 of peptide) or (peptide within 3.5 of receptor)")
        # Further refine to only include N and O atoms
        cmd.select("_hb_candidates", "_close_atoms and (elem N or elem O)")
    except Exception as e:
        print(f"Error creating selection: {str(e)}")
        cmd.select("_hb_candidates", "all")
        
    # Step 3: Find receptor atoms that could form hydrogen bonds
    stored.receptor_atoms = []
    try:
        print("Finding receptor atoms that could form hydrogen bonds...")
        cmd.iterate("_hb_candidates and receptor", 
                    "stored.receptor_atoms.append((model, chain, resi, resn, name, index))")
        print(f"Found {len(stored.receptor_atoms)} receptor atoms")
    except Exception as e:
        print(f"Error finding receptor atoms: {str(e)}")
        stored.receptor_atoms = []
        
    # Step 4: Find peptide atoms that could form hydrogen bonds
    stored.peptide_atoms = []
    try:
        print("Finding peptide atoms that could form hydrogen bonds...")
        cmd.iterate("_hb_candidates and peptide", 
                    "stored.peptide_atoms.append((model, chain, resi, resn, name, index))")
        print(f"Found {len(stored.peptide_atoms)} peptide atoms")
    except Exception as e:
        print(f"Error finding peptide atoms: {str(e)}")
        stored.peptide_atoms = []
        
    # Step 5: Check each pair for hydrogen bonds
    print("Checking atom pairs for hydrogen bonds...")
    
    # Create new visualization for hydrogen bonds
    try:
        cmd.distance("hbond_lines", "receptor", "peptide", cutoff=3.5, mode=2)
        cmd.set("dash_width", 2, "hbond_lines")
        cmd.set("dash_gap", 0.3, "hbond_lines")
        cmd.set("dash_color", "magenta", "hbond_lines")
    except Exception as e:
        print(f"Warning: Could not create hydrogen bond visualization: {str(e)}")
    
    # Direct approach - check each pair of atoms for potential hydrogen bonds
    hbond_count = 0
    
    # For debugging, we'll add all pairs we check to a list
    all_checked_pairs = []
    
    # Check each receptor atom against each peptide atom
    for r_atom in stored.receptor_atoms:
        r_model, r_chain, r_resi, r_resn, r_name, r_index = r_atom
        # Create a more specific selection to avoid ambiguity
        r_sel = f"receptor and chain {r_chain} and resi {r_resi} and name {r_name}"
        
        for p_atom in stored.peptide_atoms:
            p_model, p_chain, p_resi, p_resn, p_name, p_index = p_atom
            # Create a more specific selection to avoid ambiguity
            p_sel = f"peptide and chain {p_chain} and resi {p_resi} and name {p_name}"
            
            try:
                # Measure distance between atoms
                distance = cmd.get_distance(r_sel, p_sel)
                pair_info = f"{r_resn}{r_resi}.{r_name} - {p_resn}{p_resi}.{p_name}: {distance:.2f}Å"
                all_checked_pairs.append(pair_info)
                
                # Check if within hydrogen bonding range
                if 2.5 <= distance <= 3.5:
                    hbond_count += 1
                    hb_name = f"hbond_{hbond_count}"
                    
                    try:
                        # Create a visual representation of this specific bond
                        cmd.distance(hb_name, r_sel, p_sel, label=0)
                        cmd.set("dash_width", 2, hb_name)
                        cmd.set("dash_gap", 0.3, hb_name)
                        cmd.set("dash_color", "magenta", hb_name)
                    except:
                        # If visualization fails, continue with data collection
                        pass
                    
                    # Manually construct the paths using the data we already have
                    # We already have all the components from the atom selection
                    # This avoids any string formatting issues in PyMOL's iterate command
                    
                    # For receptor
                    r_path = f"/{r_model}//{r_chain}/{r_resn}'{r_resi}/{r_name}"
                    
                    # For peptide
                    p_path = f"/{p_model}//{p_chain}/{p_resn}'{p_resi}/{p_name}"
                    
                    # Print the paths for debugging
                    print(f"Receptor path: {r_path}")
                    print(f"Peptide path: {p_path}")
                    
                    # Identify pocket
                    pocket = identify_pocket(r_chain, r_resi)
                    
                    # Create record of this hydrogen bond
                    hbond_data = {
                        "receptor_path": r_path,
                        "peptide_path": p_path,
                        "receptor_chain": r_chain,
                        "receptor_resi": r_resi,
                        "receptor_resn": r_resn,
                        "receptor_atom": r_name,
                        "peptide_chain": p_chain,
                        "peptide_resi": p_resi,
                        "peptide_resn": p_resn,
                        "peptide_atom": p_name,
                        "distance": distance,
                        "pocket": pocket
                    }
                    
                    hbonds_data.append(hbond_data)
                    
                    # Print for debugging
                    print(f"Bond {hbond_count}: {r_path} -- {p_path} ({distance:.2f} Å) in {pocket if pocket else 'no pocket'}")
            except Exception as e:
                print(f"Warning: Could not measure distance between {r_sel} and {p_sel}: {str(e)}")
                continue
    
    # Debug info
    print(f"Checked {len(all_checked_pairs)} atom pairs in total")
    print(f"Sample of pairs checked (first 5): {all_checked_pairs[:5]}")
    print(f"Found {len(hbonds_data)} hydrogen bonds between receptor and peptide")
    
    # Clean up temporary selections
    try:
        cmd.delete("_close_atoms")
        cmd.delete("_hb_candidates")
    except:
        pass


    
    # No need for additional debug output - we already printed each bond as we found it
    
    # Write results to CSV with receptor and peptide paths
    try:
        print(f"\nAttempting to write hydrogen bond paths to CSV file: {output_path}")
        with open(output_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Receptor Path", "Peptide Path", "Distance (Å)", "Receptor Chain", "Receptor Residue", 
                           "Receptor Name", "Receptor Atom", "Peptide Chain", "Peptide Residue", 
                           "Peptide Name", "Peptide Atom", "Pocket"])
            
            if hbonds_data:
                # Filter out bonds with distances > 3.5 Å (same as those hidden in visualization)
                valid_hbonds = [hb for hb in hbonds_data if hb['distance'] <= 3.5]
                
                if valid_hbonds:
                    for hb in valid_hbonds:
                        writer.writerow([hb["receptor_path"], hb["peptide_path"], f"{round(hb['distance'], 1):.1f}", 
                                        hb["receptor_chain"], hb["receptor_resi"], hb["receptor_resn"], hb["receptor_atom"],
                                        hb["peptide_chain"], hb["peptide_resi"], hb["peptide_resn"], hb["peptide_atom"],
                                        hb["pocket"]])
                else:
                    writer.writerow(["NO VALID HYDROGEN BONDS FOUND (≤ 3.5 Å)", "", "", "", "", "", "", "", "", "", "", ""])
                
                # Print summary of filtered bonds
                print(f"Found {len(hbonds_data)} total hydrogen bonds, {len(valid_hbonds)} valid bonds (≤ 3.5 Å)")
                print(f"Excluded {len(hbonds_data) - len(valid_hbonds)} bonds with distances > 3.5 Å from CSV output")
                
            else:
                writer.writerow(["NO HYDROGEN BONDS FOUND", "", "", "", "", "", "", "", "", "", "", ""])
        
        # Verify file was created
        if os.path.exists(output_path):
            print(f"SUCCESS: Hydrogen bond paths saved to {output_path}")
            print(f"File size: {os.path.getsize(output_path)} bytes")
        else:
            print(f"ERROR: File creation failed. {output_path} does not exist.")
    except Exception as e:
        print(f"ERROR writing CSV file: {str(e)}")

except pymol.CmdException as e:
    print(f"Error analyzing hydrogen bonds: {str(e)}")
except Exception as e:
    print(f"Unexpected error: {str(e)}")
finally:
    # Clean up temporary selections
    cmd.delete("_hb_candidates")
