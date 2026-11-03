#!/usr/bin/env python3
# reorganize_trees.py
#
# Adapted from reorganize_trees.py by Anik Saha: the same output, but the gene trees are copied line by line
# instead of being read whole, so million-taxon replicates need little memory.
from pathlib import Path
import re
import sys
import shutil


def natural_key(p: Path):
    """Extract numeric part from path for sorting"""
    m = re.search(r'(\d+)', p.stem)
    return int(m.group(1)) if m else 0


def copy_cleaned(src: Path, dst: Path):
    """Copy the gene trees, one per line, without the _0_0 suffixes of the leaf labels."""
    with src.open(encoding="utf-8") as r, dst.open("w", encoding="utf-8") as w:
        for line in r:
            w.write(line.replace("_0_0", ""))


def reorganize_trees(root_dir: Path):
    """Reorganize tree files from numeric subdirs to R-prefixed subdirs"""
    if not root_dir.is_dir():
        print(f"Error: {root_dir} is not a directory")
        return False

    numeric_subdirs = [d for d in root_dir.iterdir() if d.is_dir() and d.name.isdigit()]
    if not numeric_subdirs:
        print(f"No numeric subdirectories found in {root_dir}")
        return False
    numeric_subdirs.sort(key=natural_key)

    processed_count = 0
    successfully_processed_dirs = []
    for subdir in numeric_subdirs:
        r_dir_path = root_dir / f"R{subdir.name}"
        gene_tree_file = subdir / "all_gt.tre"
        species_tree_file = subdir / "s_tree.trees"
        missing_files = [f.name for f in (gene_tree_file, species_tree_file) if not f.exists()]
        if missing_files:
            print(f"Skipping {subdir.name}: missing files {missing_files}")
            continue
        r_dir_path.mkdir(exist_ok=True)
        copy_cleaned(gene_tree_file, r_dir_path / "all_gt.tre")
        shutil.copy2(species_tree_file, r_dir_path / "s_tree.trees")
        print(f"Created R{subdir.name}/ with gene trees and species tree from {subdir.name}/")
        processed_count += 1
        successfully_processed_dirs.append(subdir)

    if processed_count == 0:
        print("No directories were processed")
        return False
    print(f"Processed {processed_count} replicates; removing the raw numbered folders")
    for subdir in successfully_processed_dirs:
        try:
            shutil.rmtree(subdir)
        except Exception as e:
            print(f"Warning: could not remove {subdir.name}/: {e}")
    return True


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("out")
    if not root.exists():
        print(f"Error: Directory {root} does not exist")
        return
    reorganize_trees(root)


if __name__ == "__main__":
    main()
