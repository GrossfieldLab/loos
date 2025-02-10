#!/usr/bin/env python3

from pathlib import Path
import shutil
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage- omg_copy [example|doc] destination_path")
        sys.exit(1)

    target = sys.argv[1]
    destination = Path(sys.argv[2])
    path_to_lib = Path(__file__).parent
    
    if target == "example":
        target_path = path_to_lib / "example"
    elif target == "doc":
        target_path = path_to_lib / "doc"
    else:
        raise RuntimeError(f"Unknown copy target {target} (should be either example or doc)")

    shutil.copytree(target_path, destination / target)
    
    
if __name__ == "__main__":
    main()