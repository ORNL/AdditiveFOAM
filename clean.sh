#!/bin/bash

# Function to attempt removing .nfs files, even if in use
remove_nfs_files() {
  find "$1" -type f -name '.nfs*' | while read -r nfs_file; do
    echo "Removing: $nfs_file"
    
    # Find the process holding the file and kill it
    lsof_output=$(lsof "$nfs_file" 2>/dev/null)
    if [[ -n "$lsof_output" ]]; then
      pid=$(echo "$lsof_output" | awk 'NR>1 {print $2}' | sort -u)
      if [[ -n "$pid" ]]; then
        echo "Killing process $pid holding file $nfs_file"
        kill -9 $pid
      fi
    fi

    # Attempt to remove the file again
    rm -f "$nfs_file"
  done
}

# Navigate through every directory from the current one
echo "Starting search and removal of .nfs files from: $(pwd)"
remove_nfs_files "$(pwd)"

# Traverse through all subdirectories
find . -type d | while read -r dir; do
  remove_nfs_files "$dir"
done

echo "Finished removing .nfs files."

