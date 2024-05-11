To make the folders more conveinent to work with, I will rename the folders to the following format:

Sample_1-a01va -> Sample_01-a01va
Sample_2-a01vb -> Sample_02-a01vb
...


```bash
for dir in Sample_*; do
  # Check if it's a directory
  if [[ -d $dir ]]; then
    # Extract the number part
    number=$(echo "$dir" | sed -n 's/Sample_\([0-9]*\)-.*/\1/p')

    # Skip if the number already has two digits
    if (( number >= 10 )); then
      continue
    fi

    # Pad the number with zero
    padded_number=$(printf "%02d" $number)
    
    # Replace the old number with the padded number
    new_dir=$(echo "$dir" | sed "s/Sample_$number-/Sample_$padded_number-/")
    
    # Rename the folder
    mv "$dir" "$new_dir" 
    echo "$dir -> $new_dir" 
  fi
done
```