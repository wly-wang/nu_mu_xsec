#!/bin/bash

# Prompt for the input file name
read -p "Enter the name of the input .txt file: " input_file

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file not found!"
    exit 1
fi

# Output file name
output_file="summary_output.txt"

# Initialize the row counter and running sum
row_counter=0
running_sum=0

# Clear the output file if it exists or create a new one
> "$output_file"

# Read each line from the input file
while IFS= read -r line
do
    # Increment the row counter
    row_counter=$((row_counter + 1))    

    # Create a temporary file with the current line
    echo "$line" > temp.txt

    # Run the Python command on temp.txt and capture output in log.txt
    python2.7 getDataInfo.py -v3 --format-numi --prescale --run-subrun-list ./temp.txt &> log.txt

    # Extract specific content from log.txt
    pot_val=$(awk '/tortgt_wcut/{getline; print $8}' log.txt)
    
    # Update the running sum
    running_sum=$(echo "$running_sum + $pot_val" | bc)

    # Create the output file name using the run and subrun
    run=$(echo $line | awk '{print $1}')
    subrun=$(echo $line | awk '{print $2}')

    # Write the row to the output file
    echo -e "$row_counter\t$run\t$subrun\t$pot_val\t$running_sum" >> "$output_file"

    # Clean up temporary files
    rm temp.txt log.txt

done < "$input_file"

echo "All events have been processed and results are saved in '$output_file'."

