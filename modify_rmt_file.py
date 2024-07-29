def modify_coordinates(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    with open(output_file, 'w') as f:
        for line in lines:
            if '-' in line:  # Check if the line contains a hyphen, indicating a range
                parts = line.split()
                chromosome = parts[0]  # Assuming chromosome information is at the beginning of the line
                for part in parts[1:]:
                    if '-' in part:  # Check if the part contains a hyphen, indicating a range
                        start, end = map(int, part.split('-'))
                        start_modified = start + 1
                        end_modified = end - 1
                        f.write(f"{chromosome} ({start_modified})-({end_modified}) None\n")
            else:
                f.write(line)  # Write other lines as is


# Input and output filenames
input_filename = "Arabidopsis_thaliana_genes.rmt"
output_filename = "modified_Arabidopsis_thaliana_genes.rmt"

# Modify coordinates
modify_coordinates(input_filename, output_filename)
