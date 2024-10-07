import subprocess
import sys
import os

def run_executable(executable, isovalue, input_file, output_format, output_file, options):
    """Runs the dmr executable with the provided arguments and returns success status."""
    command = [executable, str(isovalue), input_file, output_format, output_file] + options
    try:
        subprocess.run(command, check=True)
        print(f"{executable} ran successfully and generated {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {executable}: {e}")
        return False

def run_ijkmeshinfo(mesh_file, options):
    """Runs ijkmeshinfo on a given mesh file with specific options and returns the output or an error message."""
    command = ['./ijkmeshinfo'] + options + [mesh_file]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        error_message = f"Error running ijkmeshinfo on {mesh_file}: {e}\n"
        error_message += f"Standard Output:\n{e.stdout}\nStandard Error:\n{e.stderr}"
        return error_message
    
def parse_off_mesh(filepath):
    """Parses an OFF mesh file and returns a list of vertices and faces."""
    vertices = []
    faces = []
    try:
        with open(filepath, 'r') as file:
            lines = file.readlines()
            if lines[0].strip() == 'OFF':
                counts = list(map(int, lines[1].split()))
                vertex_count = counts[0]
                face_count = counts[1]

                # Read vertices
                for i in range(2, 2 + vertex_count):
                    vertex = tuple(map(float, lines[i].split()))
                    vertices.append(vertex)

                # Read faces as sets (to handle unordered comparison)
                for i in range(2 + vertex_count, 2 + vertex_count + face_count):
                    face = tuple(map(int, lines[i].split()[1:]))  # Ignore the first number, which is face count
                    faces.append(frozenset(face))  # Use frozenset to ignore order
    except FileNotFoundError:
        print(f"File {filepath} not found.")
    return vertices, faces

def parse_ply_mesh(filepath):
    """Parses a PLY mesh file and returns a list of vertices and faces."""
    vertices = []
    faces = []
    header_ended = False
    vertex_count = 0
    face_count = 0
    vertex_started = False
    face_started = False

    try:
        with open(filepath, 'r') as file:
            for line in file:
                if not header_ended:
                    # Process header
                    if line.startswith("element vertex"):
                        vertex_count = int(line.split()[-1])
                    elif line.startswith("element face"):
                        face_count = int(line.split()[-1])
                    elif line.strip() == "end_header":
                        header_ended = True
                        vertex_started = True
                    continue

                if vertex_started and vertex_count > 0:
                    vertex = tuple(map(float, line.split()))
                    vertices.append(vertex)
                    vertex_count -= 1
                    if vertex_count == 0:
                        vertex_started = False
                        face_started = True
                    continue

                if face_started and face_count > 0:
                    face = tuple(map(int, line.split()[1:]))  # Ignore the first number, which is face count
                    faces.append(frozenset(face))  # Use frozenset to ignore order
                    face_count -= 1
    except FileNotFoundError:
        print(f"File {filepath} not found.")
    return vertices, faces

def compare_lists(list1, list2, element_name):
    """Compare two lists and return the elements that are in one list but not the other."""
    set1 = set(list1)
    set2 = set(list2)

    only_in_1 = set1 - set2
    only_in_2 = set2 - set1

    result = []
    result.append(f"\n{element_name} present in the first mesh but not in the second ({len(only_in_1)} items):")
    for elem in only_in_1:
        result.append(str(sorted(elem)))  # Sort for better readability

    result.append(f"\n{element_name} present in the second mesh but not in the first ({len(only_in_2)} items):")
    for elem in only_in_2:
        result.append(str(sorted(elem)))  # Sort for better readability

    return result

def compare_meshes(mesh1, mesh2, format):
    """Compares two mesh files in terms of vertices and faces and returns the comparison as a string."""
    comparison_results = []
    if format == 'off':
        vertices1, faces1 = parse_off_mesh(mesh1)
        vertices2, faces2 = parse_off_mesh(mesh2)
    elif format == 'ply':
        vertices1, faces1 = parse_ply_mesh(mesh1)
        vertices2, faces2 = parse_ply_mesh(mesh2)
    else:
        comparison_results.append("Unsupported format")
        return comparison_results

    # Compare vertices
    comparison_results.extend(compare_lists(vertices1, vertices2, "Vertices"))

    # Compare faces
    comparison_results.extend(compare_lists(faces1, faces2, "Faces"))

    return comparison_results

def write_to_markdown(output_file, comparison_data):
    """Writes the comparison data to a markdown file."""
    with open(output_file, 'w') as md_file:
        md_file.write("# Mesh Comparison Results\n")
        for section in comparison_data:
            md_file.write("\n## Input File: " + section['input_file'] + "\n")
            md_file.write("\n".join(section['results']))
            md_file.write("\n\n---\n")

if __name__ == "__main__":
    # Check for correct usage
    if len(sys.argv) < 8:
        print("Usage: python compare_dmr.py <output format> <executable1> <executable2> <isovalue> [<input file 1> <input file 2> ...] <options...>")
        sys.exit(1)

    # Parse command-line arguments
    output_format = sys.argv[1]
    executable1 = sys.argv[2]
    executable2 = sys.argv[3]
    isovalue = float(sys.argv[4])
    
    # Extract input files (surrounded by [])
    input_files_index = sys.argv.index("[") + 1
    end_input_files_index = sys.argv.index("]")
    input_files = sys.argv[input_files_index:end_input_files_index]
    
    # Extract options after input files
    options = sys.argv[end_input_files_index + 1:]

    # Store comparison results
    comparison_data = []

    # Process each input file
    for input_file in input_files:
        base_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file1 = f"{base_filename}_{executable1[2:]}1.{output_format}"
        output_file2 = f"{base_filename}_{executable2[2:]}2.{output_format}"

        # Run the two executables
        if run_executable(executable1, isovalue, input_file, output_format, output_file1, options) and \
           run_executable(executable2, isovalue, input_file, output_format, output_file2, options):
            # Compare the outputs if both executables ran successfully
            results = compare_meshes(output_file1, output_file2, output_format)

            # Run ijkmeshinfo on both output files and collect the results
            meshinfo1 = run_ijkmeshinfo(output_file1, ['-manifold'])
            meshinfo2 = run_ijkmeshinfo(output_file2, ['-manifold'])

            # Add ijkmeshinfo comparison to the results
            results.append("\n\n# ijkmeshinfo comparison:\n")
            results.append(f"**Results for {output_file1}:**\n{meshinfo1}")
            results.append(f"**Results for {output_file2}:**\n{meshinfo2}")

            comparison_data.append({
                'input_file': input_file,
                'results': results
            })

    # Write all comparison results to a markdown file
    write_to_markdown("mesh_comparison_results.md", comparison_data)

    print(f"Comparison results have been written to mesh_comparison_results.md")
