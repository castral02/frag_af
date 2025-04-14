#!/usr/bin/env python3
# Creates fasta and text files required for AlphaPulldown
##

import os
import requests
from absl import flags, app, logging

# Define command-line flags
FLAGS = flags.FLAGS
flags.DEFINE_string('directory', None, 'Directory to create and put the file in')
flags.DEFINE_string('uniprot_id', None, 'UniProt ID of the Protein of Interest')


def main(argv):
    # Access the flag values
    directory = FLAGS.directory
    uniprot_id = FLAGS.uniprot_id

    # Grabbing contents for fasta/text files
    sequence = get_protein_sequence(uniprot_id)
    if sequence is None:
        logging.error(f"Failed to fetch protein sequence for UniProt ID {uniprot_id}.")
        return
    tiled_sequence = tile_sequence(sequence)
    title = fragment_title(tiled_sequence)  # This goes into a text file

    # Creating directory
    create_directory(directory)

    # Making fasta files
    tiled_fasta(directory, filename='tiled.fasta', fragments=tiled_sequence)
    full_fasta(directory, filename='full.fasta', uniprot_id=uniprot_id, sequence=sequence)

    # Making text files
    text_file(directory, filename='tiled.txt', content=title)
    text_file(directory, filename='full.txt', content=[uniprot_id])
    combine_fasta(directory, file1='tiled.fasta', file2='full.fasta', output_file='combined.fasta')


# Grabbing and Fragmenting protein sequence
def get_protein_sequence(uniprot_id):
    """Fetch the protein sequence from UniProt using its UniProt ID.

    Parameters:
        uniprot_id (string): UniProt
    Returns:
        sequence (string): Protein sequence"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        lines = response.text.split("\n")
        sequence = "".join(lines[1:])
        return sequence
    else:
        print(f"Failed to fetch protein sequence for UniProt ID {uniprot_id}. Status code: {response.status_code}")
        return None


def tile_sequence(sequence, chunk_size=60, window_size=10):
    """Chunk sequence into 60 amino acid chunks with 10 amino acid sliding window
    Parameters:
        sequence (string): amino acid sequence
    Returns:
        chunks (list): list of amino acid chunks"""
    chunks = []
    for i in range(0, len(sequence), window_size):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) == chunk_size:
            chunks.append(chunk)
    return chunks


def fragment_title(chunks):
    """Creates a list of all the fragment names
    Parameters:
        chunks (list): list of all the chunks of amino acid sequences
    Returns:
        title (list): list of all the fragment names"""
    title = []
    for i in range(len(chunks)):
        title.append(f"Fragment_{i + 1}")
    return title


# Creation of all the files necessary
def tiled_fasta(directory, filename, fragments):
    file_path = os.path.join(directory, filename)
    with open(file_path, "w") as file:
        for i, fragment in enumerate(fragments):
            file.write(f">Fragment_{i + 1}\n")
            file.write(fragment + "\n")


def full_fasta(directory, filename, uniprot_id, sequence):
    file_path = os.path.join(directory, filename)
    with open(file_path, "w") as file:
        file.write(f">{uniprot_id}\n")
        file.write(sequence + "\n")


def combine_fasta(directory, file1, file2, output_file):
    file1_path = os.path.join(directory, file1)
    file2_path = os.path.join(directory, file2)
    output_file_path = os.path.join(directory, output_file)

    with open(output_file_path, 'w') as outfile:
        # Read and write the contents of the first file
        with open(file1_path, 'r') as infile1:
            contents1 = infile1.read()
            outfile.write(contents1)
            # Ensure there's a newline between files if not already present
            if not contents1.endswith('\n'):
                outfile.write('\n')

        # Read and write the contents of the second file
        with open(file2_path, 'r') as infile2:
            contents2 = infile2.read()
            outfile.write(contents2)


def text_file(directory, filename, content):
    file_path = os.path.join(directory, filename)
    with open(file_path, "w") as file:
        if isinstance(content, list):
            content = "\n".join(content)
        file.write(content + "\n")


def create_directory(path):
    """Create a new directory if it does not exist."""
    try:
        os.makedirs(path, exist_ok=True)
        logging.info(f"Directory '{path}' created successfully.")
    except Exception as e:
        logging.error(f"Error creating directory '{path}': {e}")
        raise


if __name__ == '__main__':
    # Mandatory flags
    flags.mark_flag_as_required('directory')
    flags.mark_flag_as_required('uniprot_id')
    app.run(main)
