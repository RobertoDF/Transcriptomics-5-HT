import nbformat

# Load your notebook
notebook_path = '/home/roberto/Github/Transcriptomics-5-HT/Figure 2.ipynb'  # Replace with your notebook path
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# Modify all cells
for cell in nb.cells:
    if 'metadata' not in cell:
        cell['metadata'] = {}
    cell['metadata']['slideshow'] = {'slide_type': 'skip'}

# Save the modified notebook
with open(notebook_path, 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)

print("All cells have been set to slide type 'skip'")
