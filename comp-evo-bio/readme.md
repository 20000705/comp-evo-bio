# Assignment 1 – Newick Tree Scaling

## Overview
This assignment required implementing a stand-alone Python program (`assignment1.py`) that:

- Reads a Newick tree from `input.tre`.
- Scales all branch lengths by a user-provided factor.
- Writes the scaled tree to `scaled.tre`.
- Prints both the original and scaled trees to the console.

The script was written using only Python standard libraries (and `re` for regex).  
It successfully handles multifurcations, missing branch lengths, and invalid input cases.

---

## Development Process
I created the code by iteratively experimenting and learning with the help of an AI assistant.  
I asked a few targeted questions to clarify concepts and get started with code structure. For example:

1. **Question:** *"Could you explain what the assignment requires?"*  
   - Helped me understand the input/output structure and constraints.

2. **Question:** *"Can you give me starter Python code to scale branch lengths in a Newick string?"*  
   - Provided me with a regex-based function for branch length scaling.

3. **Question:** *"How do I test this in Jupyter Notebook without argparse failing?"*  
   - Learned how to simulate command-line arguments when testing interactively.

4. **Question:** *"Is my scaled output correct if I scale 0.1 → 0.05, 0.2 → 0.1, etc.?"*  
   - Confirmed that the scaling logic worked correctly.

I then adapted the starter code into a full command-line script with proper error handling and file input/output.

---

## Testing and Verification
- I tested the script on a small example Newick tree:
