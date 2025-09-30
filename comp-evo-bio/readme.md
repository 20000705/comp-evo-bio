# Assignment 1 – Newick Tree Scaling

## Development Process
I created the code by iteratively experimenting and learning with the help of an AI assistant.  
I asked a few targeted questions to clarify concepts and get started with code structure. For example:

1. **Question:** *"Could you tell me how to read Newick tree in python?"*  
   - Helped me start the assignment.

2. **Question:** *"Can you give me starter Python code to scale branch lengths in a Newick string?"*  
   - Provided me with a regex-based function for branch length scaling.

3. **Question:** *"How do I test this in Jupyter Notebook without argparse failing?"*  
   - Learned how to simulate command-line arguments when testing interactively.

4. **Question:** *"Is my scaled output correct if I scale 0.1 → 0.05, 0.2 → 0.1, etc.?"*  
   - Confirmed that the scaling logic worked correctly.

I then adapted the starter code into a full command-line script with proper error handling and file input/output.

---

## Testing and Verification
- I tested the script on a small example Newick tree: ((A:0.1,B:0.2,C:0.3):0.4,(E:0.5,F:0.6):0.7)H:0.8; with `-s 0.5`, the branch lengths scaled exactly as expected.

- I ran the script on the provided `input.tre`.  
The scaled tree (`scaled.tre`) contained all branch lengths halved when using `-s 0.5`.

- I visualized both trees in **Figtree** and exported them as PDF/PNG.  
The scaled tree displayed the same topology with proportionally shorter branches.

---

## Reflections
The exercise went well:

- The regex approach turned out to be simple and effective.  
- Using AI was helpful for clarifying tricky parts quickly (e.g., regex usage, Jupyter vs CLI testing).  
- Debugging and testing on both small and large trees gave me confidence in the correctness of the script.  
- Visualizing in Figtree made it easy to confirm the results.

Overall, I understand both the **Newick format** and practical **Python scripting for phylogenetics**.

---

## How to Run
Run the script with:

```bash
python assignment1.py -i input.tre -s 0.5 -o scaled.tre

