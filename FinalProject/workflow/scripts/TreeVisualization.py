# Tree visualisation script
import toytree  
import toyplot.pdf

# load a tree from a newick string
tre = toytree.tree(snakemake.input[0])

# root and draw the tree
if tre.treenode.is_root():
    rtre = tre.copy()  # Create a copy of the tree since it's already rooted
else:
    rtre = tre.root(wildcard="*")  # Use wildcard selector to automatically select a root

rtre.draw(tip_labels_align=True);

# draw a plot and store the Canvas object to a variable
canvas, axes, mark = rtre.draw(width=400, height=300);

# save tree to pdf 
toyplot.pdf.render(canvas, snakemake.output[0])
