This repository documents the changes I made to a script used within my research group to improve the performance and readability of the code being used. 

Originally, the code was performing the same tasks multiple times and had hard coded values that would change depending on what you wanted to examine. The first torder of business was to reorganize the code into functions and then add in an [argument parser](https://docs.python.org/3/library/argparse.html) to allow those values to be changed from the command-line. Additionally, the bulk of the code was using a nested looping structure to create spherical shells and determine the particles within those shells. The [modified](/old_shells.py) version highlights the changes made along with explanations around the sections changed.

All of the changes were then split into two files, [helpers.py](/helpers.py) for the helper functions and [dispersion.py](/dispersion.py) for the "main" function.
