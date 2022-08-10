This is the README file for the KhovanovSteenrodSquare program.

**********************************

Here is a step by step guideline on using this code to compute the desired 
Steenrod squares on Khovanov homology for a link. Also, note that this does 
not work if the link diagram has an unknot or Hopf link component. This is 
because it relies on code pulled from the KhovanovSteenrod codebase written 
by Robert Lipshitz and Sucharit Sarkar 
(https://github.com/sucharit/KhovanovSteenrod).




-Install Sage. (Visit http://www.sagemath.org/)

-Open a Sage terminal, go to the directory where you have downloaded all
the files of KhovanovSteenrodSquare.

-Update runKhovanovSteenrodSquare.sage with the link of interest, defined either 
via PD code or DT code as a list of lists. See the notes below on defining the link. 
Specify the particular Steenrod square  that we want to compute. To do this for 
the map sq^i: H^(p,q) -> H^(p+i,q), we  would use the inputs: i, p, q.

-Compute sq^i for this link by running
    > sage runKhovanovSteenrodSquare.sage
from the Sage terminal.


Notes on defining the link:
- The indexing for the PD code or DT code must start with a 1, not a 0.
- To determine the PD code for a link diagram, follow the instructions here: https://knotinfo.math.indiana.edu/descriptions/pd_notation.html .
- The way that the Sage documentation describes the PD code is NOT correct. Do not follow that.
- Because we rely on some of Lipshitz and Sarkar's function, which do not work for links with unknot components or Hopf link components, this repo does not work with links with either of those as components either.
