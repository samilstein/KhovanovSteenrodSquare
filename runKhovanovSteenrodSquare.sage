outp='output.sage'

load KhovanovSteenrodSquare.sage
from datetime import datetime

writing = open(outp,'w')

print(f"Starting at {datetime.now()}")
# knot1b = [[3, 1, 4, 6], [1, 5, 2, 4], [5, 3, 6, 2]]
# knot2b = [[elem + max(flatten(knot1b)) for elem in lst] for lst in knot1b]
# knot3b = [[elem + max(flatten(knot2b)) for elem in lst] for lst in knot1b]
# knot4b = [[elem + max(flatten(knot3b)) for elem in lst] for lst in knot1b]


# compute_khovanov_complex(knot1b+knot2b+knot3b+knot4b, "PD_code", mirror_image = False)
compute_khovanov_complex([[1,2,4,3],[3,4,6,5],[5,6,2,1]], "PD_code")

# The inputs of compute_sqi are (i,p,q), where i is the index of the sq^i map, 
# p is the cohomology degree of the domain Khovanov homology group, and q is the quantum degree for this group.
sqi=KhovanovComplex.compute_sqi(1,-3,-7)
print(f"Ending at {datetime.now()}")
writing.write('\n sq1: H^(-3.-7)(4 trefoils) -> H^(-2,-7) = '+repr(sqi)+'\n')
writing.flush()


writing.close()

#2022-07-19 21:44:17