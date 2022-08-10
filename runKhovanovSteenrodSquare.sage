outp='output.sage'

load KhovanovSteenrodSquare.sage
from datetime import datetime

writing = open(outp,'w')

print(f"Starting at {datetime.now()}")

# Step 0: Manually specify the inputs.
link = [[1,2,4,3],[3,4,6,5],[5,6,2,1]] #This is the PD code or DT code for the oriented link of interest.
input_format = "PD_code" #This must be either PD_code or DT_code.
i = 1 #index of the Steenrod square to compute
p = -3 #cohomological degree of the domain of sq^i
q = -7 #quantum degree of the domain and range of sq^i

# Step 1: Set up the Khovanov complex
compute_khovanov_complex(link, input_format)

# Step 2: Compute sq^i : H^(p,q)(link; Z/2) -> H^(p+i,q)(link; Z/2)
sqi=KhovanovComplex.compute_sqi(i,p,q)

print(f"Ending at {datetime.now()}")
writing.write('\n sq1: H^(-3.-7)(trefoil) -> H^(-2,-7) = '+repr(sqi)+'\n')
writing.flush()
writing.close()
