import os

os.system("rm *~ *.o")

for file in os.listdir(os.getcwd()):
    if '.cpp' in file or '.c' in file:
        os.system(f"g++ -c -fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations -o {file.split('.')[0]}.o {file}")
        print    (f"g++ -c -fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations -o {file.split('.')[0]}.o {file}")

print('=' * 60)

linker = "g++ -lblas -fopenmp -o gkmPWMlasso3 "
for file in os.listdir(os.getcwd()):
    if '.o' == file[-2:]:
        linker += f" {file}"
linker += " libopenblas_haswellp-r0.3.20.a"
print(linker)
os.system(linker)
  
print()

