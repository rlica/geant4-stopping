import sys

def parse_fresco(infile, outfile):
    with open(infile, 'r') as f:
        lines = f.readlines()
        
    start_idx = -1
    for i, line in enumerate(lines):
        if 'Partition=  2 Excit=  1' in line:
            start_idx = i + 1
            break
            
    if start_idx == -1:
        print('Could not find Partition 2 Excit 1 in fresco file')
        return
        
    # find where data starts
    while start_idx < len(lines) and not lines[start_idx].strip()[0].isdigit() and not lines[start_idx].strip().startswith('.'):
        start_idx += 1
        
    with open(outfile, 'w') as out:
        out.write('# Extracted from angular_distribution_fresco.txt (Partition 2, Excit 1)\n')
        out.write('# Theta_CM(deg)  dSigma/dOmega(mb/sr)\n')
        for line in lines[start_idx:]:
            if line.strip() == '&' or line.startswith('@') or line.startswith('#'):
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                # Fresnel output can have 'E' for exponent
                out.write(f'{parts[0]}  {parts[1]}\n')

parse_fresco('angular_distribution_fresco.txt', 'angular_distribution.txt')
