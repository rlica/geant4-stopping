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
        out.write(f'# Extracted from {infile} (Partition 2, Excit 1)\n')
        out.write('# Theta_CM(deg)  dSigma/dOmega(mb/sr)\n')
        for line in lines[start_idx:]:
            if line.strip() == '&' or line.startswith('@') or line.startswith('#'):
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                # Fresnel output can have 'E' for exponent
                out.write(f'{parts[0]}  {parts[1]}\n')

if __name__ == '__main__':
    # Usage: python3 extract_fresco.py [<raw_fresco_output> [<angular_distribution>]]
    # Defaults preserved from the original single-reaction layout.
    infile = sys.argv[1] if len(sys.argv) > 1 else 'angular_distribution_fresco.txt'
    outfile = sys.argv[2] if len(sys.argv) > 2 else 'angular_distribution.txt'
    parse_fresco(infile, outfile)
