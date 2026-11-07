import argparse
import json
from .core import qc, normalize

p=argparse.ArgumentParser()
s=p.add_subparsers(dest='command',required=True)
q=s.add_parser('qc'); q.add_argument('vcf'); q.add_argument('--out',required=True); q.add_argument('--sample')
n=s.add_parser('normalize'); n.add_argument('vcf'); n.add_argument('--fasta',required=True); n.add_argument('--out',required=True)
a=p.parse_args()
result=qc(a.vcf,a.out,a.sample) if a.command=='qc' else normalize(a.vcf,a.fasta,a.out)
print(json.dumps(result,indent=2))
