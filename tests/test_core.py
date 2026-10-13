import json
import pytest
from mva_explorer.core import allele_balance, genotype_type, audit_evidence, qc

def test_ab():
    assert allele_balance((10,5,5),2)==.25
    assert allele_balance((0,0),1) is None
    assert allele_balance((None,2),1) is None
    assert allele_balance((2,3),3) is None

def test_genotypes():
    assert genotype_type((1,2),2)=='heterozygous'
    assert genotype_type((1,1),1)=='homozygous_alt'
    assert genotype_type((None,1),1)=='missing_or_partial'
    assert genotype_type((0,0),1)=='not_carried'
    assert genotype_type((1,),1)=='hemizygous_or_haploid'

def test_evidence():
    rows=[dict(variant_key='synthetic',disease_id='test',dependency_group='pop',source='fixture',observation_id='x',criterion='PM2'),dict(variant_key='synthetic',disease_id='test',dependency_group='pop',source='fixture',observation_id='x',criterion='PP5')]
    a=audit_evidence(rows)
    assert 'absence_not_established' in a[0]['audit_flags']
    assert 'dependent_evidence' in a[1]['audit_flags']
    assert 'deprecated_criterion' in a[1]['audit_flags']

def test_pvs1_not_automatic():
    assert 'PVS1_lof_mechanism_unresolved' in audit_evidence([{'criterion':'PVS1'}])[0]['audit_flags']

def test_vcf(tmp_path):
    p=tmp_path/'synthetic.vcf'
    p.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n1\t10\t.\tA\tC,G\t30\tPASS\t.\tGT:AD\t1/2:10,5,5\n1\t20\t.\tC\tT\t.\t.\t.\tGT\t./.\n')
    out=tmp_path/'run'; m=qc(p,out)
    assert m['counts']['records']==2
    assert m['counts']['multiallelic_records']==1
    assert m['phenotype_read'] is False
    assert m['ranking']=='BLOCKED'
    assert m['counts']['filter_.']==1
    with pytest.raises(ValueError): qc(p,out)

def test_app():
    from fastapi.testclient import TestClient
    from mva_explorer.app import app
    c=TestClient(app)
    assert c.get('/').status_code==200
    assert c.post('/api/qc',json={'filename':'../../etc/passwd'}).status_code==400
    assert c.get('/api/runs/not-a-run').status_code==404

def test_normalization(tmp_path):
    import pysam
    from mva_explorer.core import normalize
    ref=tmp_path/'ref.fa'; ref.write_text('>1\n'+'A'*100+'\n')
    pysam.faidx(str(ref))
    vcf=tmp_path/'in.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10\t.\tAA\tA,AAA\t30\tPASS\t.\n')
    output=tmp_path/'out.vcf.gz'
    normalize(vcf,ref,output)
    with pysam.VariantFile(str(output)) as f:
        records=list(f)
    assert len(records)==2
    assert all(len(r.alts)==1 for r in records)
    assert all(r.pos==1 for r in records)

def test_normalization_rejects_ref_mismatch(tmp_path):
    import pysam
    from mva_explorer.core import normalize
    ref=tmp_path/'ref.fa'; ref.write_text('>1\nAAAA\n'); pysam.faidx(str(ref))
    vcf=tmp_path/'in.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t.\tC\tT\t30\tPASS\t.\n')
    with pytest.raises(Exception): normalize(vcf,ref,tmp_path/'out.vcf.gz')
