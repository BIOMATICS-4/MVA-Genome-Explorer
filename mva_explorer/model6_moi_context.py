"""Provenance-preserving GenCC MOI audit for Model 6.

This is an evidence annotation, NOT a replacement for phenotype-matched
condition-specific diagnosis or a calibrated ranking prior.
"""
from collections import defaultdict
import csv
from pathlib import Path

RELIABLE = {'Definitive','Strong','Moderate'}
SUPPORTED_MODELS = {'Autosomal dominant':'AD', 'Autosomal recessive':'AR', 'X-linked':'XL', 'X-linked recessive':'XLR'}

def load_gencc(path):
    grouped=defaultdict(list)
    with Path(path).open(newline='') as h:
        for r in csv.DictReader(h):
            if r.get('gene_symbol') and r.get('moi_title') in SUPPORTED_MODELS:
                grouped[r['gene_symbol']].append({
                    'disease_id':r.get('disease_curie',''),
                    'disease_original_id':r.get('disease_original_curie',''),
                    'disease_submitted_id':r.get('submitted_as_disease_id',''),
                    'disease_name':r.get('disease_title',''),
                    'moi':SUPPORTED_MODELS[r['moi_title']],
                    'classification':r.get('classification_title',''),
                    'source':r.get('submitter_title',''),
                    'submission_id':r.get('sgc_id',''),
                })
    return grouped


def audit_gene_moi(gene, ranked_moi, grouped):
    rows=grouped.get(gene, [])
    strong={x['moi'] for x in rows if x['classification'] in RELIABLE}
    weak={x['moi'] for x in rows if x['classification'] not in RELIABLE and x['classification'] in {'Supportive','Limited'}}
    if ranked_moi in strong:classification='SOME_DISEASES_HAVE_MODERATE_OR_STRONG_GENE_MOI_EVIDENCE'
    elif ranked_moi in weak:classification='ONLY_LIMITED_OR_SUPPORTIVE_GENE_MOI_EVIDENCE'
    elif strong:classification='NO_MODERATE_OR_STRONG_EVIDENCE_FOR_THIS_MOI_BUT_OTHER_MOI_SUPPORTED'
    else:classification='GENCC_MOI_EVIDENCE_NOT_ESTABLISHED'
    # Different diseases in a gene can have different inheritance; do not
    # collapse them into a single truth about the patient.
    return {'moi_evidence_status':classification,
            'known_supported_gene_moi':';'.join(sorted(str(x) for x in strong)),
            'reported_diseases_for_ranked_moi':';'.join(sorted({x['disease_id'] for x in rows if x['moi']==ranked_moi and x['classification'] in RELIABLE})),
            'moi_evidence_caveat':'Not matched to patient disease phenotype; no automatic rank change'}
