"""Disease-profile HPO similarity; only invoked explicitly AFTER G0 freezing."""
from collections import defaultdict
from functools import lru_cache
import csv
import math


def hpo_evidence(phenotype, ontology, annotations, genes, validity=()):
    parents, aliases, obsolete = defaultdict(set), {}, set()
    terms, current = set(), None
    with open(ontology) as f:
        for raw in f:
            line = raw.strip()
            if line.startswith('['):
                current = None
            elif line.startswith('id: HP:'):
                current = line[4:]; terms.add(current)
            elif current:
                if line.startswith('alt_id: '): aliases[line[8:]] = current
                elif line.startswith('is_a: '): parents[current].add(line[6:].split()[0])
                elif line == 'is_obsolete: true': obsolete.add(current)
    def canonical(t):
        t = aliases.get(t, t)
        if t not in terms or t in obsolete:
            raise ValueError('Unknown or obsolete HPO term: ' + t)
        return t
    with open(phenotype) as f:
        patient = sorted({canonical(r['hpo_id']) for r in csv.DictReader(f, delimiter='\t')
                          if r.get('role') == 'proband' and r.get('present', 'true').lower() not in ('false', '0', 'no')})
    if not patient:
        raise ValueError('At least one present proband HPO term is required')
    diseases, gene_diseases = defaultdict(set), defaultdict(set)
    with open(annotations) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            disease, gene = row['disease_id'], row['gene_symbol']
            term = aliases.get(row['hpo_id'], row['hpo_id'])
            if not disease or term not in terms or term in obsolete:
                continue
            diseases[disease].add(term)
            gene_diseases[gene].add(disease)
    @lru_cache(None)
    def ancestors(t):
        # Iterative traversal tolerates cycles in malformed input without recursion.
        seen, pending = set(), [t]
        while pending:
            v = pending.pop()
            if v in seen: continue
            seen.add(v); pending.extend(parents[v])
        return frozenset(seen)
    counts = defaultdict(int)
    for profile in diseases.values():
        for term in set().union(*(ancestors(t) for t in profile)):
            counts[term] += 1
    def ic(t):
        return -math.log((counts[t] + 1) / (len(diseases) + 1))
    @lru_cache(None)
    def similarity(disease):
        # IC-weighted Resnik-style patient->disease similarity.
        # Generic HPO terms have low IC and therefore cannot contribute as much
        # as specific phenotypic findings. This remains a prioritization metric,
        # not ACMG PP4 evidence.
        numerator = denominator = 0.0
        for p in patient:
            pic = ic(p)
            if pic <= 0:
                continue
            best = max((max((ic(a) for a in ancestors(p) & ancestors(d)), default=0)
                        for d in diseases[disease]), default=0)
            numerator += min(pic, best)
            denominator += pic
        return numerator / denominator if denominator else 0.0
    validity_by_gene = defaultdict(list)
    for v in validity:
        if not all(v.get(k) for k in ('gene', 'disease_id', 'inheritance_model', 'source', 'version', 'inheritance_support')):
            raise ValueError('Validity rows require explicit disease, model, support and provenance')
        validity_by_gene[v['gene']].append(v)
    result = []
    for gene, symbol in genes.items():
        candidate_diseases = gene_diseases[symbol] | {v['disease_id'] for v in validity_by_gene[gene]}
        for disease in sorted(candidate_diseases):
            matches = [v for v in validity_by_gene[gene] if v['disease_id'] == disease]
            for validity_row in matches or [{}]:
                result.append({'gene': gene, 'disease_id': disease, 'score': similarity(disease),
                               'source': str(annotations), 'version': 'input_sha256_in_layer_manifest',
                               'dependency_group': 'PHENOTYPE_SIMILARITY',
                               'similarity_method': 'IC_WEIGHTED_RESNIK_PATIENT_TO_DISEASE_NOT_PP4',
                               'phenotype_semantics': 'CASE_FIT_NOT_VARIANT_PATHOGENICITY',
                               'PP4': 'NOT_ASSIGNED',
                               'inheritance_model': validity_row.get('inheritance_model', ''),
                               'inheritance_support': validity_row.get('inheritance_support', 'UNKNOWN'),
                               'validity_source': validity_row.get('source', ''),
                               'validity_version': validity_row.get('version', ''),
                               'validity_support_basis': validity_row.get('support_basis', '')})
    return result
