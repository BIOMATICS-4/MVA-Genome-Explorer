from mva_explorer.acmg_preview import combine_2015


def c(code, direction='pathogenic', strength='supporting'):
    return {'code':code,'direction':direction,'strength':strength,'status':'MET'}


def test_pvs1_plus_strong_is_pathogenic():
    r=combine_2015([c('PVS1',strength='very_strong'),c('PS3',strength='strong')])
    assert r['classification'].startswith('Pathogenic')


def test_pvs1_plus_moderate_is_likely_pathogenic():
    r=combine_2015([c('PVS1',strength='very_strong'),c('PM3',strength='moderate')])
    assert r['classification'].startswith('Likely pathogenic')


def test_three_moderate_is_likely_pathogenic():
    r=combine_2015([c('PM1',strength='moderate'),c('PM2',strength='moderate'),c('PM3',strength='moderate')])
    assert r['classification'].startswith('Likely pathogenic')


def test_ba1_is_benign():
    r=combine_2015([c('BA1','benign','stand_alone')])
    assert r['classification'].startswith('Benign')


def test_bs_plus_bp_is_likely_benign():
    r=combine_2015([c('BS1','benign','strong'),c('BP4','benign','supporting')])
    assert r['classification'].startswith('Likely benign')


def test_no_met_criteria_is_uncertain():
    r=combine_2015([{'code':'PVS1','direction':'pathogenic','strength':'very_strong','status':'REVIEW'}])
    assert r['classification'].startswith('Uncertain significance')
