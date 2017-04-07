SUPPORTED_METHODS = list(
    'chipenrich' = "test_gam_fast",
    'fet' = "test_fisher_exact",
    'broadenrich' = "test_gam_ratio",
    'polyenrich' = "test_gam_nb_fast"
)

HIDDEN_METHODS = list(
    'binomial' = "test_binomial",
    'broadenrich_splineless' = "test_gam_ratio_splineless",
    'chipenrich_slow' = 'test_gam',
    'polyenrich_slow' = 'test_gam_nb',
    'chipapprox' = 'test_approx'
)

METHOD_NAMES = list(
	'fet' = "Fisher's Exact Test",
	'chipenrich' = "ChIP-Enrich",
	'broadenrich' = "Broad-Enrich",
	'chipapprox' = "ChIP-Enrich Approximate",
	'broadenrich_splineless' = "Broad-Enrich Splineless",
	'polyenrich' = "Poly-Enrich",
	'chipenrich_fast' = "ChIP-Enrich Fast"
)
