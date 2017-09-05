SUPPORTED_METHODS = list(
    'broadenrich' = "test_broadenrich",
    'chipenrich' = "test_chipenrich",
    'fet' = "test_fisher_exact",
    'polyenrich' = "test_polyenrich"
)

HIDDEN_METHODS = list(
    'binomial' = "test_binomial",
    'broadenrich_splineless' = "test_broadenrich_splineless",
    'chipenrich_slow' = 'test_chipenrich_slow',
    'chipapprox' = 'test_approx',
    'polyenrich_slow' = 'test_polyenrich_slow',
    'polyenrich_weighted' = 'test_polyenrich_weighted'
    'chipenrich_score' = 'test_chipenrich_score',
    'polyenrich_score' = 'test_polyenrich_score'
)

METHOD_NAMES = list(
    'broadenrich' = "Broad-Enrich",
    'broadenrich_splineless' = "Broad-Enrich Splineless",
    'chipenrich' = "ChIP-Enrich",
    'chipenrich_slow' = "ChIP-Enrich Original",
	'chipapprox' = "ChIP-Enrich Approximate",
    'fet' = "Fisher's Exact Test",
    'polyenrich' = "Poly-Enrich",
    'polyenrich_slow' = "Poly-Enrich Original",
    'polyenrich_weighted' = "Poly-Enrich Weighted"
)
