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
    'chipapprox_old' = 'test_approx',
    'polyenrich_slow' = 'test_polyenrich_slow',
    'polyenrich_weighted' = 'test_polyenrich_weighted',
    'chipapprox' = 'test_chipapprox',
    'polyapprox' = 'test_polyapprox'
)

METHOD_NAMES = list(
    'broadenrich' = "Broad-Enrich",
    'broadenrich_splineless' = "Broad-Enrich Splineless",
    'chipenrich' = "ChIP-Enrich",
    'chipenrich_slow' = "ChIP-Enrich Original",
	'chipapprox_old' = "ChIP-Enrich Approximate (Old)",
    'chipapprox' = "ChIP-Enrich Approximate",
    'fet' = "Fisher's Exact Test",
    'polyenrich' = "Poly-Enrich",
    'polyenrich_slow' = "Poly-Enrich Original",
    'polyenrich_weighted' = "Poly-Enrich Weighted",
    'polyapprox' = "Poly-Enrich Approximate"

)
