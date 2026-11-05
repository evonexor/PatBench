.PHONY: test verify schema-doc

test:            ## unit tests (no data needed; stage tests run when RDKit is installed)
	python -m pytest -q

verify:          ## code against the recorded hashes
	python -m patbench verify-code

schema-doc:      ## regenerate docs/DATA_SCHEMA.md from the schema dump of the r6 release
	python tools/make_schema_doc.py docs/data/schemas_r6.json > docs/DATA_SCHEMA.md
