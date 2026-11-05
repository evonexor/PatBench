"""Evidence gate for proposed repairs; unlisted, unresolved, or altered records are rejected.

This module does not change a structure or infer a nitrile from an isocyanide pattern.
Call require_approved(patent_core, structure_file, original_smiles, proposed_smiles)
before accepting a proposed repair. Pin record_decisions.parquet to this audit.
"""
import hashlib
from pathlib import Path
import pandas as pd
from rdkit import Chem
class RepairNotApproved(ValueError):
    pass
class NitrileEvidenceGuard:
    def __init__(self, decisions_path=None):
        p=Path(decisions_path) if decisions_path else Path(__file__).with_name('record_decisions.parquet')
        self.records=pd.read_parquet(p).set_index(['patent_core','structure_file'],verify_integrity=True)
    def require_approved(self,patent_core,structure_file,original_smiles,proposed_smiles):
        key=(patent_core,structure_file)
        if key not in self.records.index:
            raise RepairNotApproved('No source evidence for this record')
        r=self.records.loc[key]
        if r.decision!='supported_nitrile':
            raise RepairNotApproved('Source evidence is unresolved: '+str(r.review_note))
        if hashlib.sha256(original_smiles.encode()).hexdigest()!=r.before_smiles_sha256:
            raise RepairNotApproved('Original SMILES differs from the audited input')
        m=Chem.MolFromSmiles(proposed_smiles)
        if m is None or Chem.MolToInchiKey(m)!=r.expected_after_inchikey:
            raise RepairNotApproved('Proposed structure differs from the audited candidate')
        return True
