"""
Test QPX format compatibility: verify mokume can read both new and legacy QPX parquet formats.

New QPX format columns:
  - charge, run_file_name, intensities[{label, intensity}]

Legacy QPX format columns:
  - precursor_charge, reference_file_name, intensities[{sample_accession, channel, intensity}]
"""


import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


def _make_new_qpx_parquet(path: str) -> None:
    """Create a mock parquet file in new QPX format (matches latest QPX schema).

    Key differences from legacy:
    - pg_accessions: list<struct{accession, start, end, pre, post}> (not list<string>)
    - unique: bool (not int)
    - anchor_protein: string (new field)
    - charge: int16 (not precursor_charge)
    - run_file_name: string (not reference_file_name)
    - intensities: list<struct{label, intensity}> (not {sample_accession, channel, intensity})
    - is_decoy: bool (new field)
    """
    intensities_type = pa.list_(
        pa.struct([("label", pa.string()), ("intensity", pa.float32())])
    )
    pg_protein_type = pa.list_(
        pa.struct([
            ("accession", pa.string()),
            ("start", pa.int32()),
            ("end", pa.int32()),
            ("pre", pa.string()),
            ("post", pa.string()),
        ])
    )
    schema = pa.schema(
        [
            ("sequence", pa.string()),
            ("peptidoform", pa.string()),
            ("pg_accessions", pg_protein_type),
            ("anchor_protein", pa.string()),
            ("charge", pa.int16()),
            ("run_file_name", pa.string()),
            ("unique", pa.bool_()),
            ("is_decoy", pa.bool_()),
            ("intensities", intensities_type),
        ]
    )

    data = {
        "sequence": ["PEPTIDEK", "ANOTHERPEPTIDE", "PEPTIDEK"],
        "peptidoform": ["PEPTIDEK", "ANOTHERPEPTIDE", "PEPTIDEK"],
        "pg_accessions": [
            [{"accession": "sp|P12345|PROT_HUMAN", "start": 10, "end": 18, "pre": "K", "post": "A"}],
            [{"accession": "sp|P67890|PROT2_HUMAN", "start": 5, "end": 19, "pre": "R", "post": "L"}],
            [{"accession": "sp|P12345|PROT_HUMAN", "start": 10, "end": 18, "pre": "K", "post": "A"}],
        ],
        "anchor_protein": ["sp|P12345|PROT_HUMAN", "sp|P67890|PROT2_HUMAN", "sp|P12345|PROT_HUMAN"],
        "charge": [2, 3, 2],
        "run_file_name": ["run1", "run1", "run2"],
        "unique": [True, True, True],
        "is_decoy": [False, False, False],
        "intensities": [
            [{"label": "TMT126", "intensity": 1000.0}, {"label": "TMT127", "intensity": 2000.0}],
            [{"label": "TMT126", "intensity": 3000.0}, {"label": "TMT127", "intensity": 4000.0}],
            [{"label": "TMT126", "intensity": 1500.0}, {"label": "TMT127", "intensity": 2500.0}],
        ],
    }
    table = pa.table(data, schema=schema)
    pq.write_table(table, path)


def _make_legacy_qpx_parquet(path: str) -> None:
    """Create a mock parquet file in legacy QPX format."""
    intensities_type = pa.list_(
        pa.struct([
            ("sample_accession", pa.string()),
            ("channel", pa.string()),
            ("intensity", pa.float64()),
        ])
    )
    schema = pa.schema(
        [
            ("sequence", pa.string()),
            ("peptidoform", pa.string()),
            ("pg_accessions", pa.list_(pa.string())),
            ("precursor_charge", pa.int32()),
            ("reference_file_name", pa.string()),
            ("unique", pa.int32()),
            ("intensities", intensities_type),
        ]
    )

    data = {
        "sequence": ["PEPTIDEK", "ANOTHERPEPTIDE", "PEPTIDEK"],
        "peptidoform": ["PEPTIDEK", "ANOTHERPEPTIDE", "PEPTIDEK"],
        "pg_accessions": [["P12345"], ["P67890"], ["P12345"]],
        "precursor_charge": [2, 3, 2],
        "reference_file_name": ["run1.raw", "run1.raw", "run2.raw"],
        "unique": [1, 1, 1],
        "intensities": [
            [
                {"sample_accession": "S1", "channel": "TMT126", "intensity": 1000.0},
                {"sample_accession": "S2", "channel": "TMT127", "intensity": 2000.0},
            ],
            [
                {"sample_accession": "S1", "channel": "TMT126", "intensity": 3000.0},
                {"sample_accession": "S2", "channel": "TMT127", "intensity": 4000.0},
            ],
            [
                {"sample_accession": "S1", "channel": "TMT126", "intensity": 1500.0},
                {"sample_accession": "S2", "channel": "TMT127", "intensity": 2500.0},
            ],
        ],
    }
    table = pa.table(data, schema=schema)
    pq.write_table(table, path)


class TestNewQPXFormat:
    """Test mokume Feature reader with new QPX format."""

    def test_feature_init_new_format(self, tmp_path):
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        assert feat._is_new_qpx is True
        assert feat._charge_col == "charge"
        assert feat._run_col == "run_file_name"

    def test_feature_query_new_format(self, tmp_path):
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        df = feat.parquet_db.execute("SELECT * FROM parquet_db").df()
        assert len(df) > 0
        assert "charge" in df.columns
        assert "run_file_name" in df.columns
        assert "intensity" in df.columns
        assert "channel" in df.columns
        assert "sample_accession" in df.columns

    def test_feature_samples_new_format(self, tmp_path):
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        samples = feat.get_unique_samples()
        assert len(samples) > 0


class TestLegacyQPXFormat:
    """Test mokume Feature reader with legacy QPX format."""

    def test_feature_init_legacy_format(self, tmp_path):
        parquet_file = str(tmp_path / "legacy_qpx.feature.parquet")
        _make_legacy_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        assert feat._is_new_qpx is False
        assert feat._charge_col == "precursor_charge"
        assert feat._run_col == "reference_file_name"

    def test_feature_query_legacy_format(self, tmp_path):
        parquet_file = str(tmp_path / "legacy_qpx.feature.parquet")
        _make_legacy_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        df = feat.parquet_db.execute("SELECT * FROM parquet_db").df()
        assert len(df) > 0
        assert "charge" in df.columns
        assert "run_file_name" in df.columns
        assert "intensity" in df.columns
        assert "channel" in df.columns
        assert "sample_accession" in df.columns

    def test_feature_samples_legacy_format(self, tmp_path):
        parquet_file = str(tmp_path / "legacy_qpx.feature.parquet")
        _make_legacy_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        samples = feat.get_unique_samples()
        assert len(samples) > 0


class TestNewQPXDeepCompat:
    """Test deep compatibility: pg_accessions struct parsing, unique bool, etc."""

    def test_pg_accessions_struct_in_pandas(self, tmp_path):
        """Verify pg_accessions list<struct> can be parsed like list<string>."""
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        df = feat.parquet_db.execute("SELECT pg_accessions FROM parquet_db").df()
        # In new QPX, pg_accessions is list<struct{accession,...}>
        # mokume code does: df["pg_accessions"].str[0] then .split("|")
        first_elem = df["pg_accessions"].str[0]
        # With struct, first_elem would be a dict, not a string
        print(f"pg_accessions type: {type(first_elem.iloc[0])}")
        print(f"pg_accessions[0]: {first_elem.iloc[0]}")

        # This is what mokume ratio.py does:
        try:
            first_acc = df["pg_accessions"].str[0].fillna("")
            result = np.where(
                first_acc.str.contains("|", regex=False),
                first_acc.str.split("|").str[1],
                first_acc,
            )
            print(f"Parsed protein names: {result}")
            parsed_ok = True
        except Exception as e:
            print(f"FAILED to parse pg_accessions: {e}")
            parsed_ok = False

        assert parsed_ok, "pg_accessions struct parsing failed - needs compatibility fix"

    def test_unique_bool_filter_sql(self, tmp_path):
        """Verify 'unique = 1' SQL filter works with bool column."""
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        # SQLFilterBuilder generates: "unique" = 1
        df = feat.parquet_db.execute(
            'SELECT * FROM parquet_db WHERE "unique" = 1'
        ).df()
        assert len(df) > 0, "unique=1 filter on bool column returned no rows"

    def test_unique_bool_filter_pandas(self, tmp_path):
        """Verify unique == 1 Pandas filter works with bool column."""
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        df = feat.parquet_db.execute("SELECT * FROM parquet_db").df()
        # stages.py and peptide.py do: dataset_df[dataset_df["unique"] == 1]
        filtered = df[df["unique"] == 1]
        assert len(filtered) > 0, "unique==1 filter on bool column returned no rows in Pandas"

    def test_get_low_frequency_peptides(self, tmp_path):
        """Test get_low_frequency_peptides with struct pg_accessions."""
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature
        feat = Feature(parquet_file)

        try:
            result = feat.get_low_frequency_peptides(percentage=0.2)
            print(f"Low frequency peptides: {result}")
            lfp_ok = True
        except Exception as e:
            print(f"FAILED get_low_frequency_peptides: {e}")
            lfp_ok = False

        assert lfp_ok, "get_low_frequency_peptides failed with struct pg_accessions"

    def test_contaminant_filter_with_struct(self, tmp_path):
        """Test SQL contaminant filter (pg_accessions::text LIKE) with struct type."""
        parquet_file = str(tmp_path / "new_qpx.feature.parquet")
        _make_new_qpx_parquet(parquet_file)

        from mokume.io.feature import Feature, SQLFilterBuilder
        fb = SQLFilterBuilder(remove_contaminants=True)
        feat = Feature(parquet_file, filter_builder=fb)

        where_clause = fb.build_where_clause()
        df = feat.parquet_db.execute(
            f"SELECT * FROM parquet_db WHERE {where_clause}"
        ).df()
        # Should still return rows since our test data has no contaminants
        assert len(df) > 0, "Contaminant filter on struct pg_accessions returned no rows"


class TestBothFormatsProduceSameSchema:
    """Verify that both formats produce compatible output schemas."""

    def test_same_core_columns(self, tmp_path):
        """Both formats must produce all core columns needed by downstream code."""
        new_file = str(tmp_path / "new.parquet")
        legacy_file = str(tmp_path / "legacy.parquet")
        _make_new_qpx_parquet(new_file)
        _make_legacy_qpx_parquet(legacy_file)

        from mokume.io.feature import Feature
        feat_new = Feature(new_file)
        feat_legacy = Feature(legacy_file)

        df_new = feat_new.parquet_db.execute("SELECT * FROM parquet_db").df()
        df_legacy = feat_legacy.parquet_db.execute("SELECT * FROM parquet_db").df()

        core_columns = {
            "sequence", "peptidoform", "pg_accessions", "charge",
            "run_file_name", "unique", "sample_accession", "channel",
            "intensity", "run", "condition", "biological_replicate",
            "fraction", "mixture",
        }
        assert core_columns.issubset(set(df_new.columns)), (
            f"New format missing core columns: {core_columns - set(df_new.columns)}"
        )
        assert core_columns.issubset(set(df_legacy.columns)), (
            f"Legacy format missing core columns: {core_columns - set(df_legacy.columns)}"
        )

    def test_new_format_has_extra_columns(self, tmp_path):
        """New QPX format should expose is_decoy and anchor_protein."""
        new_file = str(tmp_path / "new.parquet")
        _make_new_qpx_parquet(new_file)

        from mokume.io.feature import Feature
        feat = Feature(new_file)
        df = feat.parquet_db.execute("SELECT * FROM parquet_db").df()

        assert "is_decoy" in df.columns, "New QPX should expose is_decoy"
        assert "anchor_protein" in df.columns, "New QPX should expose anchor_protein"
