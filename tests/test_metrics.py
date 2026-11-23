import math
import unittest
import tempfile
from pathlib import Path
import gzip

import numpy as np

from genoflowstats.cli import assembly_metrics, comparison_metrics, fasta_total_length, coverage_hist_from_stats, gc_depth_metrics, nx_profile, robust_scaled, samtools_stats_metrics


class MetricsTests(unittest.TestCase):
    def test_fasta_lengths_and_unreachable_ng(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "assembly.fa"
            p.write_text(">one\nACGTNN\n>two\nGG\n")
            vals, _, _, _ = assembly_metrics(p, genome_size=20)
            self.assertEqual(vals["total_length"], 8)
            self.assertEqual(vals["n50"], 6)
            self.assertIsNone(vals["ng50"])
            self.assertIsNone(vals["lg50"])
            p.write_text("ACGT\n>one\nACGT\n")
            with self.assertRaises(ValueError):
                assembly_metrics(p)
            p.write_text(">one\nAC GT\n")
            with self.assertRaises(ValueError):
                assembly_metrics(p)

    def test_reference_fasta_automatic_and_target_safety(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "ref.fa"
            fasta.write_text(">a\nACGTNN\n>b\nGG\n")
            gz = root / "ref.fa.gz"
            with gzip.open(gz, "wt") as out:
                out.write(fasta.read_text())
            self.assertEqual(fasta_total_length(gz), 8)
            stats = root / "stats.txt"
            cov = "COV\t[1-1]\t1\t6\n"
            stats.write_text(f"# The command line was: stats --coverage 1,1000,1 --ref-seq {fasta} reads.bam\n" + cov)
            vals, _ = samtools_stats_metrics(stats)
            self.assertEqual(vals["coverage_zero_sites"], 2)
            self.assertEqual(vals["coverage_denominator_source"], "reference_fasta")
            manifest = root / "manifest.tsv"
            manifest.write_text(f"name\tfasta\tsamtools_stats\nx\t{fasta}\t{stats}\n")
            rows, _ = comparison_metrics(manifest, None, 10, 1)
            self.assertEqual(rows[0]["stats.coverage_zero_sites"], 2)
            self.assertEqual(rows[0]["stats.coverage_denominator_source"], "manifest_fasta")
            stats.write_text(f"# The command line was: stats -t targets.bed --ref-seq {fasta} reads.bam\nSN\tbases inside the target:\t7\n" + cov)
            vals, _ = samtools_stats_metrics(stats, reference_fasta=fasta)
            self.assertEqual(vals["coverage_zero_sites"], 1)
            self.assertEqual(vals["coverage_denominator_source"], "samtools")
            stats.write_text(f"# The command line was: stats -t targets.bed --ref-seq {fasta} reads.bam\n" + cov)
            vals, _ = samtools_stats_metrics(stats, reference_fasta=fasta)
            self.assertIsNone(vals["coverage_zero_sites"])
            stats.write_text(cov)
            vals, _ = samtools_stats_metrics(stats, reference_fasta=fasta)
            self.assertEqual(vals["coverage_zero_sites"], 2)

    def test_gc_cumulative_percentiles_are_differenced(self):
        rows = [dict(gc_percent=40, unique_sequence_percent=20, depth_p50=10),
                dict(gc_percent=50, unique_sequence_percent=100, depth_p50=20)]
        self.assertAlmostEqual(gc_depth_metrics(rows)["gc_depth_weighted_mean_p50"], 18)

    def test_coverage_denominator_and_bins(self):
        data = {"sn": {}, "rfs_summary": None, "cov_rows": [("[1-1]", 1, 8), ("[2-2]", 2, 1)]}
        hist, meta = coverage_hist_from_stats(data, 10)
        self.assertEqual(hist.tolist(), [1, 8, 1])
        self.assertEqual(meta["coverage_fraction_scope"], "reference_or_target")
        with self.assertRaises(ValueError):
            coverage_hist_from_stats(data, 8)
        data["cov_rows"] = [("[1-5]", 1, 8)]
        with self.assertRaises(ValueError):
            coverage_hist_from_stats(data)

    def test_unknown_reference_is_not_breadth(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "stats"
            p.write_text("COV\t[1-1]\t1\t8\nCOV\t[2-2]\t2\t1\n")
            vals, _ = samtools_stats_metrics(p)
            self.assertIsNone(vals["coverage_zero_fraction"])
            self.assertNotIn("coverage_fraction_ge_1x", vals)
            self.assertEqual(vals["coverage_fraction_ge_1x_covered_sites"], 1)

    def test_missing_values_and_unreachable_ngx(self):
        self.assertTrue(math.isnan(robust_scaled(np.array([1., np.nan]), "higher")[1]))
        profile = nx_profile(np.array([5, 3]), 20)
        self.assertTrue(math.isnan(profile[49]))
        self.assertEqual(profile[19], 5)


if __name__ == "__main__":
    unittest.main()
