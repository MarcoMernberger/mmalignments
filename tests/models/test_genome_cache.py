import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mmalignments.models.genome import _file_specs, _find_existing


class GenomeCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.cwd = Path.cwd()
        os.chdir(self.temporary.name)
        self.addCleanup(os.chdir, self.cwd)
        self.argv = patch.object(sys, "argv", ["run.py"])
        self.argv.start()
        self.addCleanup(self.argv.stop)
        self.species = Path(self.temporary.name) / "ensembl" / "Homo_sapiens"
        self.group = self.species / "dna"
        self.spec = _file_specs("Homo_sapiens")["fasta"]

    def make_file(self, group, digest, name, content):
        path = self.species / group / "by_input" / digest / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def history(self, data, per_script=True):
        directory = (
            Path(".ppg/per_script/run.py/history")
            if per_script
            else Path(".ppg/history")
        )
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "SharedMultiFileGeneratingJobs.json").write_text(json.dumps(data))

    def test_per_script_history_accepts_deduplicated_older_url(self):
        digest, server_digest = "a" * 128, "b" * 128
        fasta = self.make_file("dna", digest, "genome.fasta", ">1\nAC\n")
        Path(str(fasta) + ".url").write_text("https://x/release-110/genome.gz")
        self.make_file("server", server_digest, "url.txt", "https://x/release-116/")
        self.history(
            {str(self.group): digest, str(self.species / "server"): server_digest}
        )
        self.assertEqual(_find_existing(self.group, self.spec, 116), fasta)

    def test_history_for_different_release_is_not_trusted(self):
        digest, server_digest = "a" * 128, "b" * 128
        self.make_file("dna", digest, "genome.fasta", ">1\nAC\n")
        self.make_file("server", server_digest, "url.txt", "https://x/release-115/")
        self.history(
            {str(self.group): digest, str(self.species / "server"): server_digest}
        )
        self.assertIsNone(_find_existing(self.group, self.spec, 116))

    def test_missing_history_finds_matching_marker(self):
        self.make_file("dna", "a" * 128, "genome.fasta", ">1\nAC\n")
        fasta = self.make_file("dna", "b" * 128, "genome.fasta", ">1\nAC\n")
        Path(str(fasta) + ".url").write_text("https://x/release-116/genome.gz")
        self.assertEqual(_find_existing(self.group, self.spec, 116), fasta)
        self.assertIsNone(_find_existing(self.group, self.spec, 115))

    def test_legacy_history_is_used(self):
        digest, server_digest = "a" * 128, "b" * 128
        fasta = self.make_file("dna", digest, "genome.fasta", ">1\nAC\n")
        self.make_file("server", server_digest, "url.txt", "https://x/release-116/")
        self.history(
            {str(self.group): digest, str(self.species / "server"): server_digest},
            per_script=False,
        )
        self.assertEqual(_find_existing(self.group, self.spec, 116), fasta)

    def test_symlinked_build_is_reused(self):
        digest, server_digest = "a" * 128, "b" * 128
        target = self.species / "dna" / "done" / "output"
        target.mkdir(parents=True)
        (target / "genome.fasta").write_text(">1\nAC\n")
        (target / "genome.fasta.url").write_text("https://x/release-110/genome.gz")
        link = self.group / "by_input" / digest
        link.parent.mkdir(parents=True)
        link.symlink_to(target, target_is_directory=True)
        self.make_file("server", server_digest, "url.txt", "https://x/release-116/")
        self.history(
            {str(self.group): digest, str(self.species / "server"): server_digest}
        )
        self.assertEqual(_find_existing(self.group, self.spec, 116), link / "genome.fasta")


if __name__ == "__main__":
    unittest.main()
